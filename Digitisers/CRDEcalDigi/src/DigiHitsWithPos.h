#ifndef _DIGIHITSWITHPOS_
#define _DIGIHITSWITHPOS_

#include "CRDEcalDigiAlg.h"

#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Vector3f.h"

#include "TRandom3.h"
#include "TVector3.h"
#include <math.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <map>

std::vector<edm4hep::ConstCalorimeterHit> CRDEcalDigiAlg::DigiHitsWithPos(std::vector<CRDEcalDigiEDM::DigiBar>& m_block){
   std::vector<edm4hep::ConstCalorimeterHit> m_digiCol; m_digiCol.clear();

   int Nbars = m_block.size();
   if(Nbars==0){ std::cout<<"WARNING: empty DigiHitsCol returned!"<<std::endl; return m_digiCol;}

	//std::vector<edm4hep::CalorimeterHit> m_digiColForClus; m_digiColForClus.clear();
   float rotAngle = -m_block[0].module*PI/4.;
   float Edes[2][50]={0};           //Deposited energy in slayer 0/1. 
   float totE[2]={0};
   for(int ibar=0;ibar<Nbars;ibar++){
      float m_En = r_cali*(m_block[ibar].Q1+m_block[ibar].Q2)/2;
      Edes[m_block[ibar].slayer][m_block[ibar].bar] = m_En;
      totE[m_block[ibar].slayer] += m_En;
   }

   TVector3 m_vec(0,0,0);
   for(int ibar=0;ibar<Nbars;ibar++){
      if(m_block[ibar].slayer==1) continue;     //Outer loop in slayer==0.
      CRDEcalDigiEDM::DigiBar bar_s0 = m_block[ibar];

      m_vec.SetXYZ(bar_s0.position.x(), bar_s0.position.y(), bar_s0.position.z());
      m_vec.RotateZ(rotAngle);
      bar_s0.position.SetXYZ(m_vec.x(), m_vec.y(), m_vec.z());

      for(int jbar=0;jbar<Nbars;jbar++){
         if(m_block[jbar].slayer==0) continue;  //Inner loop in slayer==1.
         CRDEcalDigiEDM::DigiBar bar_s1 = m_block[jbar];
         m_vec.SetXYZ(bar_s1.position.x(), bar_s1.position.y(), bar_s1.position.z());
         m_vec.RotateZ(rotAngle);
         bar_s1.position.SetXYZ(m_vec.x(), m_vec.y(), m_vec.z());

         TVector3 p_hit(bar_s1.position.x(), (bar_s0.position.y()+bar_s1.position.y())/2., bar_s0.position.z() );
         p_hit.RotateZ(-rotAngle);
         edm4hep::Vector3f m_vec3f(p_hit.x(), p_hit.y(), p_hit.z());
         float m_En = Edes[bar_s0.slayer][bar_s0.bar]*Edes[bar_s1.slayer][bar_s1.bar]/totE[bar_s1.slayer]
                  + Edes[bar_s1.slayer][bar_s1.bar]*Edes[bar_s0.slayer][bar_s0.bar]/totE[bar_s0.slayer];
         edm4hep::CalorimeterHit hit;
         hit.setCellID(0);
         hit.setPosition(m_vec3f);
         hit.setEnergy(m_En);
         m_digiCol.push_back(hit);
			//if(m_En>0.05*(totE[0]+totE[1])) m_digiColForClus.push_back(hit);
      }
   }

   return m_digiCol;
}

#endif
