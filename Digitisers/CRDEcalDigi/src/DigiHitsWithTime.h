#ifndef _DIGIHITSWITHTIME_
#define _DIGIHITSWITHTIME_

#include "CRDEcalDigiAlg.h"
#include  "CRDEcalDigiEDM.h"

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

std::vector<edm4hep::ConstCalorimeterHit> CRDEcalDigiAlg::DigiHitsWithTime(std::vector<CRDEcalDigiEDM::DigiBar>& m_block){
   int Nbars = m_block.size();
   if(Nbars==0) exit(0);

   std::vector<edm4hep::ConstCalorimeterHit> m_digiCol; m_digiCol.clear();
	//std::vector<edm4hep::CalorimeterHit> m_digiColForClus; m_digiColForClus.clear();
   float rotAngle = m_block[0].module*PI/4.;

   TVector3 m_vec(0,0,0);
   dd4hep::Position m_pos(0,0,0);
   for(int ibar=0;ibar<Nbars;ibar++){
      CRDEcalDigiEDM::DigiBar bar_s0 = m_block[ibar];

      m_vec.SetXYZ(bar_s0.position.x(), bar_s0.position.y(), bar_s0.position.z());
      m_vec.RotateZ(rotAngle);
      bar_s0.position.SetXYZ(m_vec.x(), m_vec.y(), m_vec.z());

      double dis_s0 = C*(bar_s0.T1-bar_s0.T2)/(2*nMat);

      TVector3 p_hit(0,0,0);
		if(bar_s0.slayer==0) p_hit.SetXYZ(bar_s0.position.x()+dis_s0, bar_s0.position.y(), bar_s0.position.z() );
		if(bar_s0.slayer==1) p_hit.SetXYZ(bar_s0.position.x(), bar_s0.position.y(), bar_s0.position.z()+dis_s0 );
      p_hit.RotateZ(-rotAngle);
      edm4hep::Vector3f m_vec3f(p_hit.x(), p_hit.y(), p_hit.z());
		float m_En = r_cali*(m_block[ibar].Q1+m_block[ibar].Q2)/2;

      edm4hep::CalorimeterHit hit;
      hit.setCellID(0);
      hit.setPosition(m_vec3f);
      hit.setEnergy(m_En);
      m_digiCol.push_back(hit);
   }

   return m_digiCol;
}

#endif
