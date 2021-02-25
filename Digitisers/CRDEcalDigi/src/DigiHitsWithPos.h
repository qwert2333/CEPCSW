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

std::vector<edm4hep::ConstCalorimeterHit> CRDEcalDigiAlg::DigiHitsWithPos( CRDEcalDigiEDM::BarCollection& barShowerX, CRDEcalDigiEDM::BarCollection& barShowerY){

   std::vector<edm4hep::ConstCalorimeterHit> m_digiCol; m_digiCol.clear();
   int NbarsX = barShowerX.Bars.size();
   int NbarsY = barShowerY.Bars.size();
	if(NbarsX==0 || NbarsY==0){ std::cout<<"WARNING: empty DigiHitsCol returned!"<<std::endl; return m_digiCol;}

	float rotAngle = barShowerX.Bars[0].module*PI/4.;

	TVector3 m_vec(0,0,0);
	for(int ibar=0;ibar<NbarsX;ibar++){
		CRDEcalDigiEDM::DigiBar barx = barShowerX.Bars[ibar];

		m_vec.SetXYZ(barx.position.x(), barx.position.y(), barx.position.z());
      m_vec.RotateZ(rotAngle);
		barx.position.SetXYZ(m_vec.x(), m_vec.y(), m_vec.z());

		for(int jbar=0;jbar<NbarsY;jbar++){
			CRDEcalDigiEDM::DigiBar bary = barShowerY.Bars[jbar];
			m_vec.SetXYZ(bary.position.x(), bary.position.y(), bary.position.z());
     		m_vec.RotateZ(rotAngle);
			bary.position.SetXYZ(m_vec.x(), m_vec.y(), m_vec.z());

         TVector3 p_hit(bary.position.x(), (barx.position.y()+bary.position.y())/2., barx.position.z() );
         p_hit.RotateZ(-rotAngle);
         edm4hep::Vector3f m_vec3f(p_hit.x(), p_hit.y(), p_hit.z());
         float m_En = barx.getEnergy()*bary.getEnergy()/barShowerY.getE() + barx.getEnergy()*bary.getEnergy()/barShowerX.getE();
         edm4hep::CalorimeterHit hit;
         hit.setCellID(0);
         hit.setPosition(m_vec3f);
         hit.setEnergy(m_En);
         m_digiCol.push_back(hit);
		}
	}

	return m_digiCol;
}

#endif
