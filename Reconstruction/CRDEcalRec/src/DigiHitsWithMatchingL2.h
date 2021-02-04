#ifndef _DIGIHITSWITHMATCHINGL2_
#define _DIGIHITSWITHMATCHINGL2_

#include "CRDEcalDigiAlg.h"
#include "DigiHitsWithPos.h"

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

using namespace std;

std::vector<CRDEcalDigiEDM::CRD2DShowerInLayer> CRDEcalDigiAlg::DigiHitsWithMatchingL2(std::vector<CRDEcalDigiEDM::BarCollection>& barShowerXCol, std::vector<CRDEcalDigiEDM::BarCollection>& barShowerYCol){

	std::vector<CRDEcalDigiEDM::CRD2DShowerInLayer> m_showerCol; m_showerCol.clear();

	const int NshowerX = barShowerXCol.size();
	const int NshowerY = barShowerYCol.size();

	double chi2[NshowerX][NshowerY];
	double chi2_E[NshowerX][NshowerY];
	double chi2_tx[NshowerX][NshowerY];
	double chi2_ty[NshowerX][NshowerY];

   double sigmaE = 0.05;  //Energy resolution 5%
   double sigmaPos = sqrt(10*10/12 + pow((Tres*C/(2*nMat)),2) );  //position resolution
   double wi_E = _chi2Wi_E/(_chi2Wi_E + _chi2Wi_T);
   double wi_T = _chi2Wi_T/(_chi2Wi_E + _chi2Wi_T);

   TVector3 m_vec(0,0,0);
   double rotAngle = -(barShowerXCol[0].Bars)[0].module*PI/4.;
   TVector3 Cblock((barShowerXCol[0].Bars)[0].position.x(), (barShowerXCol[0].Bars)[0].position.y(), (barShowerYCol[0].Bars)[0].position.z());
   Cblock.RotateZ(rotAngle);


	for(int ix=0;ix<NshowerX;ix++){
	for(int iy=0;iy<NshowerY;iy++){
		CRDEcalDigiEDM::BarCollection showerX = barShowerXCol[ix];
		CRDEcalDigiEDM::BarCollection showerY = barShowerYCol[iy];

		double Ex = showerX.getE();
		double Ey = showerY.getE();
		chi2_E[ix][iy] = pow(fabs(Ex-Ey)/sigmaE, 2);

      double PosTx = C*(showerY.getT1()-showerY.getT2())/(2*nMat) + showerY.getPos().z();
      chi2_tx[ix][iy] = pow( fabs(PosTx-showerX.getPos().z())/sigmaPos, 2);

      double PosTy = C*(showerX.getT1()-showerX.getT2())/(2*nMat);
      m_vec.SetXYZ(showerY.getPos().x(), showerY.getPos().y(), showerY.getPos().z());
      m_vec.RotateZ(rotAngle);
      chi2_ty[ix][iy] = pow( fabs(PosTy - (m_vec-Cblock).x() )/sigmaPos, 2);

		chi2[ix][iy] = chi2_E[ix][iy]*wi_E + (chi2_tx[ix][iy]+chi2_ty[ix][iy])*wi_T ;

      m_chi2.push_back(chi2[ix][iy]);
      m_chi2E.push_back(chi2_E[ix][iy]);
      m_chi2Tx.push_back(chi2_tx[ix][iy]);
      m_chi2Ty.push_back(chi2_ty[ix][iy]);
	}}


	for(int ix=0;ix<NshowerX;ix++){
	for(int iy=0;iy<NshowerY;iy++){
		if(chi2[ix][iy]<_th_chi2){
		CRDEcalDigiEDM::BarCollection showerX = barShowerXCol[ix];
		CRDEcalDigiEDM::BarCollection showerY = barShowerYCol[iy];
		if(showerX.Bars.size()==0 || showerY.Bars.size()==0) continue;

		CRDEcalDigiEDM::CRD2DShowerInLayer tmp_shower; tmp_shower.Clear();
      tmp_shower.barShowerX = showerX;
      tmp_shower.barShowerY = showerY;
      tmp_shower.CaloHits = DigiHitsWithPos(showerX, showerY);
      m_showerCol.push_back(tmp_shower);

		CRDEcalDigiEDM::BarCollection m_empty; m_empty.Clear();
		barShowerXCol[ix]=m_empty;
		barShowerYCol[iy]=m_empty;
		}
	}}

	for(int is=0;is<barShowerXCol.size();is++){
	for(int js=0;js<barShowerYCol.size();js++){
		if(barShowerXCol[is].Bars.size()==0 || barShowerYCol[js].Bars.size()==0) continue;
		CRDEcalDigiEDM::CRD2DShowerInLayer tmp_shower; tmp_shower.Clear();
		tmp_shower.barShowerX = barShowerXCol[is];
		tmp_shower.barShowerY = barShowerYCol[js];
		tmp_shower.CaloHits = DigiHitsWithPos(tmp_shower.barShowerX, tmp_shower.barShowerY);
		m_showerCol.push_back(tmp_shower);
	}}


	return m_showerCol;

}

#endif
