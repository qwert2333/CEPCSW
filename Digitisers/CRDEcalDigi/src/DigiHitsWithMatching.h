#ifndef _DIGIHITSWITHMATCHING_
#define _DIGIHITSWITHMATCHING_

#include "CRDEcalDigiAlg.h"
#include "DigiHitsWithPos.h"
#include "DigiHitsWithTime.h"

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

std::vector<edm4hep::ConstCalorimeterHit> CRDEcalDigiAlg::DigiHitsWithMatching(std::vector<CRDEcalDigiEDM::BarCollection>& barShowerXCol, std::vector<CRDEcalDigiEDM::BarCollection>& barShowerYCol){

	std::vector<edm4hep::ConstCalorimeterHit> m_digiCol; m_digiCol.clear();

	const int Nshower = barShowerXCol.size();
	double chi2[Nshower][Nshower];
	double chi2_E[Nshower][Nshower];
	double chi2_tx[Nshower][Nshower];
	double chi2_ty[Nshower][Nshower];

	double sigmaE = 0.1;
	double sigmaPos = sqrt(10*10 + pow((Tres*C/(2*nMat)),2) );	//const
	double wi_E = _chi2Wi_E/(_chi2Wi_E + _chi2Wi_T);
	double wi_T = _chi2Wi_T/(_chi2Wi_E + _chi2Wi_T);

	TVector3 m_vec(0,0,0);
	double rotAngle = (barShowerXCol[0].Bars)[0].module*PI/4.;

	for(int ix=0;ix<Nshower;ix++){
	for(int iy=0;iy<Nshower;iy++){
		CRDEcalDigiEDM::BarCollection showerX = barShowerXCol[ix];
		CRDEcalDigiEDM::BarCollection showerY = barShowerYCol[iy];

		double Ex = showerX.getE();
		double Ey = showerY.getE();
		chi2_E[ix][iy] = pow(fabs(Ex-Ey)/sigmaE, 2);

		double PosTx = C*(showerX.getT1()-showerX.getT2())/(2*nMat) + showerX.getPos().z();
		double PosTy = C*(showerY.getT1()-showerY.getT2())/(2*nMat);
		chi2_tx[ix][iy] = pow( fabs(PosTx-showerY.getPos().z())/sigmaPos, 2);

		m_vec.SetXYZ(showerX.getPos().x(), showerX.getPos().y(), showerX.getPos().z());
		m_vec.RotateZ(rotAngle);
		chi2_ty[ix][iy] = pow( fabs(PosTy - m_vec.x() )/sigmaPos, 2);

		chi2[ix][iy] = chi2_E[ix][iy]*wi_E + (chi2_tx[ix][iy]+chi2_ty[ix][iy])*wi_T ;
	}}


	int Ncomb=1;
	for(int i=Nshower; i>0; i--) Ncomb = Ncomb*i;

	map<double, vector<pair<int, int>> > matchingMap;
	int num[Nshower];
	int num_init[Nshower]; 
	for(int i=0;i<Nshower;i++){ num[i]=i; num_init[i]=i;}

	for(int icont=0;icont<Ncomb;icont++){
		vector<pair<int, int>> Index;
		for(int i=0;i<Nshower;i++){
         pair<int, int> p1(num_init[i], num[i]);
         Index.push_back(p1);
      }
		double chi2_tot=0;
		for(int i=0;i<Index.size();i++) chi2_tot += chi2[Index[i].first][Index[i].second];
		matchingMap[chi2_tot] = Index;

		Index.clear();
		if(!next_permutation(num, num+Nshower)) break;			
	}

	map<double, vector<pair<int, int>> >::iterator iter = matchingMap.begin();
	vector<pair<int, int>> Index = iter->second;

	for(int i=0;i<Index.size();i++){
      CRDEcalDigiEDM::BarCollection showerX = barShowerXCol[Index[i].first];
      CRDEcalDigiEDM::BarCollection showerY = barShowerYCol[Index[i].second];

		std::vector<edm4hep::ConstCalorimeterHit> m_hitsInShower; m_hitsInShower.clear();
		std::vector<CRDEcalDigiEDM::DigiBar> m_ShowerBlock; m_ShowerBlock.clear();
		m_ShowerBlock.insert(m_ShowerBlock.end(), showerX.Bars.begin(), showerX.Bars.end());
		m_ShowerBlock.insert(m_ShowerBlock.end(), showerY.Bars.begin(), showerY.Bars.end());
		m_hitsInShower = DigiHitsWithPos(m_ShowerBlock);
		m_digiCol.insert(m_digiCol.end(), m_hitsInShower.begin(), m_hitsInShower.end());
	}

	return m_digiCol;

}

#endif
