#ifndef _DIGIHITSWITHMATCHING_
#define _DIGIHITSWITHMATCHING_

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

std::vector<edm4hep::ConstCalorimeterHit> CRDEcalDigiAlg::DigiHitsWithMatching(std::vector<CRDEcalDigiEDM::BarCollection>& barShowerXCol, std::vector<CRDEcalDigiEDM::BarCollection>& barShowerYCol){

	std::vector<edm4hep::ConstCalorimeterHit> m_digiCol; m_digiCol.clear();

	const int Nshower = barShowerXCol.size();
	double chi2[Nshower][Nshower];
	double chi2_E[Nshower][Nshower];
	double chi2_tx[Nshower][Nshower];
	double chi2_ty[Nshower][Nshower];

	double sigmaE = 0.05;  //Energy resolution 5%
	double sigmaPos = sqrt(10*10/12 + pow((Tres*C/(2*nMat)),2) );	//position resolution
	double wi_E = _chi2Wi_E/(_chi2Wi_E + _chi2Wi_T);
	double wi_T = _chi2Wi_T/(_chi2Wi_E + _chi2Wi_T);

	TVector3 m_vec(0,0,0);
	double rotAngle = -(barShowerXCol[0].Bars)[0].module*PI/4.;
	TVector3 Cblock((barShowerXCol[0].Bars)[0].position.x(), (barShowerXCol[0].Bars)[0].position.y(), (barShowerYCol[0].Bars)[0].position.z());

//cout<<"Block center: "<<Cblock.x()<<'\t'<<Cblock.y()<<'\t'<<Cblock.z()<<endl;

	Cblock.RotateZ(rotAngle);

//cout<<"Block center after rotation: "<<Cblock.x()<<'\t'<<Cblock.y()<<'\t'<<Cblock.z()<<endl;

	for(int ix=0;ix<Nshower;ix++){
	for(int iy=0;iy<Nshower;iy++){
		CRDEcalDigiEDM::BarCollection showerX = barShowerXCol[ix];
		CRDEcalDigiEDM::BarCollection showerY = barShowerYCol[iy];

//cout<<"Shower information in ix: "<<ix<<"and iy: "<<iy<<" (x, y, z, E, T1, T2)"<<endl;
//cout<<"Shower X: "<<showerX.getPos().x()<<'\t'<<showerX.getPos().y()<<'\t'<<showerX.getPos().z()<<'\t'<<showerX.getE()<<'\t'<<showerX.getT1()<<'\t'<<showerX.getT2()<<endl;
//cout<<"Shower Y: "<<showerY.getPos().x()<<'\t'<<showerY.getPos().y()<<'\t'<<showerY.getPos().z()<<'\t'<<showerY.getE()<<'\t'<<showerY.getT1()<<'\t'<<showerY.getT2()<<endl;


		double Ex = showerX.getE();
		double Ey = showerY.getE();
		chi2_E[ix][iy] = pow(fabs(Ex-Ey)/sigmaE, 2);

//cout<<"chi2E: "<<chi2_E[ix][iy]<<endl;

		double PosTx = C*(showerY.getT1()-showerY.getT2())/(2*nMat) + showerY.getPos().z();
		chi2_tx[ix][iy] = pow( fabs(PosTx-showerX.getPos().z())/sigmaPos, 2);

		double PosTy = C*(showerX.getT1()-showerX.getT2())/(2*nMat);
		m_vec.SetXYZ(showerY.getPos().x(), showerY.getPos().y(), showerY.getPos().z());
		m_vec.RotateZ(rotAngle);
		chi2_ty[ix][iy] = pow( fabs(PosTy - (m_vec-Cblock).x() )/sigmaPos, 2);

//cout<<"Shower vec after rot: "<<m_vec.x()<<'\t'<<m_vec.y()<<'\t'<<m_vec.z()<<endl;
//cout<<"posTx, posTy and barx: "<<PosTx<<'\t'<<PosTy<<'\t'<<(m_vec-Cblock).x()<<endl;
//cout<<"chi2_tx and chi2_ty: "<<chi2_tx[ix][iy]<<'\t'<<chi2_ty[ix][iy]<<endl;
		chi2[ix][iy] = chi2_E[ix][iy]*wi_E + (chi2_tx[ix][iy]+chi2_ty[ix][iy])*wi_T ;

		m_chi2.push_back(chi2[ix][iy]);
		m_chi2E.push_back(chi2_E[ix][iy]);
		m_chi2Tx.push_back(chi2_tx[ix][iy]);
		m_chi2Ty.push_back(chi2_ty[ix][iy]);
	}}
//cout<<endl;

	int Ncomb=1;
	for(int i=Nshower; i>0; i--) Ncomb = Ncomb*i;

//cout<<"Nshower and Ncomb: "<<Nshower<<'\t'<<Ncomb<<endl;

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
		m_chi2comb.push_back(chi2_tot);

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
