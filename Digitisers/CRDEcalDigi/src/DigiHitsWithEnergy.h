#ifndef _DIGIHITSWITHENERGY_
#define _DIGIHITSWITHENERGY_

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

std::vector<edm4hep::ConstCalorimeterHit> CRDEcalDigiAlg::DigiHitsWithEnergy( std::vector<CRDEcalDigiEDM::DigiBar>& m_block, std::vector<CRDEcalDigiEDM::BarCollection>& barShowerX, std::vector<CRDEcalDigiEDM::BarCollection>& barShowerY){

std::cout<<"DigiHitsWithEnergy: Init #bar in block: "<<m_block.size()<<endl;
std::cout<<"ShowerX Energy: "; for(int i=0;i<barShowerX.size();i++) cout<<barShowerX[i].getE()<<'\t'; cout<<endl;
std::cout<<"ShowerY Energy: "; for(int i=0;i<barShowerY.size();i++) cout<<barShowerY[i].getE()<<'\t'; cout<<endl;
	std::vector<edm4hep::ConstCalorimeterHit> m_digiCol; m_digiCol.clear();

	std::vector<CRDEcalDigiEDM::BarCollection> showerX_signleE; showerX_signleE.clear();
	std::vector<CRDEcalDigiEDM::BarCollection> showerY_signleE; showerY_signleE.clear();
	for(int i=0;i<barShowerX.size();i++){
		bool f_push=false;
		bool EqualE=false;
		for(int j=0;j<barShowerX.size(); j++){
			if(j==i) continue;
			double E1 = barShowerX[i].getE();
			double E2 = barShowerX[j].getE();
			if(fabs(E1-E2)/min(E1, E2) < _Eth_diff  ){ EqualE = true; break;}
		}
		if(!EqualE) f_push=true;
		if(f_push) showerX_signleE.push_back(barShowerX[i]);
	}
	for(int i=0;i<barShowerY.size();i++){
		bool f_push=false;
		bool EqualE=false;
		for(int j=0;j<barShowerY.size(); j++){
			if(j==i) continue;
			double E1 = barShowerY[i].getE();
			double E2 = barShowerY[j].getE();
			if(fabs(E1-E2)/min(E1, E2) < _Eth_diff  ){ EqualE = true; break;}
		}
		if(!EqualE) f_push=true;
		if(f_push) showerY_signleE.push_back(barShowerY[i]);
	}

std::cout<<"DigiHitsWithEnergy: #Single shower X: "<<showerX_signleE.size()<<std::endl;
std::cout<<"DigiHitsWithEnergy: #Single shower Y: "<<showerY_signleE.size()<<std::endl;


	for(int i=0;i<showerX_signleE.size();i++){
		for(int j=0;j<showerY_signleE.size();j++){
			double E1 = barShowerX[i].getE();
			double E2 = barShowerY[j].getE();
			if(fabs(E1-E2)/min(E1, E2) < _Eth_diff ){
std::cout<<"DigiHitsWithEnergy: Case3. X-Y can match! "<<E1<<'\t'<<E2<<std::endl;
				std::vector<edm4hep::ConstCalorimeterHit> m_hitsInShower; m_hitsInShower.clear();
				std::vector<CRDEcalDigiEDM::DigiBar> m_ShowerBlock; m_ShowerBlock.clear();
				m_ShowerBlock.insert(m_ShowerBlock.end(), showerX_signleE[i].Bars.begin(), showerX_signleE[i].Bars.end());
				m_ShowerBlock.insert(m_ShowerBlock.end(), showerY_signleE[j].Bars.begin(), showerY_signleE[j].Bars.end());
				m_hitsInShower = DigiHitsWithPos(m_ShowerBlock);
				m_digiCol.insert(m_digiCol.end(), m_hitsInShower.begin(), m_hitsInShower.end());

				for(int ii=0;ii<m_ShowerBlock.size();ii++){
					auto iter = std::remove(m_block.begin(), m_block.end(), m_ShowerBlock[ii]);
					m_block.erase(iter, m_block.end());
				}

			}
		}
	}
std::cout<<"DigiHitsWithEnergy: #bar in block after Energy matching: "<<m_block.size()<<endl;

	std::vector<edm4hep::ConstCalorimeterHit> m_hitsLeft; m_hitsLeft.clear();
	m_hitsLeft = DigiHitsWithTime(m_block);
	m_digiCol.insert(m_digiCol.end(), m_hitsLeft.begin(), m_hitsLeft.end());

std::cout<<"End DigiHitsWithEnergy!"<<std::endl;
	return m_digiCol;

}

#endif
