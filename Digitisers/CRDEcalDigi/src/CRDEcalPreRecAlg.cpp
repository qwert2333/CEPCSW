#include <iostream>
#include <algorithm>
#include "CRDEcalPreRecAlg.h"
#include "CRDEcalDigiEDM.h"
using namespace std;

CRDEcalPreRecAlg::CRDEcalPreRecAlg(){
	Eth_SeedWithNeigh = 0.4;
	Eth_SeedWithTot = 0.15;
	Eth_ShowerWithTot = 0.05;
	Eth_ClusterWithTot = 0.05;
	Etot=0;
}

CRDEcalPreRecAlg::~CRDEcalPreRecAlg(){
}

std::vector<CRDEcalDigiEDM::BarCollection>  CRDEcalPreRecAlg::Bars2Shower( std::vector<CRDEcalDigiEDM::DigiBar>& barCol){

	for(int i=0;i<barCol.size();i++) Etot+=barCol[i].getEnergy();

	std::vector<CRDEcalDigiEDM::BarCollection> m_showers; m_showers.clear();
	std::vector<CRDEcalDigiEDM::BarCluster> m_clusCol; m_clusCol.clear();
	m_clusCol = CRDEcalPreRecAlg::Clustering(barCol); 

	for(int i=0;i<m_clusCol.size();i++){
		if(m_clusCol[i].getE()/Etot < Eth_ClusterWithTot) continue;
		std::vector<CRDEcalDigiEDM::BarCollection> showers = ClusterSplitting(m_clusCol[i]);
		m_showers.insert(m_showers.end(), showers.begin(), showers.end());
	}

	return m_showers;
}

std::vector<CRDEcalDigiEDM::BarCluster> CRDEcalPreRecAlg::Clustering(std::vector<CRDEcalDigiEDM::DigiBar>& barCol){

	for(int i=0;i<barCol.size();i++) Etot+=barCol[i].getEnergy();
	std::sort(barCol.begin(), barCol.end());

	std::vector<CRDEcalDigiEDM::BarCluster> m_clusCol; m_clusCol.clear();
	
	for(int i=0;i<barCol.size();i++){
		CRDEcalDigiEDM::DigiBar iBar = barCol[i];

		if((m_clusCol.size()!=0) && (m_clusCol[m_clusCol.size()-1].isNeighbor(iBar)) ){
			m_clusCol[m_clusCol.size()-1].Bars.push_back(iBar);
			continue;
		}

		CRDEcalDigiEDM::BarCluster clus;
		clus.Bars.push_back(iBar);
		m_clusCol.push_back(clus);

	}

	return m_clusCol;
}


std::vector<CRDEcalDigiEDM::BarCollection> CRDEcalPreRecAlg::ClusterSplitting( CRDEcalDigiEDM::BarCluster& m_cluster){

	std::vector<CRDEcalDigiEDM::BarCollection> m_showers; m_showers.clear();
	m_cluster.CRDEcalDigiEDM::BarCluster::sortByPos();
	//std::sort(m_cluster.Bars.begin(), m_cluster.Bars.end());
	std::vector<CRDEcalDigiEDM::BarCollection> m_seedVec; m_seedVec.clear(); 

	CRDEcalDigiEDM::BarCluster tmp_cl = m_cluster;
	CRDEcalDigiEDM::BarCollection m_seed = findSeed(m_cluster);
	while(m_seed.getE()/m_cluster.getE() > Eth_SeedWithTot ){
		m_seedVec.push_back(m_seed);
		CRDEcalDigiEDM::BarCollection m_neighbor = getNeighbors(tmp_cl, m_seed);
		RemoveBars(tmp_cl, m_seed);
		RemoveBars(tmp_cl, m_neighbor);
		m_seed.Clear();
		m_seed = findSeed(tmp_cl);
	}

	int Nseed = m_seedVec.size();
	m_cluster.Nseeds = Nseed;
	if(Nseed<2){ 
		CRDEcalDigiEDM::BarCollection shower;
		shower.Bars = m_cluster.Bars;
		if( shower.getE()/Etot > Eth_ShowerWithTot) m_showers.push_back(shower);
		return m_showers;
	}

	vector<int> index = getShowerEdge(m_cluster, m_seedVec); //array of shower edge. Size should be Nseed+1. 
	sort(index.begin(), index.end());

	int id0 = m_cluster.Bars[0].bar;
	for(int ish=0;ish<Nseed; ish++){
		CRDEcalDigiEDM::BarCollection shower;
		for(int i=index[ish];i<index[ish+1];i++){
			shower.Bars.push_back(m_cluster.Bars[i-id0]);
		}
		if( shower.getE()/Etot > Eth_ShowerWithTot) m_showers.push_back(shower); 
	}
	return m_showers;

}


CRDEcalDigiEDM::BarCollection CRDEcalPreRecAlg::findSeed(CRDEcalDigiEDM::BarCluster& m_cluster){

	CRDEcalDigiEDM::BarCollection m_seed;
	CRDEcalDigiEDM::DigiBar seedbar;
	double Emax=-1;
	for(int i=0;i<m_cluster.Bars.size();i++){
		if(m_cluster.Bars[i].getEnergy()>Emax){ seedbar = m_cluster.Bars[i]; Emax = m_cluster.Bars[i].getEnergy();}
	}
	m_seed.Bars.push_back(seedbar);

	CRDEcalDigiEDM::BarCollection m_neighbor = getNeighbors(m_cluster, m_seed);
	while(m_seed.getE()!=0 && m_neighbor.getE()!=0 && ( m_seed.getE()/(m_seed.getE()+m_neighbor.getE())<Eth_SeedWithNeigh )){
		CRDEcalDigiEDM::DigiBar subseed;
		double Esub=-1;
		for(int i=0;i<m_neighbor.Bars.size(); i++){
			if(m_neighbor.Bars[i].getEnergy()>Esub){ subseed = m_neighbor.Bars[i]; Esub = m_neighbor.Bars[i].getEnergy();}
		}
		if(Esub>0){ m_seed.Bars.push_back(subseed);}

		m_neighbor.Clear();
		m_neighbor = getNeighbors(m_cluster, m_seed);
	}	

	return m_seed;
}

vector<int> CRDEcalPreRecAlg::getShowerEdge(CRDEcalDigiEDM::BarCluster& m_cluster, std::vector<CRDEcalDigiEDM::BarCollection>& m_seedVec){
	vector<int> index;

	//m_cluster.CRDEcalDigiEDM::BarCluster::sortByPos();
	std::sort(m_cluster.Bars.begin(), m_cluster.Bars.end());
	vector<int> seedID;
	for(int i=0;i<m_seedVec.size();i++){	
		for(int j=0;j<m_seedVec[i].Bars.size();j++) seedID.push_back(m_seedVec[i].Bars[j].bar);
	}

	int id0 = m_cluster.Bars[0].bar;
	index.push_back(id0);
	std::sort(seedID.begin(), seedID.end());
	for(int iseed=0; iseed<seedID.size()-1; iseed++){

		double minE=999;
		double minIndex=-99;
		for(int ibar = seedID[iseed]; ibar<seedID[iseed+1]; ibar++){	
			CRDEcalDigiEDM::DigiBar m_bar = m_cluster.Bars[ibar-id0];
			if(m_bar.getEnergy()<minE){ minE=m_bar.getEnergy(); minIndex = ibar; }
		}
		if(minIndex>0) index.push_back(minIndex);
	}
	

	index.push_back(m_cluster.Bars[m_cluster.Bars.size()-1].bar);
	return index;
}

CRDEcalDigiEDM::BarCollection CRDEcalPreRecAlg::getNeighbors(CRDEcalDigiEDM::BarCluster& m_cluster, CRDEcalDigiEDM::BarCollection& seed){
	CRDEcalDigiEDM::BarCollection m_neighbor;
	for(int i=0;i<m_cluster.Bars.size();i++){
		if( seed.isNeighbor(m_cluster.Bars[i]) ) m_neighbor.Bars.push_back(m_cluster.Bars[i]);
	}
	if(m_neighbor.Bars.size()>2) std::cout<<"WARNING: more than 2 hits in neighborCol!!"<<std::endl;
	return m_neighbor;
}

void CRDEcalPreRecAlg::RemoveBars(CRDEcalDigiEDM::BarCluster& m_cluster, CRDEcalDigiEDM::BarCollection& m_barCol){
	for(int i=0;i<m_barCol.Bars.size();i++){
		auto iter = std::remove(m_cluster.Bars.begin(), m_cluster.Bars.end(), m_barCol.Bars[i]);
		m_cluster.Bars.erase(iter, m_cluster.Bars.end());
	}
}

void CRDEcalPreRecAlg::RemoveBars(std::vector<CRDEcalDigiEDM::DigiBar>& m_bars, std::vector<CRDEcalDigiEDM::DigiBar>& m_subbars){
	for(int i=0;i<m_subbars.size();i++){
		auto iter = std::remove(m_bars.begin(), m_bars.end(), m_subbars[i]);
		m_bars.erase(iter, m_bars.end());
	}
}
