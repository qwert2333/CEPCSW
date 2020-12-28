#include <iostream>
#include "CRDEcalDigiHitCreator.h"

using namespace std;

CRDEcalDigiHitCreator::CRDEcalDigiHitCreator(){
	_Eth_FormCluster = 0.005;
	_Eth_SeedWithTot = 0.1;
	_Eth_SeedWithNeigh = 0.7;
	_Etot = 0;
}

CRDEcalDigiHitCreator::~CRDEcalDigiHitCreator(){
}


std::vector<edm4hep::Cluster> CRDEcalDigiHitCreator::Clustering(std::vector<edm4hep::CalorimeterHit>& m_hitcol){

	std::vector<edm4hep::Cluster> m_clusters; m_clusters.clear();
	std::vector<edm4hep::CalorimeterHit> m_seeds; m_seeds.clear();
	std::vector<edm4hep::CalorimeterHit> m_neighbor; m_neighbor.clear();

	_Etot=0;
	for(int ihit=0;ihit<m_hitcol.size();ihit++) _Etot+=m_hitcol[ihit].getEnergy();
	if(_Etot==0){ std::cout<<"ERROR: 0 energy for DigiHit collection!!"<<std::endl; return m_clusters;}

	float Eseed = 0;
	m_seeds = findSeed(m_hitcol);
	for(int i=0;i<m_seeds.size();i++) Eseed += m_seeds.at(i).getEnergy();


	while(Eseed/_Etot > _Eth_SeedWithTot){
		m_neighbor.clear();
		//std::vector<edm4hep::CalorimeterHit>(m_neighbor).swap(m_neighbor);
		m_neighbor = getNeighbors(m_hitcol, m_seeds); 

		edm4hep::Cluster clus;
		for(int ihit=0;ihit<m_seeds.size();ihit++) clus.addToHits(m_seeds.at(ihit));
		for(int ihit=0;ihit<m_neighbor.size();ihit++) clus.addToHits(m_neighbor[ihit]);
		m_clusters.push_back(clus);

		RemoveHits(m_hitcol, clus);
		//_Etot=0;
		//for(int ihit=0;ihit<m_hitcol.size();ihit++) _Etot+=m_hitcol[ihit].getEnergy();

		Eseed=0; m_seeds.clear();
		//std::vector<edm4hep::CalorimeterHit>(m_seeds).swap(m_seeds);
		m_seeds = findSeed(m_hitcol);

		for(int i=0;i<m_seeds.size();i++) Eseed += m_seeds.at(i).getEnergy(); 
	}

	return m_clusters;
}


std::vector<edm4hep::CalorimeterHit> CRDEcalDigiHitCreator::findSeed(std::vector<edm4hep::CalorimeterHit>& m_hitcol){

	std::vector<edm4hep::CalorimeterHit> m_seeds;
	edm4hep::CalorimeterHit seed;
	float Emax=-1;
	for(int i=0;i<m_hitcol.size();i++){
		if(m_hitcol[i].getEnergy()>Emax) {seed = m_hitcol[i]; Emax=seed.getEnergy();}
	}
	m_seeds.push_back(seed);

/*	std::vector<edm4hep::CalorimeterHit> neighbors;
	neighbors = getNeighbors(m_hitcol, m_seeds);
	float Eseed=0;
	float Eneigh=0;
	for(int i=0;i<m_seeds.size();i++) Eseed += m_seeds[i].getEnergy();
	for(int i=0;i<neighbors.size();i++) Eneigh += neighbors[i].getEnergy();


	while(Eseed!=0 && Eneigh!=0 && (Eseed/(Eseed+Eneigh)<_Eth_SeedWithNeigh)){ //FIXME: should find in neighbor
		edm4hep::CalorimeterHit subseed;
		float Esub=-1;
		for(int i=0;i<m_hitcol.size();i++){
			if(neighbors[i].getEnergy()>Emax) subseed = m_hitcol[i]; //FIXME: debug: forgot to refresh Emax.
		}
		m_seeds.push_back(subseed);

		neighbors.clear();
		std::vector<edm4hep::CalorimeterHit>(neighbors).swap(neighbors);
		neighbors = getNeighbors(m_hitcol, m_seeds);
		Eseed=0; Eneigh=0;
		for(int i=0;i<m_seeds.size();i++) Eseed += m_seeds[i].getEnergy();
		for(int i=0;i<neighbors.size();i++) Eneigh += neighbors[i].getEnergy();
	}
*/
	return m_seeds;

}

std::vector<edm4hep::CalorimeterHit> CRDEcalDigiHitCreator::getNeighbors(	std::vector<edm4hep::CalorimeterHit>& m_hitcol, 
																																					std::vector<edm4hep::CalorimeterHit>& m_seeds  ){

	std::vector<edm4hep::CalorimeterHit> m_neighbors; m_neighbors.clear();

	for(int ihit=0;ihit<m_hitcol.size();ihit++){
		TVector3 poshit; poshit.SetXYZ(m_hitcol[ihit].getPosition()[0], m_hitcol[ihit].getPosition()[1], m_hitcol[ihit].getPosition()[2]);

		for(int iseed=0;iseed<m_seeds.size();iseed++){
			TVector3 posseed; posseed.SetXYZ(m_seeds[iseed].getPosition()[0], m_seeds[iseed].getPosition()[1], m_seeds[iseed].getPosition()[2]);

			if((poshit-posseed).Mag()>10*sqrt(2.01)) continue;

			for(int jseed=0;jseed<m_seeds.size();jseed++){
				if( m_hitcol[ihit].getEnergy() != m_seeds[jseed].getEnergy() ) m_neighbors.push_back(m_hitcol[ihit]);
			}
		}
	}
	return m_neighbors;
}

std::vector<edm4hep::ConstCalorimeterHit> CRDEcalDigiHitCreator::SplitToHits(std::vector<edm4hep::Cluster>& m_clusters){
	std::vector<edm4hep::ConstCalorimeterHit> m_digihits; m_digihits.clear();
	for(int icl=0;icl<m_clusters.size();icl++){
		edm4hep::Cluster clus = m_clusters[icl];
		for(int ihit=0;ihit<clus.hits_size();ihit++) m_digihits.push_back(clus.getHits(ihit));
	}
	return m_digihits;
}


void CRDEcalDigiHitCreator::RemoveHits(std::vector<edm4hep::CalorimeterHit>& m_hitcol, edm4hep::Cluster& clus){
	std::vector<edm4hep::Cluster> clusvec; clusvec.push_back(clus);
	std::vector<edm4hep::ConstCalorimeterHit> hitincl = SplitToHits(clusvec);
	for(int i=0;i<hitincl.size();i++){
		auto iter = std::remove(m_hitcol.begin(), m_hitcol.end(), hitincl[i]);
		m_hitcol.erase(iter, m_hitcol.end());
	}

}
