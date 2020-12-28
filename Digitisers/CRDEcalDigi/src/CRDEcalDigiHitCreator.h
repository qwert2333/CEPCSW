#ifndef _CRDECAL_DIGIHITCREATOR_
#define _CRDECAL_DIGIHITCREATOR_

#include "CRDEcalDigiAlg.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/Cluster.h"

#include "TVector3.h"
using namespace std;

class CRDEcalDigiHitCreator{

public: 
	CRDEcalDigiHitCreator();
	~CRDEcalDigiHitCreator();
	

	std::vector<edm4hep::Cluster>        Clustering(std::vector<edm4hep::CalorimeterHit>& m_hitcol);
	std::vector<edm4hep::CalorimeterHit> findSeed(std::vector<edm4hep::CalorimeterHit>& m_hitcol);
	std::vector<edm4hep::CalorimeterHit> getNeighbors(std::vector<edm4hep::CalorimeterHit>& m_hitcol, std::vector<edm4hep::CalorimeterHit>& m_seeds);
	std::vector<edm4hep::ConstCalorimeterHit> SplitToHits(std::vector<edm4hep::Cluster>& m_clusters);
	void RemoveHits(std::vector<edm4hep::CalorimeterHit>& m_hitcol, edm4hep::Cluster& clus);

private: 

	mutable float _Etot;
	float _Eth_FormCluster;
	float _Eth_SeedWithTot;
	float _Eth_SeedWithNeigh;
};
#endif
