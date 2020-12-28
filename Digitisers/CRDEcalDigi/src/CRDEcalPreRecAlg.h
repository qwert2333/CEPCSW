#ifndef _CRDECAL_PRERECALG_
#define _CRDECAL_PRERECALG_
#include "CRDEcalDigiEDM.h"
#include "CRDEcalDigiAlg.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/Cluster.h"

#include "TVector3.h"
using namespace std;

class CRDEcalPreRecAlg{

public: 
	CRDEcalPreRecAlg();
	~CRDEcalPreRecAlg();

	std::vector<CRDEcalDigiEDM::BarCollection>  Bars2Shower( std::vector<CRDEcalDigiEDM::DigiBar>& barCol );
	std::vector<CRDEcalDigiEDM::BarCluster> Clustering(std::vector<CRDEcalDigiEDM::DigiBar>& barCol); 
	CRDEcalDigiEDM::BarCollection  findSeed( CRDEcalDigiEDM::BarCluster& m_cluster );
	CRDEcalDigiEDM::BarCollection  getNeighbors(CRDEcalDigiEDM::BarCluster& m_cluster, CRDEcalDigiEDM::BarCollection& seed);
	std::vector<CRDEcalDigiEDM::BarCollection>  ClusterSplitting( CRDEcalDigiEDM::BarCluster& m_cluster);
	vector<int> getShowerEdge(CRDEcalDigiEDM::BarCluster& m_cluster, std::vector<CRDEcalDigiEDM::BarCollection>& m_seeds);
	void RemoveBars(CRDEcalDigiEDM::BarCluster& m_cluster, CRDEcalDigiEDM::BarCollection& m_barCol);
	void RemoveBars(std::vector<CRDEcalDigiEDM::DigiBar>& m_bars, std::vector<CRDEcalDigiEDM::DigiBar>& m_subbars);
	
	double Eth_SeedWithNeigh;
	double Eth_SeedWithTot;
	double Eth_ShowerWithTot;
	double Eth_ClusterWithTot;
	double Etot;

private: 

};
#endif
