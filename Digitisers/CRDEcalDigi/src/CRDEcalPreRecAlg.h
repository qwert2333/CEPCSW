#ifndef _CRDECAL_PRERECALG_
#define _CRDECAL_PRERECALG_
#include "CRDEcalDigiEDM.h"
#include "CRDEcalDigiAlg.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/Cluster.h"

#include "TVector3.h"
#include "TVector.h"
#include "TMatrix.h"
using namespace std;

class CRDEcalPreRecAlg{

public: 
	CRDEcalPreRecAlg();
	~CRDEcalPreRecAlg();

	std::vector<CRDEcalDigiEDM::BarCollection>  Bars2Shower( std::vector<CRDEcalDigiEDM::DigiBar>& barCol );
	std::vector<CRDEcalDigiEDM::BarCluster> Clustering(std::vector<CRDEcalDigiEDM::DigiBar>& barCol); 
	std::vector<CRDEcalDigiEDM::DigiBar>  findSeeds( CRDEcalDigiEDM::BarCluster& m_cluster );
	std::vector<CRDEcalDigiEDM::DigiBar>  getNeighbors(CRDEcalDigiEDM::BarCluster& m_cluster, CRDEcalDigiEDM::DigiBar& seed);
	std::vector<CRDEcalDigiEDM::BarCollection>  ClusterSplitting( CRDEcalDigiEDM::BarCluster& m_cluster);
	void RemoveBars(CRDEcalDigiEDM::BarCluster& m_cluster, CRDEcalDigiEDM::BarCollection& m_barCol);
	void RemoveBars(std::vector<CRDEcalDigiEDM::DigiBar>& m_bars, std::vector<CRDEcalDigiEDM::DigiBar>& m_subbars);
	void CalculateInitialEseed(std::vector<CRDEcalDigiEDM::DigiBar>& Seeds, dd4hep::Position* pos, double* Eseed);
	double GetShowerProfile(dd4hep::Position& p_bar, dd4hep::Position& p_seed );

	double Sth_split;
	double Eth_SeedAbs;
	double Eth_ShowerAbs;
	double Eth_ClusAbs;
	double Eth_SeedWithNeigh;
	double Eth_SeedWithTot;
	double Eth_ShowerWithTot;
	double Eth_ClusterWithTot;
	double Etot;

private: 

};
#endif
