#ifndef _ENERGYSPLITTING_ALG_H
#define _ENERGYSPLITTING_ALG_H

#include "PandoraPlusDataCol.h"

#include "TVector3.h"
#include "TVector.h"
#include "TMatrix.h"

using namespace CRDEcalEDM;

class EnergySplittingAlg {
public: 

  class Settings {
  public: 
    Settings() {};
    void SetInitialValue();
    void SetOnlyFindMax(bool _fl) { fl_findMaxOnly=_fl; }
    void SetDebug(int _debug) { Debug=_debug; }

    bool fl_findMaxOnly; 
    double Sth_split;
    double Eth_localMax;
    double Eth_SeedAbs;
    double Eth_ShowerAbs;
    double Eth_ClusAbs;
    double Eth_SeedWithNeigh;
    double Eth_MaxWithNeigh;
    double Eth_SeedWithTot;
    double Eth_ShowerWithTot;
    double Eth_ClusterWithTot;
    double Etot;
    int    Debug;

  };

  EnergySplittingAlg(){};
  ~EnergySplittingAlg(){};

  StatusCode Initialize();
  StatusCode RunAlgorithm( EnergySplittingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);
  StatusCode ClearAlgorithm(); 

  StatusCode GetLocalMaxinLayer( CRDEcalEDM::CRDCaloBlock& m_block,CRDEcalEDM::CRDCaloLayer& m_layer );
  StatusCode GetLocalMaxBar( std::vector<CRDEcalEDM::CRDCaloBar>& barCol, std::vector<CRDEcalEDM::CRDCaloBar>& localMaxCol );

  StatusCode ClusteringinLayer( CRDEcalEDM::CRDCaloBlock& m_blocks,CRDEcalEDM::CRDCaloLayer& m_layer );
  bool Clustering( std::vector<CRDEcalEDM::CRDCaloBar>& barCol,  std::vector<CRDEcalEDM::CRDCaloBarCluster>& outClus);
  bool Clustering( std::vector<CRDEcalEDM::CRDCaloBar>& barCol,  std::vector<CRDEcalEDM::CRDCaloBarCluster>& outClus, const std::vector<CRDEcalEDM::CRDShadowCluster>& m_CandiCol );
  bool ClusterSplitting( CRDEcalEDM::CRDCaloBarCluster& m_cluster, std::vector<CRDEcalEDM::CRDCaloBarShower>& outshCol ); 
  bool ClusterSplitting( CRDEcalEDM::CRDCaloBarCluster& m_cluster, std::vector<CRDEcalEDM::CRDCaloBarShower>& outshCol, const std::vector<CRDEcalEDM::CRDShadowCluster>& m_EMcandidates, const std::vector<CRDEcalEDM::CRDShadowCluster>& m_Trkcandidates );
  std::vector<CRDEcalEDM::CRDCaloBar>  getNeighbors( CRDEcalEDM::CRDCaloBar& seed, std::vector<CRDEcalEDM::CRDCaloBar>& barCol);
  std::vector<CRDEcalEDM::CRDCaloBar>  getNeighbors( CRDEcalEDM::CRDCaloBarCluster& m_cluster, CRDEcalEDM::CRDCaloBar& seed);
  std::vector<CRDEcalEDM::CRDCaloBar>  findSeeds( CRDEcalEDM::CRDCaloBarCluster& m_cluster);
  std::vector<CRDEcalEDM::CRDCaloBar>  findSeeds( CRDEcalEDM::CRDCaloBarCluster& m_cluster, const std::vector<CRDEcalEDM::CRDShadowCluster>& m_CandiCol);
  bool GetMatchedSeeds( std::vector<CRDEcalEDM::CRDCaloBar>& seeds, const std::vector<CRDEcalEDM::CRDShadowCluster>& m_CandiCol, std::vector<CRDEcalEDM::CRDCaloBar>& m_matchedCol, std::vector<CRDEcalEDM::CRDCaloBar>& m_unmatchedCol );

  bool MergeToClosestCluster( CRDEcalEDM::CRDCaloBarCluster& iclus, std::vector<CRDEcalEDM::CRDCaloBarCluster>& clusvec );
  void CalculateInitialEseed( const std::vector<CRDEcalEDM::CRDCaloBar>& Seeds, const dd4hep::Position* pos, double* Eseed);
  double GetShowerProfile(const dd4hep::Position& p_bar, const dd4hep::Position& p_seed );
  static bool compE( CRDEcalEDM::CRDCaloBar& bar1, CRDEcalEDM::CRDCaloBar& bar2 ) { return bar1.getEnergy()>bar2.getEnergy(); }



  Settings settings;

private: 

};


#endif
