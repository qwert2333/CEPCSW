#ifndef ECAL_CLUSTER_REC_H
#define ECAL_CLUSTER_REC_H

#include "GaudiAlg/GaudiAlgorithm.h"
#include "PandoraPlusDataCol.h"
#include "Algorithm/EnergySplittingAlg.h"
#include "Algorithm/EnergyTimeMatchingAlg.h"
#include "Algorithm/ConeClusteringAlg.h"
#include "Algorithm/CandidateMakingAlg.h"
#include "Algorithm/ClusterMergingAlg.h"
#include "Algorithm/ArborClusteringAlg.h"

class EcalClusterReconstruction {
public: 

  class Settings{
  public: 
    Settings(){};
  };

  EcalClusterReconstruction() {};
  ~EcalClusterReconstruction() {};

  void Initialize();
  StatusCode RunAlgorithm( Settings& settings,  PandoraPlusDataCol& dataCol );
  StatusCode ClearAlgorithm();

  Settings   m_settings;


  EnergySplittingAlg     *m_energysplittingAlg;
  EnergyTimeMatchingAlg  *m_etmatchingAlg;
  ConeClusteringAlg      *m_coneclusterAlg;
  CandidateMakingAlg     *m_candidatemakingAlg;
  ClusterMergingAlg      *m_clustermergingAlg; 
  ArborClusteringAlg     *m_arborclusteringAlg; 


  EnergySplittingAlg::Settings    *m_ESAlgSettings;
  EnergyTimeMatchingAlg::Settings *m_ETAlgSettings;
  ConeClusteringAlg::Settings     *m_CCAlgSettings;
  CandidateMakingAlg::Settings    *m_CMAlgSettings; 
  ClusterMergingAlg::Settings     *m_CLAlgSettings;  
  ArborClusteringAlg::Settings    *m_ACAlgSettings;

private: 


};
#endif
