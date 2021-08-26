#ifndef ECAL_CLUSTER_REC_H
#define ECAL_CLUSTER_REC_H

#include "GaudiAlg/GaudiAlgorithm.h"
#include "PandoraPlusDataCol.h"
#include "Algorithm/EnergySplittingAlg.h"
#include "Algorithm/EnergyTimeMatchingAlg.h"
#include "Algorithm/ConeClusteringAlg.h"
#include "Algorithm/ShadowMakingAlg.h"
#include "Algorithm/ClusterMergingAlg.h"
#include "Algorithm/ArborClusteringAlg.h"
#include "Algorithm/ArborTreeMergingAlg.h"
#include "Algorithm/BasicClusterIDAlg.h"
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
  ShadowMakingAlg        *m_shadowmakingAlg;
  ClusterMergingAlg      *m_clustermergingAlg; 
  ArborClusteringAlg     *m_arborclusteringAlg; 
  ArborTreeMergingAlg    *m_arbortreemergingAlg; 
  BasicClusterIDAlg      *m_clusteridAlg; 

  EnergySplittingAlg::Settings    *m_ESAlgSettings;
  EnergyTimeMatchingAlg::Settings *m_ETAlgSettings;
  ConeClusteringAlg::Settings     *m_CCAlgSettings;
  ShadowMakingAlg::Settings       *m_CMAlgSettings; 
  ClusterMergingAlg::Settings     *m_CLAlgSettings;  
  ArborClusteringAlg::Settings    *m_ACAlgSettings;
  ArborTreeMergingAlg::Settings   *m_AMAlgSettings; 
  BasicClusterIDAlg::Settings     *m_CIDAlgSettings; 

private: 


};
#endif
