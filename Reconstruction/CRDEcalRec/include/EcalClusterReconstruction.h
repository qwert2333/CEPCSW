#ifndef ECAL_CLUSTER_REC_H
#define ECAL_CLUSTER_REC_H

#include "GaudiAlg/GaudiAlgorithm.h"
#include "PandoraPlusDataCol.h"
#include "Algorithm/EnergySplittingAlg.h"
#include "Algorithm/ConeClustering2DAlg.h"
#include "Algorithm/HoughClusteringAlg.h"
#include "Algorithm/ShadowMakingAlg.h"
#include "Algorithm/EnergyTimeMatchingAlg.h"
#include "Algorithm/ConeClusteringAlg.h"
#include "Algorithm/ClusterMergingAlg.h"
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
  ConeClustering2DAlg    *m_coneclus2DAlg;
  HoughClusteringAlg     *m_houghclusAlg; 
  ShadowMakingAlg        *m_shadowmakingAlg;
  EnergyTimeMatchingAlg  *m_etmatchingAlg;
  ConeClusteringAlg      *m_coneclusterAlg;
  ClusterMergingAlg      *m_clustermergingAlg;
  BasicClusterIDAlg      *m_clusteridAlg;


  EnergySplittingAlg::Settings    *m_ESAlgSettings;
  ConeClustering2DAlg::Settings   *m_CC2AlgSettings;
  HoughClusteringAlg::Settings    *m_HCAlgSettings; 
  ShadowMakingAlg::Settings       *m_CMAlgSettings;
  EnergyTimeMatchingAlg::Settings *m_ETAlgSettings;
  ConeClusteringAlg::Settings     *m_CCAlgSettings;
  ClusterMergingAlg::Settings     *m_CLAlgSettings;
  BasicClusterIDAlg::Settings     *m_CIDAlgSettings;


private: 


};
#endif
