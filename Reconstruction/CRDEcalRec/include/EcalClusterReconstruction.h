#ifndef ECAL_CLUSTER_REC_H
#define ECAL_CLUSTER_REC_H

#include "GaudiAlg/GaudiAlgorithm.h"
#include "PandoraPlusDataCol.h"
#include "Algorithm/EnergySplittingAlg.h"
#include "Algorithm/EnergyTimeMatchingAlg.h"
#include "Algorithm/ConeClusteringAlg.h"

class EcalClusterReconstruction {
public: 

  class Settings{
  public: 
    Settings(){};
  };

  EcalClusterReconstruction() {};
  ~EcalClusterReconstruction();

  void Initialize();
  StatusCode RunAlgorithm( Settings& settings,  PandoraPlusDataCol& dataCol );
  StatusCode ClearAlgorithm();

  Settings   m_settings;


  EnergySplittingAlg     *m_energysplittingAlg;
  EnergyTimeMatchingAlg  *m_etmatchingAlg;
  ConeClusteringAlg      *m_coneclusterAlg;

  EnergySplittingAlg::Settings    *m_ESAlgSettings;
  EnergyTimeMatchingAlg::Settings *m_ETAlgSettings;
  ConeClusteringAlg::Settings     *m_CCAlgSettings;

private: 

};
#endif
