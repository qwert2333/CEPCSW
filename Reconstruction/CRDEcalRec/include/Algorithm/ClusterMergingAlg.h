#ifndef _CLUSTERMERGING_ALG_H
#define _CLUSTERMERGING_ALG_H

#include "PandoraPlusDataCol.h"
#include "Algorithm/EnergySplittingAlg.h"
#include "Algorithm/EnergyTimeMatchingAlg.h"
#include "Algorithm/ConeClusteringAlg.h"

using namespace CRDEcalEDM;

class ClusterMergingAlg{

public: 
  class Settings{
  public: 
    Settings() {};
    void SetInitialValue();
    void SetDebug(int _level) { Debug  =_level; }
    int Debug;
  };

  ClusterMergingAlg() {};
  StatusCode Initialize();
  StatusCode RunAlgorithm( ClusterMergingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);
  StatusCode ClearAlgorithm(); 

  std::vector<int> GetGhostHitsLayer( std::vector<CRDEcalEDM::CRDCaloHit3DShower>& m_goodClus, std::vector<CRDEcalEDM::CRDCaloHit3DShower>& m_badClus );
  std::vector<int> GetLayerNeedModification( std::vector<CRDEcalEDM::CRDCaloHit3DShower>& m_goodClus, std::vector<CRDEcalEDM::CRDCaloHit3DShower>& m_badClus );

  Settings settings;


};

#endif
