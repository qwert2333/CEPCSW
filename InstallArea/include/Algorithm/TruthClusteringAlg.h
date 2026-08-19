#ifndef _TRUTHCLUS_ALG_H
#define _TRUTHCLUS_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;
class TruthClusteringAlg: public PandoraPlus::Algorithm{
public: 

  TruthClusteringAlg(){};
  ~TruthClusteringAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new TruthClusteringAlg(); } 

  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode HalfClusterToTowers( std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& m_halfClusU,
                                  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& m_halfClusV,
                                  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>& m_towers );

private: 

  PandoraPlusDataCol m_bkCol;

  std::vector<std::shared_ptr<PandoraPlus::Track>> m_TrackCol;
  //For ECAL
  std::vector<std::shared_ptr<PandoraPlus::CaloUnit>> m_bars;
  std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>> m_1dclusterUCol;
  std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>> m_1dclusterVCol;
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_halfclusterU;
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_halfclusterV;
  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> m_towers;
  //For HCAL
  std::vector<std::shared_ptr<PandoraPlus::CaloHit>> m_hits;
  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> m_clusters;

};
#endif
