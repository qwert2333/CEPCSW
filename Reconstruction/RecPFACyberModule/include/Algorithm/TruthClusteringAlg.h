#ifndef _TRUTHCLUS_ALG_H
#define _TRUTHCLUS_ALG_H

#include "Tools/Algorithm.h"

using namespace Cyber;
class TruthClusteringAlg: public Cyber::Algorithm{
public: 

  TruthClusteringAlg(){};
  ~TruthClusteringAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
  public: 
    Cyber::Algorithm* CreateAlgorithm() const{ return new TruthClusteringAlg(); } 

  };

  StatusCode ReadSettings(Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode HalfClusterToTowers( std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& m_halfClusU,
                                  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& m_halfClusV,
                                  std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_towers );

private: 

  CyberDataCol m_bkCol;

  std::vector<std::shared_ptr<Cyber::Track>> m_TrackCol;
  //For ECAL
  std::vector<std::shared_ptr<Cyber::CaloUnit>> m_bars;
  std::vector<std::shared_ptr<Cyber::Calo1DCluster>> m_1dclusterUCol;
  std::vector<std::shared_ptr<Cyber::Calo1DCluster>> m_1dclusterVCol;
  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>> m_halfclusterU;
  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>> m_halfclusterV;
  std::vector<std::shared_ptr<Cyber::Calo3DCluster>> m_towers;
  //For HCAL
  std::vector<std::shared_ptr<Cyber::CaloHit>> m_hits;
  std::vector<std::shared_ptr<Cyber::Calo3DCluster>> m_clusters;

};
#endif
