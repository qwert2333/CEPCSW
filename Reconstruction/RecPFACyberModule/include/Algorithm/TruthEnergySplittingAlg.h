#ifndef _TRUTHENERGYSPLITTING_ALG_H
#define _TRUTHENERGYSPLITTING_ALG_H

#include "Tools/Algorithm.h"

using namespace Cyber;
class TruthEnergySplittingAlg: public Cyber::Algorithm{
public: 

  TruthEnergySplittingAlg(){};
  ~TruthEnergySplittingAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
  public: 
    Cyber::Algorithm* CreateAlgorithm() const{ return new TruthEnergySplittingAlg(); } 

  };

  StatusCode ReadSettings(Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode HalfClusterToTowers( std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& m_halfClusU,
                                  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& m_halfClusV,
                                  std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_towers );

private: 
  std::vector<Cyber::CaloHalfCluster*> p_HalfClusterU;
  std::vector<Cyber::CaloHalfCluster*> p_HalfClusterV;

  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>> m_newClusUCol;
  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>> m_newClusVCol;
  std::vector<std::shared_ptr<Cyber::Calo3DCluster>> m_towerCol;

  CyberDataCol m_bkCol;
};
#endif
