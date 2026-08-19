#ifndef _TRUTHENERGYSPLITTING_ALG_H
#define _TRUTHENERGYSPLITTING_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;
class TruthEnergySplittingAlg: public PandoraPlus::Algorithm{
public: 

  TruthEnergySplittingAlg(){};
  ~TruthEnergySplittingAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new TruthEnergySplittingAlg(); } 

  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode HalfClusterToTowers( std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& m_halfClusU,
                                  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& m_halfClusV,
                                  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>& m_towers );

private: 
  std::vector<PandoraPlus::CaloHalfCluster*> p_HalfClusterU;
  std::vector<PandoraPlus::CaloHalfCluster*> p_HalfClusterV;

  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_newClusUCol;
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_newClusVCol;
  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> m_towerCol;

  PandoraPlusDataCol m_bkCol;
};
#endif
