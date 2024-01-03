#ifndef _LOCALMAXFINDING_ALG_H
#define _LOCALMAXFINDING_ALG_H

#include <set>
#include "Tools/Algorithm.h"
using namespace PandoraPlus;

class LocalMaxFindingAlg: public PandoraPlus::Algorithm{
public: 

  LocalMaxFindingAlg(){};
  ~LocalMaxFindingAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const { return new LocalMaxFindingAlg(); }

  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode GetLocalMax( const PandoraPlus::Calo1DCluster* m_1dClus, std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>>& m_output);
  StatusCode GetLocalMaxBar( std::vector<const PandoraPlus::CaloUnit*>& barCol, std::vector<const PandoraPlus::CaloUnit*>& localMaxCol );
  std::vector<const PandoraPlus::CaloUnit*>  getNeighbors(const PandoraPlus::CaloUnit* seed, std::vector<const PandoraPlus::CaloUnit*>& barCol);

private: 

  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>* p_HalfClusU = nullptr;
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>* p_HalfClusV = nullptr;

  static bool compBar( const PandoraPlus::CaloUnit* bar1, const PandoraPlus::CaloUnit* bar2 )
      { return bar1->getBar() < bar2->getBar(); }

};
#endif
