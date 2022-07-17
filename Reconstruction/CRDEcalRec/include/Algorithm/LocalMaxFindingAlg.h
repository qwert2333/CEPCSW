#ifndef _LOCALMAXFINDING_ALG_H
#define _LOCALMAXFINDING_ALG_H

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
  StatusCode Initialize();
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode GetLocalMax( PandoraPlus::Calo2DCluster* m_2dClus);
  StatusCode GetLocalMaxBar( std::vector<const PandoraPlus::CaloUnit*>& barCol, std::vector<const PandoraPlus::CaloUnit*>& localMaxCol );
  std::vector<const PandoraPlus::CaloUnit*>  getNeighbors(const PandoraPlus::CaloUnit* seed, std::vector<const PandoraPlus::CaloUnit*>& barCol);

private: 

  static bool compBar( const PandoraPlus::CaloUnit* bar1, const PandoraPlus::CaloUnit* bar2 )
      { return bar1->getBar() < bar2->getBar(); }

};
#endif
