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

  StatusCode GetLocalMax( PandoraPlus::CaloBlock* m_block);
  StatusCode GetLocalMaxBar( std::vector<const PandoraPlus::CaloBar*>& barCol, std::vector<const PandoraPlus::CaloBar*>& localMaxCol );
  std::vector<const PandoraPlus::CaloBar*>  getNeighbors(const PandoraPlus::CaloBar* seed, std::vector<const PandoraPlus::CaloBar*>& barCol);

private: 

};
#endif
