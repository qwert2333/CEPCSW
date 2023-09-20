#ifndef _TRUTHMATCHING_ALG_H
#define _TRUTHMATCHING_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;
class TruthMatchingAlg: public PandoraPlus::Algorithm{
public: 

  TruthMatchingAlg(){};
  ~TruthMatchingAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new TruthMatchingAlg(); } 

  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

private: 

};
#endif
