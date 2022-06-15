#ifndef _EXAMPLE_ALG_H
#define _EXAMPLE_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;
class ExampleAlg: public PandoraPlus::Algorithm{
public: 

  ExampleAlg(){};
  ~ExampleAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new ExampleAlg(); } 

  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize();
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  //Self defined algorithms
  StatusCode SelfAlg1(); 

private: 

};
#endif
