#ifndef _EXAMPLE_ALG_H
#define _EXAMPLE_ALG_H

#include "Tools/Algorithm.h"

using namespace Cyber;
class ExampleAlg: public Cyber::Algorithm{
public: 

  ExampleAlg(){};
  ~ExampleAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
  public: 
    Cyber::Algorithm* CreateAlgorithm() const{ return new ExampleAlg(); } 

  };

  StatusCode ReadSettings(Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  //Self defined algorithms
  StatusCode SelfAlg1(); 

private: 

};
#endif
