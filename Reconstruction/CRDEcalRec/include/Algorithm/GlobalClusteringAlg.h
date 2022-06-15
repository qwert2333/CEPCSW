#ifndef GLOBALCLUSTERING_ALG_H
#define GLOBALCLUSTERING_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;

class GlobalClusteringAlg : public PandoraPlus::Algorithm{
public: 

  GlobalClusteringAlg(){};
  ~GlobalClusteringAlg(){};
 
  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public:
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new GlobalClusteringAlg(); }
 
  };
 
  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize();
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

private:

};
#endif
