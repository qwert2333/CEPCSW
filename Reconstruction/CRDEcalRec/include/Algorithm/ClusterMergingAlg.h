#ifndef _CLUSTERMERGING_ALG_H
#define _CLUSTERMERGING_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;
class ClusterMergingAlg: public PandoraPlus::Algorithm{
public: 

  ClusterMergingAlg(){};
  ~ClusterMergingAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new ClusterMergingAlg(); } 

  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  //Self defined algorithms
  StatusCode SelfAlg1(); 

private: 

};
#endif
