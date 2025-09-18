#ifndef _CLUSTERMERGING_ALG_H
#define _CLUSTERMERGING_ALG_H

#include "Tools/Algorithm.h"

using namespace Cyber;
class ClusterMergingAlg: public Cyber::Algorithm{
public: 

  ClusterMergingAlg(){};
  ~ClusterMergingAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
  public: 
    Cyber::Algorithm* CreateAlgorithm() const{ return new ClusterMergingAlg(); } 

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
