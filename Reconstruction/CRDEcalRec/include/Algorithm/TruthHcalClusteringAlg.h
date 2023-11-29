#ifndef _TRUTHHCALCLUS_ALG_H
#define _TRUTHHCALCLUS_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;
class TruthHcalClusteringAlg: public PandoraPlus::Algorithm{
public: 

  TruthHcalClusteringAlg(){};
  ~TruthHcalClusteringAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new TruthHcalClusteringAlg(); } 

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
