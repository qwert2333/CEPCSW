#ifndef _ARBORCLUSTERING_ALG_H
#define _ARBORCLUSTERING_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;
class ArborClusteringAlg: public PandoraPlus::Algorithm{
public: 

  ArborClusteringAlg(){};
  ~ArborClusteringAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new ArborClusteringAlg(); } 

  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize();
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  //Self defined algorithms
  double GetTrkCaloHitDistance(const PandoraPlus::Track* trk, const PandoraPlus::CaloHit* hit); 


private: 

};
#endif
