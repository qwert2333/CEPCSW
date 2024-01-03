#ifndef _TRUTHTRKMATCHING_ALG_H
#define _TRUTHTRKMATCHING_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;
class TruthTrackMatchingAlg: public PandoraPlus::Algorithm{
public: 

  TruthTrackMatchingAlg(){};
  ~TruthTrackMatchingAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new TruthTrackMatchingAlg(); } 

  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  //Self defined algorithms
  StatusCode SelfAlg1(); 

private: 

  std::vector<PandoraPlus::Track*> m_TrackCol;
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>* p_HalfClusterV;
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>* p_HalfClusterU;  



};
#endif
