#ifndef _TRUTHTRKMATCHING_ALG_H
#define _TRUTHTRKMATCHING_ALG_H

#include "Tools/Algorithm.h"

using namespace Cyber;
class TruthTrackMatchingAlg: public Cyber::Algorithm{
public: 

  TruthTrackMatchingAlg(){};
  ~TruthTrackMatchingAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
  public: 
    Cyber::Algorithm* CreateAlgorithm() const{ return new TruthTrackMatchingAlg(); } 

  };

  StatusCode ReadSettings(Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  //Self defined algorithms
  StatusCode SelfAlg1(); 

private: 

  std::vector<Cyber::Track*> m_TrackCol;
  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>* p_HalfClusterV;
  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>* p_HalfClusterU;  



};
#endif
