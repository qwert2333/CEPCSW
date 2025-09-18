#ifndef _TRUTHPATTERNREC_ALG_H
#define _TRUTHPATTERNREC_ALG_H

#include "Tools/Algorithm.h"

using namespace Cyber;
class TruthPatternRecAlg: public Cyber::Algorithm{
public: 

  TruthPatternRecAlg(){};
  ~TruthPatternRecAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
  public: 
    Cyber::Algorithm* CreateAlgorithm() const{ return new TruthPatternRecAlg(); } 

  };

  StatusCode ReadSettings(Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode OverlapMerging   ( std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& m_axisCol );

private: 

  std::vector<Cyber::CaloHalfCluster*> p_HalfClusterV;
  std::vector<Cyber::CaloHalfCluster*> p_HalfClusterU;

};
#endif
