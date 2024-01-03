#ifndef _TRUTHPATTERNREC_ALG_H
#define _TRUTHPATTERNREC_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;
class TruthPatternRecAlg: public PandoraPlus::Algorithm{
public: 

  TruthPatternRecAlg(){};
  ~TruthPatternRecAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new TruthPatternRecAlg(); } 

  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode OverlapMerging   ( std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& m_axisCol );

private: 

  std::vector<PandoraPlus::CaloHalfCluster*> p_HalfClusterV;
  std::vector<PandoraPlus::CaloHalfCluster*> p_HalfClusterU;

};
#endif
