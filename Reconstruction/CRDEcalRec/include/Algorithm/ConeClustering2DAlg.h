#ifndef _CONECLUSTERING2D_ALG_H
#define _CONECLUSTERING2D_ALG_H

#include "PandoraPlusDataCol.h"
#include "Tools/Algorithm.h"
#include "TMath.h"

using namespace PandoraPlus;
class ConeClustering2DAlg: public PandoraPlus::Algorithm{
public: 

  ConeClustering2DAlg(){};
  ~ConeClustering2DAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new ConeClustering2DAlg(); } 

  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize();
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  //Self defined algorithms
  StatusCode LongiConeLinking( std::map<int, std::vector<const PandoraPlus::CaloBarShower*> >& orderedShower, std::vector<PandoraPlus::LongiCluster*>& ClusterCol);

private: 

};
#endif
