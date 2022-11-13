#ifndef _CONECLUSTERING_ALG_H
#define _CONECLUSTERING_ALG_H

#include "PandoraPlusDataCol.h"
#include "Tools/Algorithm.h"
#include "TMath.h"
using namespace PandoraPlus;

class ConeClusteringAlg: public PandoraPlus::Algorithm {
public: 

  ConeClusteringAlg(){};
  ~ConeClusteringAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new ConeClusteringAlg(); }
  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize();
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol);
  StatusCode ClearAlgorithm(); 

  StatusCode LongiConeLinking( const std::map<int, std::vector<PandoraPlus::CaloHit*> >& orderedShower, std::vector<PandoraPlus::Calo3DCluster*>& ClusterCol );
  StatusCode MergeGoodClusters( std::vector<PandoraPlus::Calo3DCluster*>& m_clusCol); 
  StatusCode MergeBadToGoodCluster( std::vector<PandoraPlus::Calo3DCluster*>& m_goodClusCol, PandoraPlus::Calo3DCluster* m_badClus );
  PandoraPlus::Calo3DCluster* GetClosestGoodCluster( std::vector< PandoraPlus::Calo3DCluster* >& m_goodClusCol, PandoraPlus::Calo3DCluster* m_badClus );


private: 


};
#endif
