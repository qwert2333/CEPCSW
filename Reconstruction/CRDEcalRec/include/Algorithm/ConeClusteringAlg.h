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

  StatusCode LongiConeLinking( const std::map<int, std::vector<PandoraPlus::CaloHit*> >& orderedShower, std::vector<PandoraPlus::CaloCluster*>& ClusterCol );
  StatusCode MergeGoodClusters( std::vector<PandoraPlus::CaloCluster*>& m_clusCol); 
  StatusCode MergeBadToGoodCluster( std::vector<PandoraPlus::CaloCluster*>& m_goodClusCol, PandoraPlus::CaloCluster* m_badClus );
  PandoraPlus::CaloCluster* GetClosestGoodCluster( std::vector< PandoraPlus::CaloCluster* >& m_goodClusCol, PandoraPlus::CaloCluster* m_badClus );

  static bool compBegin( PandoraPlus::CaloCluster* clus1, PandoraPlus::CaloCluster* clus2 ) { return clus1->getBeginningDlayer() < clus2->getBeginningDlayer(); }

private: 


};
#endif
