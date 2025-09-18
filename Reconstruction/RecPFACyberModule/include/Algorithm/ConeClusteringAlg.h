#ifndef _CONECLUSTERING_ALG_H
#define _CONECLUSTERING_ALG_H

#include "CyberDataCol.h"
#include "Tools/Algorithm.h"
#include "TMath.h"
using namespace Cyber;

class ConeClusteringAlg: public Cyber::Algorithm {
public: 

  ConeClusteringAlg(){};
  ~ConeClusteringAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
    Cyber::Algorithm* CreateAlgorithm() const{ return new ConeClusteringAlg(); }
  };

  StatusCode ReadSettings(Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol);
  StatusCode ClearAlgorithm(); 

  StatusCode LongiConeLinking( const std::map<int, std::vector<const Cyber::CaloHit*> >& orderedShower, std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& ClusterCol );
  //StatusCode MergeGoodClusters( std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_clusCol); 
  //StatusCode MergeBadToGoodCluster( std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_goodClusCol, std::shared_ptr<Cyber::Calo3DCluster> m_badClus );
  //Cyber::Calo3DCluster* GetClosestGoodCluster( std::vector< Cyber::Calo3DCluster* >& m_goodClusCol, Cyber::Calo3DCluster* m_badClus );

  //static bool compBegin( Cyber::Calo3DCluster* clus1, Cyber::Calo3DCluster* clus2 ) { return clus1->getBeginningDlayer() < clus2->getBeginningDlayer(); }

private: 


};
#endif
