#ifndef _ENERGYSPLITTING_ALG_H
#define _ENERGYSPLITTING_ALG_H

#include "PandoraPlusDataCol.h"
#include "Tools/Algorithm.h"
#include "TVector3.h"
#include "TVector.h"
#include "TMatrix.h"

using namespace PandoraPlus;
class EnergySplittingAlg: public PandoraPlus::Algorithm{
public: 

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public:
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new EnergySplittingAlg(); }

  };

  EnergySplittingAlg(){};
  ~EnergySplittingAlg(){};

  StatusCode ReadSettings( PandoraPlus::Settings& m_settings );
  StatusCode Initialize();
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm(); 


  StatusCode LongiClusterLinking( std::vector<const PandoraPlus::CaloBlock*>& m_blocks, std::vector<const PandoraPlus::LongiCluster*>& m_oldClusCol, std::vector<const PandoraPlus::LongiCluster*>& m_outClusCol );

  StatusCode Clustering( std::vector<const PandoraPlus::CaloUnit*>& barCol, std::vector<PandoraPlus::CaloBarCluster*>& outClus, std::vector<const PandoraPlus::LongiCluster*>& m_longiClusCol );

  StatusCode Clustering( std::vector<const PandoraPlus::CaloUnit*>& barCol, std::vector<PandoraPlus::CaloBarCluster*>& outClus );

  StatusCode ClusterSplitting( PandoraPlus::CaloBarCluster* m_cluster, std::vector<const PandoraPlus::CaloBarShower*>& outshCol );

  StatusCode MergeToClosestCluster( PandoraPlus::CaloBarCluster* iclus, std::vector<PandoraPlus::CaloBarCluster*>& clusvec );

  StatusCode MergeToClosestCluster( const PandoraPlus::CaloBarShower* m_shower, std::vector<PandoraPlus::LongiCluster*>& m_clusters );

  StatusCode findSeeds( PandoraPlus::CaloBarCluster* m_cluster, std::vector<const PandoraPlus::CaloUnit*>& seedCol );

  std::vector<const PandoraPlus::CaloUnit*>  getNeighbors(PandoraPlus::CaloBarCluster* m_cluster, const PandoraPlus::CaloUnit* seed); 

  void CalculateInitialEseed( const std::vector<const PandoraPlus::CaloUnit*>& Seeds, const TVector3* pos, double* Eseed);

  double GetShowerProfile(const TVector3& p_bar, const TVector3& p_seed );


private: 

  static bool compBar( const PandoraPlus::CaloUnit* bar1, const PandoraPlus::CaloUnit* bar2 )
    { return bar1->getBar() < bar2->getBar(); }
  static bool compLayer( const PandoraPlus::CaloBarShower* sh1, const PandoraPlus::CaloBarShower* sh2 )
    { return sh1->getDlayer() < sh2->getDlayer(); }

};

#endif
