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
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm(); 


  StatusCode LongitudinalLinking( std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>>& m_showers, 
                                  std::vector<const PandoraPlus::CaloHalfCluster*>& m_oldClusCol, 
                                  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& m_newClusCol );

  StatusCode HalfClusterToTowers( std::vector<PandoraPlus::CaloHalfCluster*>& m_halfClusU, 
                                  std::vector<PandoraPlus::CaloHalfCluster*>& m_halfClusV, 
                                  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>& m_towers,  
                                  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& bk_HFclus,
                                  std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>>& bk_1Dclus,
                                  std::vector<std::shared_ptr<PandoraPlus::Calo2DCluster>>& bk_2Dclus ); 

  StatusCode ClusterSplitting( const PandoraPlus::Calo1DCluster* m_cluster, 
                               std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>>& outshCol, 
                               std::vector<std::shared_ptr<PandoraPlus::CaloUnit>>& bk_bars );

  StatusCode MergeToClosestCluster( PandoraPlus::Calo1DCluster* iclus, std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>>& clusvec );

  StatusCode MergeToClosestCluster( const PandoraPlus::Calo1DCluster* m_shower, std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& m_clusters );

  void CalculateInitialEseed( const std::vector<const PandoraPlus::CaloUnit*>& Seeds, const TVector3* pos, double* Eseed);

  double GetShowerProfile(const TVector3& p_bar, const TVector3& p_seed );


private: 
  std::vector<PandoraPlus::CaloHalfCluster*> p_HalfClusterU;
  std::vector<PandoraPlus::CaloHalfCluster*> p_HalfClusterV;
  std::vector<const PandoraPlus::CaloHalfCluster*> m_axisUCol; 
  std::vector<const PandoraPlus::CaloHalfCluster*> m_axisVCol; 

  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_newClusUCol; 
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_newClusVCol; 
  std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>> m_1dShowerUCol; 
  std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>> m_1dShowerVCol; 
  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> m_towerCol;

  //static bool compBar( const PandoraPlus::CaloUnit* bar1, const PandoraPlus::CaloUnit* bar2 )
  //  { return bar1->getBar() < bar2->getBar(); }
  static bool compLayer( const PandoraPlus::Calo1DCluster* sh1, const PandoraPlus::Calo1DCluster* sh2 )
    { return sh1->getDlayer() < sh2->getDlayer(); }

};

#endif
