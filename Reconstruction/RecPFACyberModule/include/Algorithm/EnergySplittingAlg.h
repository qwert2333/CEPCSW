#ifndef _ENERGYSPLITTING_ALG_H
#define _ENERGYSPLITTING_ALG_H

#include "CyberDataCol.h"
#include "Tools/Algorithm.h"
#include "TVector3.h"
#include "TVector.h"
#include "TMatrix.h"

using namespace Cyber;
class EnergySplittingAlg: public Cyber::Algorithm{
public: 

  class Factory : public Cyber::AlgorithmFactory
  {
  public:
    Cyber::Algorithm* CreateAlgorithm() const{ return new EnergySplittingAlg(); }

  };

  EnergySplittingAlg(){};
  ~EnergySplittingAlg(){};

  StatusCode ReadSettings( Cyber::Settings& m_settings );
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm(); 


  StatusCode LongitudinalLinking( std::vector<std::shared_ptr<Cyber::Calo1DCluster>>& m_showers, 
                                  std::vector<const Cyber::CaloHalfCluster*>& m_oldClusCol, 
                                  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& m_newClusCol );

  StatusCode HalfClusterToTowers( std::vector<Cyber::CaloHalfCluster*>& m_halfClusU, 
                                  std::vector<Cyber::CaloHalfCluster*>& m_halfClusV, 
                                  std::vector<Cyber::CaloHalfCluster*>& m_emptyClusU,
                                  std::vector<Cyber::CaloHalfCluster*>& m_emptyClusV,
                                  std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_towers );

  StatusCode ClusterSplitting( const Cyber::Calo1DCluster* m_cluster, 
                               std::vector<std::shared_ptr<Cyber::Calo1DCluster>>& outshCol );

  StatusCode SplitOverlapCluster( std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& m_HFClusCol );

  StatusCode MergeToClosestCluster( Cyber::Calo1DCluster* iclus, std::vector<std::shared_ptr<Cyber::Calo1DCluster>>& clusvec );

  StatusCode MergeToClosestCluster( const Cyber::Calo1DCluster* m_shower, std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& m_clusters );

  void CalculateInitialEseed( const std::vector<const Cyber::CaloUnit*>& Seeds, const TVector3* pos, double* Eseed);

  double GetShowerProfile(const TVector3& p_bar, const TVector3& p_seed );


private: 
  std::vector<Cyber::CaloHalfCluster*> p_HalfClusterU;
  std::vector<Cyber::CaloHalfCluster*> p_HalfClusterV;
  std::vector<Cyber::CaloHalfCluster*> p_emptyHalfClusterU;
  std::vector<Cyber::CaloHalfCluster*> p_emptyHalfClusterV;  
  std::vector<const Cyber::CaloHalfCluster*> m_axisUCol; 
  std::vector<const Cyber::CaloHalfCluster*> m_axisVCol; 

  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>> m_newClusUCol; 
  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>> m_newClusVCol; 
  std::vector<std::shared_ptr<Cyber::Calo1DCluster>> m_1dShowerUCol; 
  std::vector<std::shared_ptr<Cyber::Calo1DCluster>> m_1dShowerVCol; 
  std::vector<std::shared_ptr<Cyber::Calo3DCluster>> m_towerCol;

  CyberDataCol m_bkCol;

  //static bool compBar( const Cyber::CaloUnit* bar1, const Cyber::CaloUnit* bar2 )
  //  { return bar1->getBar() < bar2->getBar(); }
  static bool compLayer( const Cyber::Calo1DCluster* sh1, const Cyber::Calo1DCluster* sh2 )
    { return sh1->getDlayer() < sh2->getDlayer(); }

};

#endif
