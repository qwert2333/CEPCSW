#ifndef _AXISMERGING_ALG_H
#define _AXISMERGING_ALG_H

#include "Tools/Algorithm.h"

using namespace Cyber;
class AxisMergingAlg: public Cyber::Algorithm{
public: 

  AxisMergingAlg(){};
  ~AxisMergingAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
  public: 
    Cyber::Algorithm* CreateAlgorithm() const{ return new AxisMergingAlg(); } 

  };

  StatusCode ReadSettings( Cyber::Settings& m_settings );
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  //Self defined algorithms
  StatusCode TrkMatchedMerging( std::vector<Cyber::CaloHalfCluster*>& m_axisCol );
  StatusCode OverlapMerging   ( std::vector<Cyber::CaloHalfCluster*>& m_axisCol );
  StatusCode BranchMerging   ( std::vector<Cyber::CaloHalfCluster*>& m_axisCol );     // yyy: trying to merge fake photon to track axis
  StatusCode ConeMerging      ( std::vector<Cyber::CaloHalfCluster*>& m_axisCol );
  StatusCode FragmentsMerging ( std::vector<Cyber::CaloHalfCluster*>& m_axisCol );
  bool MergeToClosestCluster( Cyber::CaloHalfCluster* m_badaxis, std::vector<Cyber::CaloHalfCluster*>& m_axisCol );

private: 

  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>* p_HalfClusterU = nullptr;
  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>* p_HalfClusterV = nullptr;

  std::vector<const Cyber::CaloHalfCluster*> m_axisUCol;
  std::vector<const Cyber::CaloHalfCluster*> m_axisVCol;
  std::vector<Cyber::CaloHalfCluster*> m_newAxisUCol;
  std::vector<Cyber::CaloHalfCluster*> m_newAxisVCol;

  static bool compLayer( const Cyber::CaloHalfCluster* sh1, const Cyber::CaloHalfCluster* sh2 )
    { if( sh1->getBeginningDlayer() != sh2->getBeginningDlayer() )
        return sh1->getBeginningDlayer() < sh2->getBeginningDlayer();
      else
        return sh1->getEnergy() > sh2->getEnergy();
    }
  static bool compLayerPtr( const std::shared_ptr<Cyber::CaloHalfCluster> sh1, const std::shared_ptr<Cyber::CaloHalfCluster> sh2 )
    { if( sh1->getBeginningDlayer() != sh2->getBeginningDlayer() )
        return sh1->getBeginningDlayer() < sh2->getBeginningDlayer();
      else
        return sh1->getEnergy() > sh2->getEnergy();
    }


};

#endif
