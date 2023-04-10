#ifndef _AXISMERGING_ALG_H
#define _AXISMERGING_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;
class AxisMergingAlg: public PandoraPlus::Algorithm{
public: 

  AxisMergingAlg(){};
  ~AxisMergingAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new AxisMergingAlg(); } 

  };

  StatusCode ReadSettings( PandoraPlus::Settings& m_settings );
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  //Self defined algorithms
  StatusCode TrkMatchedMerging( std::vector<PandoraPlus::CaloHalfCluster*>& m_axisCol );
  StatusCode OverlapMerging( std::vector<PandoraPlus::CaloHalfCluster*>& m_axisCol );
  StatusCode ConeMerging( std::vector<PandoraPlus::CaloHalfCluster*>& m_axisCol );
  StatusCode FragmentsMerging( std::vector<PandoraPlus::CaloHalfCluster*>& m_axisCol );
  StatusCode MergeToClosestCluster( PandoraPlus::CaloHalfCluster* m_badaxis, std::vector<PandoraPlus::CaloHalfCluster*>& m_axisCol );

private: 

  std::vector<PandoraPlus::CaloHalfCluster*>* p_HalfClustersU = NULL;
  std::vector<PandoraPlus::CaloHalfCluster*>* p_HalfClustersV = NULL;
  std::vector<const PandoraPlus::CaloHalfCluster*> m_axisUCol;
  std::vector<const PandoraPlus::CaloHalfCluster*> m_axisVCol;
  std::vector<PandoraPlus::CaloHalfCluster*> m_newAxisUCol;
  std::vector<PandoraPlus::CaloHalfCluster*> m_newAxisVCol;

  static bool compLayer( const PandoraPlus::CaloHalfCluster* sh1, const PandoraPlus::CaloHalfCluster* sh2 )
    { return sh1->getBeginningDlayer() < sh2->getBeginningDlayer(); }

};

#endif
