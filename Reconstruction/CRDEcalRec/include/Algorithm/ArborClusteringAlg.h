#ifndef _ARBORCLUSTERING_ALG_H
#define _ARBORCLUSTERING_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;
class ArborClusteringAlg: public PandoraPlus::Algorithm{
public: 

  ArborClusteringAlg(){};
  ~ArborClusteringAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new ArborClusteringAlg(); } 

  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize();
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  //Self defined algorithms
  double GetTrkCaloHitDistance(const PandoraPlus::Track* trk, const PandoraPlus::CaloHit* hit); 

  StatusCode InitArborTree( std::map<int, std::vector<ArborNode*> >& m_orderedNodes,
                            std::vector<ArborTree>& m_treeCol,
                            std::vector<ArborNode*>& m_isoNodes );

  StatusCode MergeConnectedTrees( std::vector<ArborTree>& m_inTreeCol, std::vector<ArborTree>& m_outTreeCol );

  StatusCode CleanConnection( ArborTree& m_tree );

  StatusCode DepartArborTree( ArborTree& m_tree, std::vector<ArborTree>& m_departedTrees, std::vector<ArborNode*>& m_isoNodes );

  StatusCode MergeNeighborTree( std::vector<ArborTree>& m_inTreeCol, std::vector<ArborTree>& m_outTreeCol );

private: 

};
#endif
