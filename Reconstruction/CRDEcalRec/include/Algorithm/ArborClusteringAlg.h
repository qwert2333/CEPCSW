#ifndef _ARBORCLUSTERING_ALG_H
#define _ARBORCLUSTERING_ALG_H

#include "PandoraPlusDataCol.h"

using namespace CRDEcalEDM;

class ArborClusteringAlg{

public: 
  class Settings{
  public: 
    Settings() {};
    void SetInitialValue(); 

    double Rth_link;  // Node distance threshold when linking nodes. 

    double wiB;
    double wiF;
    double pTheta;
    double pR;
    double pE; 

    double Rth_nbrRoot; // Root node distance threshold when merging neighbor trees. 

  };

  ArborClusteringAlg() {};
  ~ArborClusteringAlg() {};
  StatusCode Initialize();
  StatusCode RunAlgorithm( ArborClusteringAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);
  StatusCode ClearAlgorithm(); 

  StatusCode InitArborTree( std::map<int, std::vector<CRDEcalEDM::CRDArborNode> >& m_orderedNodes,
                            std::vector<CRDEcalEDM::CRDArborTree>& m_treeCol,
                            std::vector<CRDEcalEDM::CRDArborNode*>& m_isoNodes );
  StatusCode MergeConnectedTrees( std::vector<CRDEcalEDM::CRDArborTree>& m_inTreeCol,
                                  std::vector<CRDEcalEDM::CRDArborTree>& m_outTreeCol );

  StatusCode CleanConnection( CRDEcalEDM::CRDArborTree& m_tree );

//  StatusCode DepartArborTree( CRDEcalEDM::CRDArborTree& m_tree, std::vector<CRDEcalEDM::CRDArborTree>& m_departedTrees, std::vector<CRDEcalEDM::CRDArborNode*>& m_isoNodes);

//  StatusCode MergeNeighborTree( std::vector<CRDEcalEDM::CRDArborTree>& m_inTreeCol, std::vector<CRDEcalEDM::CRDArborTree>& m_outTreeCol );

//  StatusCode MergeBranches( std::vector<CRDEcalEDM::CRDArborTree>& m_inTreeCol, std::vector<CRDEcalEDM::CRDArborTree>& m_outTreeCol );

//  CRDEcalEDM::CRDArborTree GetClosestTree( CRDEcalEDM::CRDArborTree m_badTree, std::vector<CRDEcalEDM::CRDArborTree> m_goodTreeCol );

  Settings settings;

private: 

};

#endif
