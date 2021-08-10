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
    void SetRnodeThres( double _value, double _slope ) { Rth_value=_value; Rth_slope=_slope; }
    void SetRefWeight( double _wiF, double _wiB ) { wiF=_wiF; wiB=_wiB; }
    void SetKappaOrderWeight( double _w1, double _w2, double _w3 ) { pTheta=_w1; pR=_w2; pE=_w3; }
    void SetRseedThres( double _th) { Rth_nbrRoot=_th; }
    void SetDebug( int _debug) { Debug=_debug; }
    void PrintSettings() const; 

    double Rth_value;  
    double Rth_slope;  // Node distance pars when linking nodes. R = Rth_value + Rth_slope * Dlayer

    double wiB;
    double wiF;
    double pTheta;
    double pR;
    double pE; 

    double Rth_nbrRoot; // Root node distance threshold when merging neighbor trees. 
    int Debug;

  };

  ArborClusteringAlg() {};
  ~ArborClusteringAlg() {};
  StatusCode Initialize();
  StatusCode RunAlgorithm( ArborClusteringAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);
  StatusCode ClearAlgorithm(); 

  StatusCode InitArborTree( std::map<int, std::vector<CRDEcalEDM::CRDArborNode*> >& m_orderedNodes,
                            std::vector<CRDEcalEDM::CRDArborTree>& m_treeCol,
                            std::vector<CRDEcalEDM::CRDArborNode*>& m_isoNodes );
  StatusCode MergeConnectedTrees( std::vector<CRDEcalEDM::CRDArborTree>& m_inTreeCol,
                                  std::vector<CRDEcalEDM::CRDArborTree>& m_outTreeCol );

  StatusCode CleanConnection( CRDEcalEDM::CRDArborTree& m_tree );

  StatusCode DepartArborTree( CRDEcalEDM::CRDArborTree& m_tree, std::vector<CRDEcalEDM::CRDArborTree>& m_departedTrees, std::vector<CRDEcalEDM::CRDArborNode*>& m_isoNodes);

  StatusCode MergeNeighborTree( std::vector<CRDEcalEDM::CRDArborTree>& m_inTreeCol, std::vector<CRDEcalEDM::CRDArborTree>& m_outTreeCol );

  //StatusCode ClassifyArborTree( std::vector<CRDEcalEDM::CRDArborTree>& m_inTreeCol, std::vector<CRDEcalEDM::CRDArborTree>& m_GoodTreeCol, std::vector<CRDEcalEDM::CRDArborTree>& m_BadTreeCol );

  //CRDEcalEDM::CRDArborTree GetClosestTree( CRDEcalEDM::CRDArborTree m_badTree, std::vector<CRDEcalEDM::CRDArborTree> m_goodTreeCol );

  Settings settings;

private: 

};

#endif
