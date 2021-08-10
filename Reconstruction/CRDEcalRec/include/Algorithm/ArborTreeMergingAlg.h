#ifndef _ARBORTREEMERGING_ALG_H
#define _ARBORTREEMERGING_ALG_H

#include "PandoraPlusDataCol.h"
using namespace CRDEcalEDM; 

class ArborTreeMergingAlg{

public: 

  class Settings{
  public:
    Settings() {}; 
    ~Settings() {};

    void SetInitialValue(); 
    void SetGoodTreeLevel(int _tree, int _node) { th_GoodTreeLevel=_tree; th_GoodTreeNodes=_node; }
    void SetMergeTrees(bool _fl) { fl_MergeTrees=_fl; }

    int th_Nsigma = 2;
    double th_daughterR = 40.;
    double th_MergeR = 60;
    double th_MergeTheta = PI/18.;
    int th_GoodTreeLevel;
    int th_GoodTreeNodes;
    bool fl_MergeTrees = true; 
  };

  ArborTreeMergingAlg() {};
  ~ArborTreeMergingAlg() {};

  StatusCode Initialize();
  StatusCode RunAlgorithm( ArborTreeMergingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);
  StatusCode ClearAlgorithm();

  StatusCode GetSecondaryVtxT1(  std::vector<CRDEcalEDM::CRDArborTree>& m_ArborTreeCol,
                                 std::vector<CRDEcalEDM::CRDCaloHit2DShower>& m_2DshowerCol,
                                 std::map<CRDEcalEDM::CRDArborNode*, CRDEcalEDM::CRDArborTree*>& m_svtx );

  StatusCode GetSecondaryVtxT2(  std::vector<CRDEcalEDM::CRDArborTree>& m_ArborTreeCol,
                                 std::map<CRDEcalEDM::CRDArborNode*, CRDEcalEDM::CRDArborTree*>& m_svtx );


  Settings settings; 

};
#endif
