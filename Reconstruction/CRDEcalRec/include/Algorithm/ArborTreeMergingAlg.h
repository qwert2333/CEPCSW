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
    void SetGoodTreeLevel(int _ly1, int _ly2, int _node) { th_GoodTreeLayer1=_ly1; th_GoodTreeLayer2=_ly2; th_GoodTreeNodes=_node; }
    void SetMergeTrees(bool _fl) { fl_MergeTrees=_fl; }
    void SetMergeTreePars(double _Rth, double _Tth) { th_MergeR=_Rth; th_MergeTheta=_Tth; }
    void SetOverwrite(bool _fl) { fl_overwrite=_fl; }
    void SetClusterType(string _st) { clusType=_st; }

    int th_Nsigma = 2;
    double th_daughterR = 40.;
    double th_MergeR = 60;
    double th_MergeTheta = PI/18.;
    int th_GoodTreeLayer1;
    int th_GoodTreeLayer2;
    int th_GoodTreeNodes;    //GoodTree criteria: Nlayer>=th_GoodTreeLayer1 || (Nlayer>=th_GoodTreeLayer2 && Nnodes>=th_GoodTreeNodes). 
    bool fl_MergeTrees = true; 
    bool fl_overwrite = true; 

    string clusType = "";
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
