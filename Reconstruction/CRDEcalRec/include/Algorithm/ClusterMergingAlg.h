#ifndef _CLUSTERMERGING_ALG_H
#define _CLUSTERMERGING_ALG_H

#include "PandoraPlusDataCol.h"

using namespace CRDEcalEDM;

class ClusterMergingAlg{

public: 
  class Settings{
  public: 
    Settings() {};
    void SetInitialValue();
    void SetValue(double _dR1, double _dR2, int _ly) { axis_Angle=_dR1; relP_Angle=_dR2; skipLayer=_ly; }
    void SetDebug(int _level)  { Debug  =_level; }
    void SetMergeGoodCluster( bool _fl ) { fl_MergeGoodClus=_fl; }
    void SetMergeBadCluster( bool _fl ) { fl_MergeBadClus=_fl; }
    void SetMergeEMTailShower( bool _fl ) { fl_MergeEMTailShower=_fl; }
    void SetMergeEMTailCluster( bool _fl ) { fl_MergeEMTailCluster=_fl; }

    double axis_Angle;  //Delta R between axis of two cluster
    double relP_Angle;  //Delta R between axis(cluster1) and r(cluster1-cluster2). 
    int skipLayer;      //Allowed layer between two cluster. 
    bool fl_MergeGoodClus ;
    bool fl_MergeBadClus ;
    bool fl_MergeEMTailShower; 
    bool fl_MergeEMTailCluster; 
    int  Debug;
  };

  ClusterMergingAlg() {};
  ~ClusterMergingAlg() {};
  StatusCode Initialize();
  StatusCode RunAlgorithm( ClusterMergingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);
  StatusCode ClearAlgorithm(); 

  bool MergeToGoodCluster( std::vector<CRDEcalEDM::CRDCaloHit3DCluster>& goodClusCol,  CRDEcalEDM::CRDCaloHit3DCluster& badClus, bool ForceMerging);

  CRDEcalEDM::CRDCaloHit3DCluster GetClosestGoodCluster( std::vector<CRDEcalEDM::CRDCaloHit3DCluster>& goodClusCol,  CRDEcalEDM::CRDCaloHit3DCluster& badClus );

  static bool compBegin( CRDEcalEDM::CRDCaloHit3DCluster& clus1, CRDEcalEDM::CRDCaloHit3DCluster& clus2 ) { return clus1.getBeginningDlayer()<clus2.getBeginningDlayer(); }

  Settings settings;

private: 

};

#endif
