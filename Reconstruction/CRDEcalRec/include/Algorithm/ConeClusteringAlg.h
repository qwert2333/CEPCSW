#ifndef _CONECLUSTERING_ALG_H
#define _CONECLUSTERING_ALG_H
#include "PandoraPlusDataCol.h"

class ConeClusteringAlg {

public: 
  class Settings{
  public: 
    Settings(){};

    void SetInitialValue();
    void SetMergeBadClus(bool flag) { fl_MergeBadClus = flag; }
    void SetGoodClusLevel(int level) { fl_GoodClusLevel = level; }

    double th_ConeTheta_l1;
    double th_ConeR_l1;
    double th_ConeTheta_l2;
    double th_ConeR_l2;
    double th_ClusChi2;
    int  fl_GoodClusLevel;
    bool fl_MergeBadClus; 
  };

  ConeClusteringAlg();

  StatusCode Initialize();

  StatusCode RunAlgorithm( ConeClusteringAlg::Settings& settings, PandoraPlusDataCol& m_datacol);

  int CheckClusterQuality(CRDEcalEDM::CRDCaloHit3DShower& clus);

  bool MergeToGoodCluster( std::vector<CRDEcalEDM::CRDCaloHit3DShower>& goodClusCol,  CRDEcalEDM::CRDCaloHit3DShower& badClus, bool ForceMerging);

  CRDEcalEDM::CRDCaloHit3DShower GetClosestGoodCluster( std::vector<CRDEcalEDM::CRDCaloHit3DShower>& goodClusCol,  CRDEcalEDM::CRDCaloHit3DShower& badClus );

  Settings settings;

private: 


};
#endif
