#ifndef _CONECLUSTERING_ALG_H
#define _CONECLUSTERING_ALG_H
#include "PandoraPlusDataCol.h"

class ConeClusteringAlg {

public: 
  class Settings{
  public: 
    Settings(){};
    void SetInitialValue();
    double th_ConeTheta;
    double th_ConeR;
    double th_ClusChi2;
  };

  ConeClusteringAlg();

  StatusCode Initialize();

  StatusCode RunAlgorithm( ConeClusteringAlg::Settings& settings, PandoraPlusDataCol& m_datacol);

  bool CheckClusterQuality(CRDEcalEDM::CRDCaloHit3DShower& clus);

  bool MergeToGoodCluster( std::vector<CRDEcalEDM::CRDCaloHit3DShower>& goodClusCol,  CRDEcalEDM::CRDCaloHit3DShower& badClus, bool ForceMerging);

  CRDEcalEDM::CRDCaloHit3DShower GetClosestGoodCluster( std::vector<CRDEcalEDM::CRDCaloHit3DShower>& goodClusCol,  CRDEcalEDM::CRDCaloHit3DShower& badClus );

  Settings settings;

private: 


};
#endif
