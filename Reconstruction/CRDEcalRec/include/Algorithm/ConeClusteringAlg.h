#ifndef _CONECLUSTERING_ALG_H
#define _CONECLUSTERING_ALG_H
#include "PandoraPlusDataCol.h"

class ConeClusteringAlg {

public: 
  class Settings{
  public: 
    Settings(){};

    void SetInitialValue();
    void SetConeValue( double _coneTheta_l1, double _coneR_l1, double _coneTheta_l2, double _coneR_l2 );
    void SetGoodClusLevel(int level) { fl_GoodClusLevel = level; }
    void SetUseCandidate(int level) { fl_UseCandidate = level; }


    double th_ConeTheta_l1;
    double th_ConeR_l1;
    double th_ConeTheta_l2;
    double th_ConeR_l2;
    double th_ClusChi2;
    int  fl_GoodClusLevel;
    int  fl_UseCandidate; 
  };

  ConeClusteringAlg();

  StatusCode Initialize();

  StatusCode RunAlgorithm( ConeClusteringAlg::Settings& settings, PandoraPlusDataCol& m_datacol);

  StatusCode LongiConeLinking( std::map<int, std::vector<CRDEcalEDM::CRDCaloHit2DShower> >& orderedShower, std::vector<CRDEcalEDM::CRDCaloHit3DShower>& ClusterCol );

  int CheckClusterQuality(CRDEcalEDM::CRDCaloHit3DShower& clus);


  Settings settings;

private: 


};
#endif
