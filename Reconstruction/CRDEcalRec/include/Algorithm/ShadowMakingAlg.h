#ifndef _SHADOWMAKING_ALG_H
#define _SHADOWMAKING_ALG_H

#include "PandoraPlusDataCol.h"

using namespace CRDEcalEDM;

class ShadowMakingAlg{
public: 

  class Settings{
  public: 
    Settings() {};
    void SetInitialValue();
    bool fl_UseTrack; 
    int EndLayer; 
    int th_GoodLayer; 

  };

  ShadowMakingAlg(){};
  ~ShadowMakingAlg(){};

  StatusCode Initialize();
  StatusCode RunAlgorithm( ShadowMakingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);
  StatusCode ClearAlgorithm();

  std::vector<CRDEcalEDM::CRDCaloBlock> GetBlocksNeedModification( CRDEcalEDM::CRDCaloHitLongiCluster& m_clus );

  Settings settings; 


};

#endif
