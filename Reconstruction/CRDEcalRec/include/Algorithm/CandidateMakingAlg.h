#ifndef _CANDIDATEMAKING_ALG_H
#define _CANDIDATEMAKING_ALG_H

#include "PandoraPlusDataCol.h"
#include "Algorithm/EnergySplittingAlg.h"
#include "Algorithm/EnergyTimeMatchingAlg.h"
#include "Algorithm/ConeClusteringAlg.h"

using namespace CRDEcalEDM;

class CandidateMakingAlg{

public: 
  class Settings{
  public: 
    Settings() {};
    void SetInitialValue();
    void SetDebug(int _level)  { Debug  =_level; }
    void SetUseTrk( bool _fl ) { UseTrk = _fl;   }

    int  Debug;
    bool UseTrk; 
  };

  CandidateMakingAlg() {};
  StatusCode Initialize();
  StatusCode RunAlgorithm( CandidateMakingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);
  StatusCode ClearAlgorithm(); 

  std::vector<CRDEcalEDM::CRDCaloBlock> GetBlocksNeedModification( std::vector<CRDEcalEDM::CRDCaloHit3DShower>& m_goodClus, std::vector<CRDEcalEDM::CRDCaloHit3DShower>& m_badClus );

  Settings settings;


};

#endif
