#ifndef ECALHIT_CREATOR_H
#define ECALHIT_CREATOR_H

#include "k4FWCore/DataHandle.h"
#include "PandoraPlusDataCol.h"
#include "CRDEcalEDMSvc/ICRDEcalEDMSvc.h"

class EcalHitsCreator{

public: 

  class Settings{
  public: 
    Settings(){};

  };
  
  //initialize a CaloHitCreator
  EcalHitsCreator( Settings& settings ){};
  ~EcalHitsCreator() {};

  StatusCode ReadSettings( Settings& settings ){ return StatusCode::SUCCESS; };

  StatusCode GetEcalBars(PandoraPlusDataCol& dataCol , ICRDEcalEDMSvc& m_svc ){ 
    m_svc.getDigiSystem( dataCol.BlockVec );  
    return StatusCode::SUCCESS;
  };

  void Reset(){};

  Settings  settings; 

private: 


};
#endif
