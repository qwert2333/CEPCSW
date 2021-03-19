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
  ~EcalHitsCreator();

  StatusCode ReadSettings( Settings& settings ){ return StatusCode::SUCCESS; };

  StatusCode GetEcalBars(PandoraPlusDataCol& dataCol , ICRDEcalEDMSvc& m_svc ){ 
    dataCol.blockVec = m_svc.getDigiSystem();  
    return StatusCode::SUCCESS;
  };

  void Reset(){};

private: 


};
#endif
