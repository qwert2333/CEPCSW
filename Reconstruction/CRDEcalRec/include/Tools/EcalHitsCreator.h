#ifndef ECALHIT_CREATOR_H
#define ECALHIT_CREATOR_H

#include "k4FWCore/DataHandle.h"
#include "PandoraPlusDataCol.h"
#include "CRDEcalSvc/ICRDEcalSvc.h"

class EcalHitsCreator{

public: 

  class Settings{
  public: 
    Settings(){};

    std::vector<std::string> m_EcalCaloHitCollections; 
    std::vector<std::string> m_EcalReadouts; 

  };
  
  //initialize a CaloHitCreator
  EcalHitsCreator( const Settings& m_settings );
  ~EcalHitsCreator() {};

  StatusCode CreateEcalHits( PandoraPlusDataCol& m_DataCol ); 

  StatusCode Reset() { return StatusCode::SUCCESS; };

private: 
  const Settings  settings; 


};
#endif
