#ifndef HCALHIT_CREATOR_H
#define HCALHIT_CREATOR_H

#include "PandoraPlusDataCol.h"


class HcalHitsCreator{

public: 

  class Settings{
  public:
    Settings(){};

  };
  
  //initialize a CaloHitCreator
  HcalHitsCreator( Settings& settings ){};
  ~HcalHitsCreator() {};

  StatusCode ReadSettings( Settings& settings ){ return StatusCode::SUCCESS; };

  StatusCode GetHcalHits(PandoraPlusDataCol& dataCol ){ return StatusCode::SUCCESS; };



  void Reset(){};

private: 

  



};
#endif
