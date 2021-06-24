#ifndef PFO_CREATOR_H
#define PFO_CREATOR_H

#include "PandoraPlusDataCol.h"

#include "TVector3.h"
#include "TLorentzVector.h"

class PFOCreator{

public: 

  class Settings{
  public:
    Settings(){};

  };
  
  //initialize a CaloHitCreator
  PFOCreator( Settings& settings );
  ~PFOCreator() {};

  StatusCode ReadSettings( Settings& settings ){ return StatusCode::SUCCESS; };

  StatusCode CreatePFO(PandoraPlusDataCol& dataCol );

  void Reset(){};

private: 



};
#endif
