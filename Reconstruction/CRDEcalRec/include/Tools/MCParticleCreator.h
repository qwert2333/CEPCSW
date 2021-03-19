#ifndef MCPARTICLE_CREATOR_H
#define MCPARTICLE_CREATOR_H

#include "PandoraPlusDataCol.h"


class MCParticleCreator{

public: 

  class Settings{
  public:
    Settings(){};

  };
  
  //initialize a CaloHitCreator
  MCParticleCreator( Settings& settings ){};
  ~MCParticleCreator();

  StatusCode ReadSettings( Settings& settings ){ return StatusCode::SUCCESS; };

  StatusCode GetMCParticle(PandoraPlusDataCol& dataCol ){ return StatusCode::SUCCESS; };

  StatusCode CreateTrackMCParticleRelation(){ return StatusCode::SUCCESS; };

  StatusCode CreateEcalBarMCParticleRelation(){ return StatusCode::SUCCESS; };

  StatusCode CreateHcalHitsMCParticleRelation(){ return StatusCode::SUCCESS; };


  void Reset(){};

private: 

  



};
#endif
