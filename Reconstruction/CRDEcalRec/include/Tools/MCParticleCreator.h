#ifndef MCPARTICLE_CREATOR_H
#define MCPARTICLE_CREATOR_H

#include "PandoraPlusDataCol.h"
#include "Tools/Algorithm.h"

namespace PandoraPlus{

  class MCParticleCreator{

  public: 

    //initialize a CaloHitCreator
    MCParticleCreator( const Settings& m_settings );
    ~MCParticleCreator() {};

    StatusCode CreateMCParticle( PandoraPlusDataCol& m_DataCol, 
                                 DataHandle<edm4hep::MCParticleCollection>& r_MCParticleCol ); 

    //StatusCode CreateTrackMCParticleRelation(){ return StatusCode::SUCCESS; };

    //StatusCode CreateEcalBarMCParticleRelation(){ return StatusCode::SUCCESS; };

    //StatusCode CreateHcalHitsMCParticleRelation(){ return StatusCode::SUCCESS; };

    StatusCode Reset() { return StatusCode::SUCCESS; };

  private: 
    const PandoraPlus::Settings   settings;  

  };

};
#endif
