#ifndef MCPARTICLE_CREATOR_H
#define MCPARTICLE_CREATOR_H

#include "PandoraPlusDataCol.h"

namespace PandoraPlus{

  class MCParticleCreator{

  public: 

    class Settings{
    public:
      Settings(){};

      std::string m_mcParticleCollections;        // MC particle collection. 
      std::string m_CaloHitRelationCollections;   // SimCaloHit to CaloHit particle relations

    };
    
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

    const Settings   settings;  

  };

};
#endif
