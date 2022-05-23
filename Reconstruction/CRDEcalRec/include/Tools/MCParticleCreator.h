#ifndef MCPARTICLE_CREATOR_H
#define MCPARTICLE_CREATOR_H

#include "PandoraPlusDataCol.h"


class MCParticleCreator{

public: 

  class Settings{
  public:
    Settings(){};

    std::string m_mcParticleCollections;        // MC particle collection. 
    std::string m_CaloHitRelationCollections;   // SimCaloHit to CaloHit particle relations
    float m_BField; 

  };
  
  //initialize a CaloHitCreator
  MCParticleCreator( const Settings& m_settings );
  ~MCParticleCreator() {};

  StatusCode CreateMCParticles( PandoraPlusDataCol& m_DataCol ); 

  //StatusCode CreateTrackMCParticleRelation(){ return StatusCode::SUCCESS; };

  //StatusCode CreateEcalBarMCParticleRelation(){ return StatusCode::SUCCESS; };

  //StatusCode CreateHcalHitsMCParticleRelation(){ return StatusCode::SUCCESS; };

  StatusCode Reset() { return StatusCode::SUCCESS; };


private: 

  const Settings   settings;  



};
#endif
