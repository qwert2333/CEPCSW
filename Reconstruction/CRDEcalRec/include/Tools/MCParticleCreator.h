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

  StatusCode GetMCParticle( PandoraPlusDataCol& dataCol, 
                            const edm4hep::MCParticleCollection& const_MCPCol){

    std::vector<edm4hep::MCParticle> m_MCPcol; m_MCPcol.clear(); 
    for(int imc=0;imc<const_MCPCol.size(); imc++){ 
      edm4hep::MCParticle m_MCp = const_MCPCol[imc];
      //if( m_MCp.daughters_size()!=0 ) continue; 
      m_MCPcol.push_back( const_MCPCol[imc] );

    }
    dataCol.MCParticleCol = m_MCPcol;
    return StatusCode::SUCCESS; 
  };

  StatusCode CreateTrackMCParticleRelation(){ return StatusCode::SUCCESS; };

  StatusCode CreateEcalBarMCParticleRelation(){ return StatusCode::SUCCESS; };

  StatusCode CreateHcalHitsMCParticleRelation(){ return StatusCode::SUCCESS; };


  void Reset(){};

private: 

  



};
#endif
