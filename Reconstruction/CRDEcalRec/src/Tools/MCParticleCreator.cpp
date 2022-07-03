#ifndef MCPARTICLE_CREATOR_C
#define MCPARTICLE_CREATOR_C

#include "Tools/MCParticleCreator.h"

namespace PandoraPlus{
  MCParticleCreator::MCParticleCreator( const Settings& m_settings ): settings(m_settings){

  }


  StatusCode MCParticleCreator::CreateMCParticle( PandoraPlusDataCol& m_DataCol, DataHandle<edm4hep::MCParticleCollection>& r_MCParticleCol ){
    if(settings.map_stringPars.at("MCParticleCollections").empty()) return StatusCode::SUCCESS;
    m_DataCol.collectionMap_MC.clear();

    const edm4hep::MCParticleCollection* const_MCPCol = r_MCParticleCol.get();

    std::vector<edm4hep::MCParticle> m_MCPvec; m_MCPvec.clear(); 
    for(int imc=0; imc<const_MCPCol->size(); imc++){
      edm4hep::MCParticle m_MCp = const_MCPCol->at(imc);
      m_MCPvec.push_back(m_MCp);
    }

    m_DataCol.collectionMap_MC[ settings.map_stringPars.at("MCParticleCollections") ] = m_MCPvec; 

    return StatusCode::SUCCESS;
  };

}
#endif
