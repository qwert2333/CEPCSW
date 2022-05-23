#ifndef MCPARTICLE_CREATOR_C
#define MCPARTICLE_CREATOR_C

#include "k4FWCore/DataHandle.h"
#include "Tools/MCParticleCreator.h"

MCParticleCreator::MCParticleCreator( const Settings& m_settings ): settings(m_settings){

}


StatusCode MCParticleCreator::CreateMCParticles( PandoraPlusDataCol& m_DataCol ){
  if(settings.m_mcParticleCollections.empty()) return StatusCode::SUCCESS;
  m_DataCol.collectionMap_MC.clear();

  //std::string col_name = settings.m_mcParticleCollections; 

  //DataHandle<edm4hep::MCParticleCollection>* r_MCParticleCol = new DataHandle<edm4hep::MCParticleCollection>(settings.m_mcParticleCollections, Gaudi::DataHandle::Reader, this);
  //DataHandle<edm4hep::MCParticleCollection> r_MCParticleCol{"MCParticle", Gaudi::DataHandle::Reader, this};
  const edm4hep::MCParticleCollection* const_MCPCol = r_MCParticleCol.get();
  auto const_MCPCol = r_mcParticle.get();


  std::vector<edm4hep::MCParticle> m_MCPvec; m_MCPvec.clear(); 
  for(int imc=0; imc<const_MCPCol->size(); imc++){
    edm4hep::MCParticle m_MCp = const_MCPCol->at(imc);
    m_MCPvec.push_back(m_MCp);
  }

  m_DataCol.collectionMap_MC[settings.m_mcParticleCollections] = m_MCPvec; 

  return StatusCode::SUCCESS;
}



#endif
