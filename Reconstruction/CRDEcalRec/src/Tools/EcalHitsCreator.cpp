#ifndef TRACK_CREATOR_C
#define TRACK_CREATOR_C

#include "Tools/EcalHitsCreator.h"

EcalHitsCreator::EcalHitsCreator(const Settings& m_settings) : settings( m_settings ){

} 


StatusCode EcalHitsCreator::CreateEcalHits( PandoraPlusDataCol& m_DataCol ){

  if(settings.m_EcalCaloHitCollections.size()==0 || settings.m_EcalReadouts.size()==0) return StatusCode::SUCCESS;

  m_DataCol.collectionMap_CaloHit.clear(); 
  for(unsigned int icol=0; icol<settings.m_EcalCaloHitCollections.size(); icol++){
    std::string col_name = settings.m_EcalCaloHitCollections[icol];

    DataHandle<edm4hep::CalorimeterHit> r_CaloHitCol{col_name, Gaudi::DataHandle::Reader, this};
    const edm4hep::CalorimeterHitCollection* const_CaloHitCol = r_CaloHitCol.get(); 

    std::vector<edm4hep::CalorimeterHit> m_HitCol; m_HitCol.clear(); 
    for(unsigned int ihit=0; ihit<const_CaloHitCol->size(); ihit++){
      edm4hep::CalorimeterHit m_hit = const_CaloHitCol->at(ihit);
      m_HitCol.push_back(m_hit);
    }

    m_DataCol.collectionMap_CaloHit[col_name] = m_HitCol; 
  }


  return StatusCode::SUCCESS;
}


#endif
