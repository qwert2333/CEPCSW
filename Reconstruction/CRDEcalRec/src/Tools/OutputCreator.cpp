#ifndef OUTPUT_CREATOR_C
#define OUTPUT_CREATOR_C

#include "Tools/OutputCreator.h"

namespace PandoraPlus{

  OutputCreator::OutputCreator( const Settings& m_settings ): settings(m_settings){

  }

  StatusCode OutputCreator::CreateRecCaloHits( PandoraPlusDataCol& m_DataCol, DataHandle<edm4hep::CalorimeterHitCollection>& m_outRecHitsHandler ){
    edm4hep::CalorimeterHitCollection* m_calohitCol = m_outRecHitsHandler.createAndPut();

    std::vector<PandoraPlus::CaloCluster*> p_clusCol = m_DataCol.map_CaloCluster["EcalCluster"];
    for(int ic=0; ic<p_clusCol.size(); ic++){
      std::vector<const PandoraPlus::CaloHit*> p_hits = p_clusCol[ic]->getCaloHits();
      for(int ih=0; ih<p_hits.size(); ih++){
        auto _hit = m_calohitCol->create();
        _hit.setCellID(0);
        _hit.setEnergy( p_hits[ih]->getEnergy() );
        _hit.setPosition( edm4hep::Vector3f(p_hits[ih]->getPosition().x(), p_hits[ih]->getPosition().y(), p_hits[ih]->getPosition().z()) );
      }
    }
    return StatusCode::SUCCESS;
  }

  StatusCode OutputCreator::CreateCluster( PandoraPlusDataCol& m_DataCol, DataHandle<edm4hep::ClusterCollection>& m_outClusterColHandler ){
    edm4hep::ClusterCollection* m_clusCol = m_outClusterColHandler.createAndPut();

    std::vector<PandoraPlus::CaloCluster*> p_clusCol = m_DataCol.map_CaloCluster["EcalCluster"]; 
    for(int ic=0; ic<p_clusCol.size(); ic++){
      auto _clus = m_clusCol->create();

      for(int ih=0; ih<p_clusCol[ic]->getCaloHits().size(); ih++){
        edm4hep::MutableCalorimeterHit _hit;
        _hit.setCellID(0); //NOTE: Need a cellID coder if want to pass CaloHit to Pandora/Arbor. 
        _hit.setEnergy( p_clusCol[ic]->getCaloHits()[ih]->getEnergy() );
        edm4hep::Vector3f pos(p_clusCol[ic]->getCaloHits()[ih]->getPosition().x(), p_clusCol[ic]->getCaloHits()[ih]->getPosition().y(), p_clusCol[ic]->getCaloHits()[ih]->getPosition().z());
        _hit.setPosition( pos );
        //_hit.setType();
        _clus.addToHits(_hit);
      }
      _clus.setEnergy( p_clusCol[ic]->getShowerE() );
      edm4hep::Vector3f pos( p_clusCol[ic]->getShowerCenter().x(), p_clusCol[ic]->getShowerCenter().y(), p_clusCol[ic]->getShowerCenter().z() );
      _clus.setPosition( pos );
    }

    return StatusCode::SUCCESS;
  }

  StatusCode OutputCreator::CreatePFO( PandoraPlusDataCol& m_DataCol, DataHandle<edm4hep::ReconstructedParticleCollection>& m_recPFOHandler ){
    //edm4hep::ReconstructedParticleCollection* m_pfoCol = m_recPFOHandler.createAndPut();


    return StatusCode::SUCCESS;
  }

}
#endif
