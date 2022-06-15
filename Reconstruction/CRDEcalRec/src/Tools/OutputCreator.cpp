#ifndef OUTPUT_CREATOR_C
#define OUTPUT_CREATOR_C

#include "Tools/OutputCreator.h"

namespace PandoraPlus{

  OutputCreator::OutputCreator( const Settings& m_settings ): settings(m_settings){

  }

  StatusCode OutputCreator::CreateRecCaloHits( PandoraPlusDataCol& m_DataCol, DataHandle<edm4hep::CalorimeterHitCollection>& m_outRecHitsHandler ){

    return StatusCode::SUCCESS;
  }

  StatusCode OutputCreator::CreateCluster( PandoraPlusDataCol& m_DataCol, DataHandle<edm4hep::ClusterCollection>& m_outClusterColHandler ){
    edm4hep::ClusterCollection* m_clusCol = m_outClusterColHandler.createAndPut();

    std::vector<PandoraPlus::CaloCluster*> p_clusCol = m_DataCol.ClusterCol; 
    for(int ic=0; ic<p_clusCol.size(); ic++){
      auto _clus = m_clusCol->create();

      for(int ih=0; ih<p_clusCol[ic]->getCaloHits().size(); ih++){
        //MutableCalorimeterHit _hit;
        //_hit.setCellID(0);
        //_hit.setEnergy( p_clusCol[ic]->getCaloHits()[ih]->getEnergy() );
        //_hit.setPosition(  p_clusCol[ic]->getCaloHits()[ih]->getPosition() );
        //_hit.setType();
        //_clus.addToHits(_hit);
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
