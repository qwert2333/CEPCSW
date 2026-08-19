#ifndef OUTPUT_CREATOR_H
#define OUTPUT_CREATOR_H

#include "k4FWCore/DataHandle.h"
#include "edm4hep/MutableCalorimeterHit.h"
#include "edm4hep/Vector3f.h"
#include "PandoraPlusDataCol.h"
#include "Tools/Algorithm.h"

namespace PandoraPlus{

  class OutputCreator{
  public: 

    OutputCreator( const Settings& m_settings);
    ~OutputCreator() {};

    StatusCode CreateRecCaloHits( PandoraPlusDataCol& m_DataCol,  DataHandle<edm4hep::CalorimeterHitCollection>& m_outRecHitsHandler); 

    StatusCode CreateCluster( PandoraPlusDataCol& m_DataCol, DataHandle<edm4hep::ClusterCollection>& m_outClusterColHandler );

    StatusCode CreatePFO( PandoraPlusDataCol& m_DataCol, DataHandle<edm4hep::ReconstructedParticleCollection>& m_recPFOHandler );    

    StatusCode Reset() { return StatusCode::SUCCESS; }

  private: 
    const PandoraPlus::Settings   settings;
  
  };
};
#endif
