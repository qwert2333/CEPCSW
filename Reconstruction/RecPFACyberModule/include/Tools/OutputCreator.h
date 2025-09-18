#ifndef OUTPUT_CREATOR_H
#define OUTPUT_CREATOR_H

#include "k4FWCore/DataHandle.h"
#include "edm4hep/MutableCalorimeterHit.h"
#include "edm4hep/Vector3f.h"
#include "CyberDataCol.h"
#include "Tools/Algorithm.h"
//#include "Tools/HelixClassD.h"

namespace Cyber{

  class OutputCreator{
  public: 

    OutputCreator( const Settings& m_settings);
    ~OutputCreator() {};

    StatusCode CreateOutputCollections( CyberDataCol& m_DataCol, 
                                        DataHandle<edm4hep::CalorimeterHitCollection>& m_outRecHitsHandler, 
                                        DataHandle<edm4hep::CalorimeterHitCollection>& m_outRecCoreHandler, 
                                        DataHandle<edm4hep::CalorimeterHitCollection>& m_outRecHcalHitsHandler, 
                                        DataHandle<edm4hep::TrackCollection>& m_outTrkHandler, 
                                        std::map<std::string, DataHandle<edm4hep::ClusterCollection>*>& m_outClusterColHandler, 
                                        DataHandle<edm4hep::ReconstructedParticleCollection>& m_recPFOHandler );


    StatusCode Reset() { return StatusCode::SUCCESS; }

    edm4hep::Track TruthTrack(edm4hep::MCParticle _mcp, edm4hep::TrackCollection* _trkCol );

  private: 
    const Cyber::Settings   settings;
  
  };
};
#endif
