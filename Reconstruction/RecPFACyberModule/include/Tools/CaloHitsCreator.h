#ifndef ECALHIT_CREATOR_H
#define ECALHIT_CREATOR_H

#include "k4FWCore/DataHandle.h"
#include "CyberDataCol.h"
#include "Tools/Algorithm.h"
#include "DetInterface/IGeomSvc.h"
#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/Segmentations.h>

namespace Cyber{

  class CaloHitsCreator{

  public: 
    //initialize a CaloHitCreator
    CaloHitsCreator( const Settings& m_settings );
    ~CaloHitsCreator() {};
   
    StatusCode CreateCaloHits( CyberDataCol& m_DataCol,
                               std::vector<DataHandle<edm4hep::CalorimeterHitCollection>*>& r_CaloHitCols );

    //StatusCode CreateMCParticleCaloHitsAsso( std::vector<DataHandle<edm4hep::CalorimeterHitCollection>*>& r_CaloHitCols, 
    //                                         DataHandle<edm4hep::MCRecoCaloParticleAssociationCollection>* r_MCParticleRecoCaloCol );

    //StatusCode Clustering( CyberDataCol& m_DataCol ) { return StatusCode::SUCCESS; };
   
    StatusCode Reset() { return StatusCode::SUCCESS; };

  private: 
    const Cyber::Settings  settings; 

  };
};
#endif
