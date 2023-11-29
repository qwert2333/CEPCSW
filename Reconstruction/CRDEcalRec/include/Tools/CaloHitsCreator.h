#ifndef ECALHIT_CREATOR_H
#define ECALHIT_CREATOR_H

#include "k4FWCore/DataHandle.h"
#include "PandoraPlusDataCol.h"
#include "Tools/Algorithm.h"
#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/Segmentations.h>

namespace PandoraPlus{

  class CaloHitsCreator{

  public: 
    //initialize a CaloHitCreator
    CaloHitsCreator( const Settings& m_settings );
    ~CaloHitsCreator() {};
   
    StatusCode CreateCaloHits( PandoraPlusDataCol& m_DataCol,
                               std::vector<DataHandle<edm4hep::CalorimeterHitCollection>*>& r_CaloHitCols, 
                               std::map<std::string, dd4hep::DDSegmentation::BitFieldCoder*>& map_decoder, 
                               std::map<std::string, DataHandle<edm4hep::MCRecoCaloParticleAssociationCollection>*>& map_CaloParticleAssoCol ); 

    //StatusCode CreateMCParticleCaloHitsAsso( std::vector<DataHandle<edm4hep::CalorimeterHitCollection>*>& r_CaloHitCols, 
    //                                         DataHandle<edm4hep::MCRecoCaloParticleAssociationCollection>* r_MCParticleRecoCaloCol );

    //StatusCode Clustering( PandoraPlusDataCol& m_DataCol ) { return StatusCode::SUCCESS; };
   
    StatusCode Reset() { return StatusCode::SUCCESS; };

  private: 
    const PandoraPlus::Settings  settings; 

  };
};
#endif
