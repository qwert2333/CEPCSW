#ifndef ECALHIT_CREATOR_H
#define ECALHIT_CREATOR_H

#include "k4FWCore/DataHandle.h"
#include "PandoraPlusDataCol.h"
#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/Segmentations.h>

namespace PandoraPlus{

  class CaloHitsCreator{

  public: 

    class Settings{
    public: 
      Settings(){};
   
      std::vector<std::string> m_CaloHitCollections; 
      //std::vector<std::string> m_CaloReadouts; 
   
    };
    
    //initialize a CaloHitCreator
    CaloHitsCreator( const Settings& m_settings );
    ~CaloHitsCreator() {};
   
    StatusCode CreateCaloHits( PandoraPlusDataCol& m_DataCol,
                               std::vector<DataHandle<edm4hep::CalorimeterHitCollection>*>& r_ECalHitCols, 
                               std::map<std::string, dd4hep::DDSegmentation::BitFieldCoder*>& map_decoder  ); 
   
    StatusCode Clustering( PandoraPlusDataCol& m_DataCol ) { return StatusCode::SUCCESS; };
   
    StatusCode Reset() { return StatusCode::SUCCESS; };

  private: 
    const Settings  settings; 

  };
};
#endif
