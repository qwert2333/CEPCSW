#ifndef TRACK_CREATOR_C
#define TRACK_CREATOR_C

#include "Tools/CaloHitsCreator.h"

namespace Cyber{
  CaloHitsCreator::CaloHitsCreator(const Settings& m_settings) : settings( m_settings ){

  }; 

  StatusCode CaloHitsCreator::CreateCaloHits( CyberDataCol& m_DataCol, 
                                              std::vector<DataHandle<edm4hep::CalorimeterHitCollection>*>& r_CaloHitCols)
//                                              std::map<std::string, dd4hep::DDSegmentation::BitFieldCoder*>& map_decoder,
//                                              std::map<std::string, DataHandle<edm4hep::MCRecoCaloParticleAssociationCollection>*>& map_CaloParticleAssoCol,
//                                              SmartIF<IGeomSvc>& m_geosvc,
//                                              std::map<std::tuple<int, int, int, int, int>, int>& barNumberMapEndcapMap )
  {
    if(r_CaloHitCols.size()==0 || settings.map_stringVecPars.at("CaloHitCollections").size()==0) StatusCode::SUCCESS;

    //Save readin collections
    m_DataCol.collectionMap_CaloHit.clear(); 
    for(unsigned int icol=0; icol<r_CaloHitCols.size(); icol++){

      const edm4hep::CalorimeterHitCollection* const_CaloHitCol = r_CaloHitCols[icol]->get(); 

      std::vector<edm4hep::CalorimeterHit> m_HitCol; m_HitCol.clear(); 
      for(unsigned int ihit=0; ihit<const_CaloHitCol->size(); ihit++){
        edm4hep::CalorimeterHit m_hit = const_CaloHitCol->at(ihit);
        m_HitCol.push_back(m_hit);
      }
cout<<"  Hit size: "<<m_HitCol.size()<<endl;

      m_DataCol.collectionMap_CaloHit[settings.map_stringVecPars.at("CaloHitCollections")[icol]] = m_HitCol; 
    }

    //Convert to local objects: CalorimeterHit to CaloUnit (For ECALBarrel only)
    if(settings.map_stringPars.at("EcalType")=="BarEcal"){
      std::vector<std::shared_ptr<Cyber::CaloUnit>> m_barCol; m_barCol.clear(); 

      //Readin ECAL barrel hits      
      if( m_DataCol.collectionMap_CaloHit.find("CalibCaloHit") != m_DataCol.collectionMap_CaloHit.end() ){
        auto CaloHits = m_DataCol.collectionMap_CaloHit["CalibCaloHit"]; 
        std::map<int, std::vector<Cyber::CaloUnit>> map_cellID_hits; map_cellID_hits.clear();
        for(auto& hit : CaloHits){ 
          Cyber::CaloUnit m_bar; 
          m_bar.setcellID(hit.getCellID());
          m_bar.setPosition( TVector3(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z) );
          m_bar.setQ(hit.getEnergy()/2., hit.getEnergy()/2.);
          m_bar.setT(hit.getTime(), hit.getTime());

          //unsigned long long tmp_id = hit.getCellID();
          //dd4hep::PlacedVolume ipv = m_volumeManager.lookupVolumePlacement(tmp_id);
          //dd4hep::Volume ivol = ipv.volume();
          //std::vector< double > iVolParam = ivol.solid().dimensions();
          //auto maxElement = std::max_element(iVolParam.begin(), iVolParam.end());
          //iVolParam.clear();
          //m_bar.setBarLength(*maxElement * 20);
          m_bar.setBarLength(120.);

          map_cellID_hits[ hit.getCellID()/100 ].push_back(m_bar);
        }
        for(auto& hit : map_cellID_hits){
          if(hit.second.size()!=2){ std::cout<<"WARNING: didn't find correct hit pairs! "<<std::endl; continue; }
   
          //Cyber::CaloUnit* m_bar = new Cyber::CaloUnit(); 
          std::shared_ptr<Cyber::CaloUnit> m_bar = std::make_shared<Cyber::CaloUnit>();
   
          unsigned long long id = hit.second[0].getcellID();
          int layer = id/10000;
          int index = (id%10000)/100;
          int side = id%100;
          m_bar->setcellID( id );
          m_bar->setcellID( Cyber::CaloUnit::System_Barrel,
                            25,
                            7,
                            -1,               
                            (layer+1)/2,
                            layer%2,
                            index);
          m_bar->setPosition(hit.second[0].getPosition());
          m_bar->setBarLength(hit.second[0].getBarLength());  
          m_bar->setQ( hit.second[0].getEnergy(), hit.second[1].getEnergy() );
          m_bar->setT( hit.second[0].getT1(), hit.second[1].getT1() );

          m_barCol.push_back(m_bar);  //Save for later use in algorithms
        }
   
        m_DataCol.map_BarCol["BarCol"] = m_barCol; 
      }
std::cout<<"  After calohit creactor: bar size "<<m_barCol.size()<<std::endl;
      m_barCol.clear();

    }
    return StatusCode::SUCCESS;
  };

};

#endif
