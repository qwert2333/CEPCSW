#ifndef TRACK_CREATOR_C
#define TRACK_CREATOR_C

#include "Tools/CaloHitsCreator.h"

namespace PandoraPlus{
  CaloHitsCreator::CaloHitsCreator(const Settings& m_settings) : settings( m_settings ){

  }; 

  StatusCode CaloHitsCreator::CreateCaloHits( PandoraPlusDataCol& m_DataCol, 
                                              std::vector<DataHandle<edm4hep::CalorimeterHitCollection>*>& r_CaloHitCols, 
                                              std::map<std::string, dd4hep::DDSegmentation::BitFieldCoder*>& map_decoder )
  {
    if(r_CaloHitCols.size()==0 || settings.m_CaloHitCollections.size()==0) StatusCode::SUCCESS;

    //Save readin collections
    m_DataCol.collectionMap_CaloHit.clear(); 
    for(unsigned int icol=0; icol<r_CaloHitCols.size(); icol++){

      const edm4hep::CalorimeterHitCollection* const_CaloHitCol = r_CaloHitCols[icol]->get(); 

      std::vector<edm4hep::CalorimeterHit> m_HitCol; m_HitCol.clear(); 
      for(unsigned int ihit=0; ihit<const_CaloHitCol->size(); ihit++){
        edm4hep::CalorimeterHit m_hit = const_CaloHitCol->at(ihit);
        m_HitCol.push_back(m_hit);
      }

      m_DataCol.collectionMap_CaloHit[settings.m_CaloHitCollections[icol]] = m_HitCol; 
    }

    //Convert to local objects: CalorimeterHit to CaloBar (For ECALBarrel only)
    std::vector<PandoraPlus::CaloBar*> m_barCol; m_barCol.clear(); 

    auto CaloHits = m_DataCol.collectionMap_CaloHit["ECALBarrel"]; 
    std::map<std::uint64_t, std::vector<PandoraPlus::CaloBar> > map_cellID_hits; map_cellID_hits.clear();
    for(auto& hit : CaloHits){ 
      PandoraPlus::CaloBar m_bar; 
      m_bar.setcellID(hit.getCellID());
      m_bar.setPosition( TVector3(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z) );
      m_bar.setQ(hit.getEnergy(), hit.getEnergy());
      m_bar.setT(hit.getTime(), hit.getTime());
      map_cellID_hits[hit.getCellID()].push_back(m_bar);
    }

    for(auto& hit : map_cellID_hits){
      if(hit.second.size()!=2){ std::cout<<"WARNING: didn't find correct hit pairs! "<<std::endl; continue; }

      PandoraPlus::CaloBar* m_bar = new PandoraPlus::CaloBar(); 
      unsigned long long id = hit.first; 
      m_bar->setcellID( id );
      m_bar->setcellID( map_decoder["ECALBarrel"]->get(id, "system"),
                        map_decoder["ECALBarrel"]->get(id, "module"),
                        map_decoder["ECALBarrel"]->get(id, "stave"),
                        map_decoder["ECALBarrel"]->get(id, "dlayer"),
                        map_decoder["ECALBarrel"]->get(id, "part"),
                        map_decoder["ECALBarrel"]->get(id, "slayer"),
                        map_decoder["ECALBarrel"]->get(id, "bar"));
      m_bar->setPosition(hit.second[0].getPosition());
      m_bar->setQ( hit.second[0].getEnergy(), hit.second[1].getEnergy() );
      m_bar->setT( hit.second[0].getT1(), hit.second[1].getT1() );
      m_barCol.push_back(m_bar);
    }

    m_DataCol.BarCol = m_barCol; 
    //Group bars to Blocks

    return StatusCode::SUCCESS;
  };
};

#endif
