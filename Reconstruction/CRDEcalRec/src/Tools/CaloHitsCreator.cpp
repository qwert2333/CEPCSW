#ifndef TRACK_CREATOR_C
#define TRACK_CREATOR_C

#include "Tools/CaloHitsCreator.h"

namespace PandoraPlus{
  CaloHitsCreator::CaloHitsCreator(const Settings& m_settings) : settings( m_settings ){

  }; 

  StatusCode CaloHitsCreator::CreateCaloHits( PandoraPlusDataCol& m_DataCol, 
                                              std::vector<DataHandle<edm4hep::CalorimeterHitCollection>*>& r_CaloHitCols, 
                                              std::map<std::string, dd4hep::DDSegmentation::BitFieldCoder*>& map_decoder,
                                              DataHandle<edm4hep::MCRecoCaloParticleAssociationCollection>* r_MCParticleRecoCaloCol )
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

      m_DataCol.collectionMap_CaloHit[settings.map_stringVecPars.at("CaloHitCollections")[icol]] = m_HitCol; 
    }

    //Convert to local objects: 
    for(auto iter : m_DataCol.collectionMap_CaloHit){
      if( settings.map_stringPars.at("EcalType")=="BarEcal" && iter.first == "ECALBarrel") continue; 
      
      std::vector<std::shared_ptr<PandoraPlus::CaloHit>> m_hitCol; m_hitCol.clear();

      for(int ihit=0; ihit<iter.second.size(); ihit++){
        //PandoraPlus::CaloHit* m_hit = new PandoraPlus::CaloHit();
        std::shared_ptr<PandoraPlus::CaloHit> m_hit = std::make_shared<PandoraPlus::CaloHit>();

        m_hit->setcellID( iter.second[ihit].getCellID() );
        m_hit->setLayer( map_decoder[iter.first]->get(iter.second[ihit].getCellID(), "layer") );

        TVector3 pos( iter.second[ihit].getPosition().x, iter.second[ihit].getPosition().y, iter.second[ihit].getPosition().z );
        m_hit->setPosition( pos );
        m_hit->setEnergy( iter.second[ihit].getEnergy() );

        m_hitCol.push_back( m_hit );
        //m_DataCol.bk_HitCol.push_back( m_hit );
      }
      m_DataCol.map_CaloHit[iter.first] = m_hitCol;
    }


    //Convert to local objects: CalorimeterHit to CaloUnit (For ECALBarrel only)
    if(settings.map_stringPars.at("EcalType")=="BarEcal"){
      std::vector<std::shared_ptr<PandoraPlus::CaloUnit>> m_barCol; m_barCol.clear(); 

      const edm4hep::MCRecoCaloParticleAssociationCollection* const_MCPCaloAssoCol = r_MCParticleRecoCaloCol->get();   
      auto CaloHits = m_DataCol.collectionMap_CaloHit["ECALBarrel"]; 
      std::map<std::uint64_t, std::vector<PandoraPlus::CaloUnit> > map_cellID_hits; map_cellID_hits.clear();
      for(auto& hit : CaloHits){ 
        PandoraPlus::CaloUnit m_bar; 
        m_bar.setcellID(hit.getCellID());
        m_bar.setPosition( TVector3(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z) );
        m_bar.setQ(hit.getEnergy(), hit.getEnergy());
        m_bar.setT(hit.getTime(), hit.getTime());
        for(int ilink=0; ilink<const_MCPCaloAssoCol->size(); ilink++){
          if( hit == const_MCPCaloAssoCol->at(ilink).getRec() ) m_bar.addLinkedMCP( std::make_pair(const_MCPCaloAssoCol->at(ilink).getSim(), const_MCPCaloAssoCol->at(ilink).getWeight()) );
        } 
        map_cellID_hits[hit.getCellID()].push_back(m_bar);
      }
   
      for(auto& hit : map_cellID_hits){
        if(hit.second.size()!=2){ std::cout<<"WARNING: didn't find correct hit pairs! "<<std::endl; continue; }
   
        //PandoraPlus::CaloUnit* m_bar = new PandoraPlus::CaloUnit(); 
        std::shared_ptr<PandoraPlus::CaloUnit> m_bar = std::make_shared<PandoraPlus::CaloUnit>();
   
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

        //---oooOOO000OOOooo---set bar length---oooOOO000OOOooo---
        //TODO: reading bar length from geosvc. 
        double t_bar_length;
        if(m_bar->getSlayer()==1) { t_bar_length = 600.; }
        else { t_bar_length = 480.-m_bar->getDlayer()*10.; }
        m_bar->setBarLength(t_bar_length);
        //---oooOOO000OOOooo---set bar length---oooOOO000OOOooo---

        //add MCParticle link
        for(int ilink=0; ilink<hit.second[0].getLinkedMCP().size(); ilink++) m_bar->addLinkedMCP( hit.second[0].getLinkedMCP()[ilink] );

        m_barCol.push_back(m_bar);  //Save for later use in algorithms
      }

   
      m_DataCol.map_BarCol["BarCol"] = m_barCol; 
      //Group bars to Blocks
    }
    return StatusCode::SUCCESS;
  };

};

#endif
