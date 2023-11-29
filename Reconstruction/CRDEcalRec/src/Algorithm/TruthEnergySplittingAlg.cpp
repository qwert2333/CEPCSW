#ifndef _TRUTHENERGYSPLITTING_ALG_C
#define _TRUTHENERGYSPLITTING_ALG_C

#include "Algorithm/TruthEnergySplittingAlg.h"

StatusCode TruthEnergySplittingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_stringPars.find("ReadinAxisName")==settings.map_stringPars.end())  settings.map_stringPars["ReadinAxisName"] = "MergedAxis";
  if(settings.map_stringPars.find("OutputClusName")==settings.map_stringPars.end())  settings.map_stringPars["OutputClusName"] = "ESHalfCluster";
  if(settings.map_stringPars.find("OutputTowerName")==settings.map_stringPars.end()) settings.map_stringPars["OutputTowerName"] = "TruthESTower";

  if(settings.map_floatPars.find("Eth_HFClus")==settings.map_floatPars.end())        settings.map_floatPars["Eth_HFClus"] = 0.05;
  if(settings.map_intPars.find("th_Nhit")==settings.map_intPars.end())               settings.map_intPars["th_Nhit"] = 2;
  if(settings.map_boolPars.find("CompactHFCluster")==settings.map_boolPars.end())    settings.map_boolPars["CompactHFCluster"] = true;

  return StatusCode::SUCCESS;
};

StatusCode TruthEnergySplittingAlg::Initialize( PandoraPlusDataCol& m_datacol ){

  p_HalfClusterU.clear();
  p_HalfClusterV.clear();
  m_newClusUCol.clear();
  m_newClusVCol.clear();


  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColU"].size(); ih++)
    p_HalfClusterU.push_back( m_datacol.map_HalfCluster["HalfClusterColU"][ih].get() );
  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColV"].size(); ih++)
    p_HalfClusterV.push_back( m_datacol.map_HalfCluster["HalfClusterColV"][ih].get() );

  return StatusCode::SUCCESS;
};

StatusCode TruthEnergySplittingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  for(int ih=0; ih<p_HalfClusterU.size(); ih++){
    std::vector<const PandoraPlus::CaloHalfCluster*> m_axisUCol;
    if( settings.map_stringPars["ReadinAxisName"] == "AllAxis" ) m_axisUCol = p_HalfClusterU[ih]->getAllHalfClusterCol();
    else m_axisUCol = p_HalfClusterU[ih]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"]);    

    std::vector<const PandoraPlus::Calo1DCluster*> m_1dclusCol = p_HalfClusterU[ih]->getCluster();
    
    for(int iax=0; iax<m_axisUCol.size(); iax++){
      edm4hep::MCParticle truthMCP_axis = m_axisUCol[iax]->getLeadingMCP();
      std::shared_ptr<PandoraPlus::CaloHalfCluster> m_newClus = std::make_shared<PandoraPlus::CaloHalfCluster>();

      for(int ish=0; ish<m_1dclusCol.size(); ish++){
        std::shared_ptr<PandoraPlus::Calo1DCluster> m_shower = std::make_shared<PandoraPlus::Calo1DCluster>(); 
        for(int ibar=0; ibar<m_1dclusCol[ish]->getBars().size(); ibar++){
          auto truthMap = m_1dclusCol[ish]->getBars()[ibar]->getLinkedMCP();
          for(auto& iter: truthMap){
            if(!(iter.first==truthMCP_axis)) continue;
            if(iter.second<0.05) continue;

            //This bar has contribution from MCP. Split it. 
            auto m_bar = m_1dclusCol[ish]->getBars()[ibar]->Clone();
            m_bar->setQ(m_bar->getQ1()*iter.second, m_bar->getQ2()*iter.second);

            m_shower->addUnit(m_bar.get());
            m_datacol.map_BarCol["bkBars"].push_back(m_bar);
          }
        }
        m_shower->getLinkedMCPfromUnit();
        m_shower->setSeed();
        m_shower->setIDInfo();
        m_newClus->addUnit(m_shower.get());
        m_datacol.map_1DCluster["bk1DCluster"].push_back(m_shower);
      }//End loop 1DClusters in HFCluster

      m_newClusUCol.push_back(m_newClus);
      m_datacol.map_HalfCluster["bkHalfCluster"].push_back(m_newClus);
    }//End loop axis
  }//End loop HalfClusters.

  //Merge new HFClusters linked to the same MCP
  for(int iax=0; iax<m_newClusUCol.size(); iax++){
    for(int jax=iax+1; jax<m_newClusUCol.size(); jax++){
      if(m_newClusUCol[iax]->getLeadingMCP() == m_newClusUCol[jax]->getLeadingMCP()){
        m_newClusUCol[iax]->mergeHalfCluster( m_newClusUCol[jax].get() );
        m_newClusUCol.erase(m_newClusUCol.begin()+jax);
        jax--;
        if(iax>jax+1) iax--;
      }
    }
  }

  for(int ih=0; ih<p_HalfClusterV.size(); ih++){
    std::vector<const PandoraPlus::CaloHalfCluster*> m_axisVCol;
    if( settings.map_stringPars["ReadinAxisName"] == "AllAxis" ) m_axisVCol = p_HalfClusterV[ih]->getAllHalfClusterCol();
    else m_axisVCol = p_HalfClusterV[ih]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"]);

    std::vector<const PandoraPlus::Calo1DCluster*> m_1dclusCol = p_HalfClusterV[ih]->getCluster();

    for(int iax=0; iax<m_axisVCol.size(); iax++){
      edm4hep::MCParticle truthMCP_axis = m_axisVCol[iax]->getLeadingMCP();
      std::shared_ptr<PandoraPlus::CaloHalfCluster> m_newClus = std::make_shared<PandoraPlus::CaloHalfCluster>();

      for(int ish=0; ish<m_1dclusCol.size(); ish++){
        std::shared_ptr<PandoraPlus::Calo1DCluster> m_shower = std::make_shared<PandoraPlus::Calo1DCluster>();
        for(int ibar=0; ibar<m_1dclusCol[ish]->getBars().size(); ibar++){
          auto truthMap = m_1dclusCol[ish]->getBars()[ibar]->getLinkedMCP();
          for(auto& iter: truthMap){
            if(!(iter.first==truthMCP_axis)) continue;
            if(iter.second<0.05) continue;

            //This bar has contribution from MCP. Split it.
            auto m_bar = m_1dclusCol[ish]->getBars()[ibar]->Clone();
            m_bar->setQ(m_bar->getQ1()*iter.second, m_bar->getQ2()*iter.second);

            m_shower->addUnit(m_bar.get());
            m_datacol.map_BarCol["bkBars"].push_back(m_bar);
          }
        }
        m_shower->getLinkedMCPfromUnit();
        m_shower->setSeed();
        m_shower->setIDInfo();
        m_newClus->addUnit(m_shower.get());
        m_datacol.map_1DCluster["bk1DCluster"].push_back(m_shower);
      }//End loop 1DClusters in HFCluster

      m_newClusVCol.push_back(m_newClus);
      m_datacol.map_HalfCluster["bkHalfCluster"].push_back(m_newClus);
    }//End loop axis
  }//End loop HalfClusters.

  //Merge new HFClusters linked to the same MCP
  for(int iax=0; iax<m_newClusVCol.size(); iax++){
    for(int jax=iax+1; jax<m_newClusVCol.size(); jax++){
      if(m_newClusVCol[iax]->getLeadingMCP() == m_newClusVCol[jax]->getLeadingMCP()){
        m_newClusVCol[iax]->mergeHalfCluster( m_newClusVCol[jax].get() );
        m_newClusVCol.erase(m_newClusVCol.begin()+jax);
        jax--;
        if(iax>jax+1) iax--;
      }
    }
  }

  m_datacol.map_HalfCluster[settings.map_stringPars["OutputClusName"]+"U"] = m_newClusUCol;
  m_datacol.map_HalfCluster[settings.map_stringPars["OutputClusName"]+"V"] = m_newClusVCol;  

  //Make tower
  m_towerCol.clear();
  HalfClusterToTowers(m_newClusUCol, m_newClusVCol, m_towerCol);
  m_datacol.map_CaloCluster[settings.map_stringPars["OutputTowerName"]] = m_towerCol;

  return StatusCode::SUCCESS;
};

StatusCode TruthEnergySplittingAlg::ClearAlgorithm(){
  p_HalfClusterU.clear();
  p_HalfClusterV.clear();
  m_newClusUCol.clear();
  m_newClusVCol.clear();

  return StatusCode::SUCCESS;
};


StatusCode TruthEnergySplittingAlg::HalfClusterToTowers( std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& m_halfClusU,
                                                         std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& m_halfClusV,
                                                         std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>& m_towers ){


  m_towers.clear();

  std::map<std::vector<int>, std::vector<const PandoraPlus::Calo2DCluster*> > map_2DCluster;
  std::map<std::vector<int>, std::vector<PandoraPlus::CaloHalfCluster*> > map_HalfClusterU;
  std::map<std::vector<int>, std::vector<PandoraPlus::CaloHalfCluster*> > map_HalfClusterV;

  //Split CaloHalfClusterU
  for(int il=0; il<m_halfClusU.size(); il++){
    if(m_halfClusU[il]->getCluster().size()==0) {std::cout<<"WARNING: Have an empty CaloHalfCluster! Skip it! "<<std::endl; continue;}

    //HalfCluster does not cover tower:
    if( m_halfClusU[il]->getTowerID().size()==1 && m_halfClusU[il]->getCluster().size()>=settings.map_intPars["th_Nhit"]){
      std::vector<int> cl_towerID =  m_halfClusU[il]->getTowerID()[0];
      if(settings.map_boolPars["CompactHFCluster"]) m_halfClusU[il]->mergeClusterInLayer();
      map_HalfClusterU[cl_towerID].push_back(m_halfClusU[il].get());
      continue;
    }

    //CaloHalfCluster covers towers: Loop check showers.
    std::map<std::vector<int>, PandoraPlus::CaloHalfCluster* > tmp_LongiClusMaps; tmp_LongiClusMaps.clear();
    for(int is=0; is<m_halfClusU[il]->getCluster().size(); is++){
      const PandoraPlus::Calo1DCluster* p_shower = m_halfClusU[il]->getCluster()[is];

      if(p_shower->getSeeds().size()==0){
        std::cout<<"  HalfClusterToTowers ERROR: No Seed in 1DShower, Check! "<<std::endl;
        continue;
      }
      std::vector<int> seedID(3);
      seedID[0] = p_shower->getSeeds()[0]->getModule();
      seedID[1] = p_shower->getSeeds()[0]->getPart();
      seedID[2] = p_shower->getSeeds()[0]->getStave();

      if( tmp_LongiClusMaps.find( seedID )!=tmp_LongiClusMaps.end() ){
        tmp_LongiClusMaps[seedID]->addUnit( p_shower );
        tmp_LongiClusMaps[seedID]->setTowerID( seedID );
      }
      else{
        std::shared_ptr<PandoraPlus::CaloHalfCluster> tmp_clus = std::make_shared<PandoraPlus::CaloHalfCluster>();
        tmp_clus->addUnit( p_shower );
        tmp_clus->setTowerID( seedID );
        tmp_LongiClusMaps[seedID] = tmp_clus.get();
        m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(tmp_clus);
      }
      p_shower = nullptr;

    }

    //Connect cousins
    if(tmp_LongiClusMaps.size()>1){
      for(auto &iter : tmp_LongiClusMaps){
        for(auto &iter1 : tmp_LongiClusMaps){
          if(iter!= iter1 &&
             iter.second->getEnergy()>settings.map_floatPars["Eth_HFClus"] &&
             iter1.second->getEnergy()>settings.map_floatPars["Eth_HFClus"] &&
             iter.second->getCluster().size()>=settings.map_intPars["th_Nhit"] &&
             iter1.second->getCluster().size()>=settings.map_intPars["th_Nhit"] ){ iter.second->addCousinCluster(iter1.second); }
        }
      }
    }
    for(auto &iter : tmp_LongiClusMaps){
      if(iter.second->getEnergy()<settings.map_floatPars["Eth_HFClus"]) continue;
      iter.second->addHalfCluster("ParentCluster", m_halfClusU[il].get());
      for(int itrk=0; itrk<m_halfClusU[il]->getAssociatedTracks().size(); itrk++)
        iter.second->addAssociatedTrack( m_halfClusU[il]->getAssociatedTracks()[itrk] );
      if(iter.second->getCluster().size()>=settings.map_intPars["th_Nhit"]){
        if(settings.map_boolPars["CompactHFCluster"]){
          iter.second->mergeClusterInLayer();
          iter.second->setTowerID(iter.first);
        }
        map_HalfClusterU[iter.first].push_back(iter.second);
      }
    }

  }

  //Split CaloHalfClusterV
  for(int il=0; il<m_halfClusV.size(); il++){
    if(m_halfClusV[il]->getCluster().size()==0) {std::cout<<"WARNING: Have an empty CaloHalfCluster! Skip it! "<<std::endl; continue;}

    //HalfCluster does not cover tower:
    if( m_halfClusV[il]->getTowerID().size()==1 && m_halfClusV[il]->getCluster().size()>=settings.map_intPars["th_Nhit"]){
      std::vector<int> cl_towerID =  m_halfClusV[il]->getTowerID()[0];
      if(settings.map_boolPars["CompactHFCluster"]) m_halfClusV[il]->mergeClusterInLayer();
      map_HalfClusterV[cl_towerID].push_back(m_halfClusV[il].get());
      continue;
    }

    //CaloHalfCluster covers towers: Loop check showers.
    std::map<std::vector<int>, PandoraPlus::CaloHalfCluster* > tmp_LongiClusMaps; tmp_LongiClusMaps.clear();
    for(int is=0; is<m_halfClusV[il]->getCluster().size(); is++){
      const PandoraPlus::Calo1DCluster* p_shower = m_halfClusV[il]->getCluster()[is];
      if(p_shower->getSeeds().size()==0){
        std::cout<<"  HalfClusterToTowers ERROR: No Seed in 1DShower, Check! "<<std::endl;
        continue;
      }
      std::vector<int> seedID(3);
      seedID[0] = p_shower->getSeeds()[0]->getModule();
      seedID[1] = p_shower->getSeeds()[0]->getPart();
      seedID[2] = p_shower->getSeeds()[0]->getStave();

      if( tmp_LongiClusMaps.find( seedID )!=tmp_LongiClusMaps.end() ){
        tmp_LongiClusMaps[seedID]->addUnit( p_shower );
        tmp_LongiClusMaps[seedID]->setTowerID( seedID );
      }
      else{
        std::shared_ptr<PandoraPlus::CaloHalfCluster> tmp_clus = std::make_shared<PandoraPlus::CaloHalfCluster>();
        tmp_clus->addUnit( p_shower );
        tmp_clus->setTowerID( seedID );
        tmp_LongiClusMaps[seedID] = tmp_clus.get();
        m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(tmp_clus);
      }
      p_shower = nullptr;

    }

    //Connect cousins
    if(tmp_LongiClusMaps.size()>1){
      for(auto &iter : tmp_LongiClusMaps){
        for(auto &iter1 : tmp_LongiClusMaps){
          if(iter!= iter1 &&
             iter.second->getEnergy()>settings.map_floatPars["Eth_HFClus"] &&
             iter1.second->getEnergy()>settings.map_floatPars["Eth_HFClus"] &&
             iter.second->getCluster().size()>=settings.map_intPars["th_Nhit"] &&
             iter1.second->getCluster().size()>=settings.map_intPars["th_Nhit"] ){ iter.second->addCousinCluster(iter1.second); }
        }
      }
    }
    for(auto &iter : tmp_LongiClusMaps){
      if(iter.second->getEnergy()<settings.map_floatPars["Eth_HFClus"]) continue;
      iter.second->addHalfCluster("ParentCluster", m_halfClusV[il].get());
      for(int itrk=0; itrk<m_halfClusV[il]->getAssociatedTracks().size(); itrk++)
        iter.second->addAssociatedTrack( m_halfClusV[il]->getAssociatedTracks()[itrk] );
      if(iter.second->getCluster().size()>=settings.map_intPars["th_Nhit"]){
        if(settings.map_boolPars["CompactHFCluster"]){
          iter.second->mergeClusterInLayer();
          iter.second->setTowerID(iter.first);
        }
        map_HalfClusterV[iter.first].push_back(iter.second);
      }
    }

  }

  //Build 2DCluster
  for(auto &iterU : map_HalfClusterU){
    if( map_HalfClusterV.find(iterU.first)==map_HalfClusterV.end() ){
      iterU.second.clear();
      continue;
    }

    std::vector<PandoraPlus::CaloHalfCluster*> p_halfClusU = iterU.second;
    std::vector<PandoraPlus::CaloHalfCluster*> p_halfClusV = map_HalfClusterV[iterU.first];

    //Get ordered showers for looping in layers.
    std::map<int, std::vector<const PandoraPlus::Calo1DCluster*>> m_orderedShowerU; m_orderedShowerU.clear();
    std::map<int, std::vector<const PandoraPlus::Calo1DCluster*>> m_orderedShowerV; m_orderedShowerV.clear();

    for(int ic=0; ic<p_halfClusU.size(); ic++){
      for(int is=0; is<p_halfClusU.at(ic)->getCluster().size(); is++)
        m_orderedShowerU[p_halfClusU.at(ic)->getCluster()[is]->getDlayer()].push_back( p_halfClusU.at(ic)->getCluster()[is] );
    }
    for(int ic=0; ic<p_halfClusV.size(); ic++){
      for(int is=0; is<p_halfClusV.at(ic)->getCluster().size(); is++)
        m_orderedShowerV[p_halfClusV.at(ic)->getCluster()[is]->getDlayer()].push_back( p_halfClusV.at(ic)->getCluster()[is] );
    }
    p_halfClusU.clear(); p_halfClusV.clear();

    //Create super-layers (block)
    std::vector<const PandoraPlus::Calo2DCluster*> m_blocks; m_blocks.clear();
    for(auto &iter1 : m_orderedShowerU){
      if( m_orderedShowerV.find( iter1.first )==m_orderedShowerV.end() ) continue;
      std::shared_ptr<PandoraPlus::Calo2DCluster> tmp_block = std::make_shared<PandoraPlus::Calo2DCluster>();
      for(int is=0; is<iter1.second.size(); is++) tmp_block->addUnit( iter1.second.at(is) );
      for(int is=0; is<m_orderedShowerV[iter1.first].size(); is++) tmp_block->addUnit( m_orderedShowerV[iter1.first].at(is) );
      tmp_block->setTowerID( iterU.first );
      m_blocks.push_back( tmp_block.get() );
      m_bkCol.map_2DCluster["bk2DCluster"].push_back(tmp_block);
    }
    map_2DCluster[iterU.first] = m_blocks;
  }
  for(auto &iterV : map_HalfClusterV){
    if( map_HalfClusterU.find(iterV.first)==map_HalfClusterU.end() ){
      iterV.second.clear();
    }
  }

  //Form a tower:
  for(auto &iter : map_2DCluster){
    std::vector<int> m_towerID = iter.first;
    //Check cousin clusters:
    std::vector<PandoraPlus::CaloHalfCluster*> m_HFClusUInTower = map_HalfClusterU[m_towerID];
    for(auto &m_HFclus : m_HFClusUInTower){
      std::vector<const CaloHalfCluster*> tmp_delClus; tmp_delClus.clear();
      for(int ics=0; ics<m_HFclus->getHalfClusterCol("CousinCluster").size(); ics++){
        std::vector<int> tmp_towerID = m_HFclus->getHalfClusterCol("CousinCluster")[ics]->getTowerID()[0];

        if( map_2DCluster.find( tmp_towerID )==map_2DCluster.end() )
          tmp_delClus.push_back( m_HFclus->getHalfClusterCol("CousinCluster")[ics] );
      }
      for(int ics=0; ics<tmp_delClus.size(); ics++) m_HFclus->deleteCousinCluster( tmp_delClus[ics] );
    }

    std::vector<PandoraPlus::CaloHalfCluster*> m_HFClusVInTower = map_HalfClusterV[m_towerID];
    for(auto &m_HFclus : m_HFClusVInTower){
      std::vector<const CaloHalfCluster*> tmp_delClus; tmp_delClus.clear();
      for(int ics=0; ics<m_HFclus->getHalfClusterCol("CousinCluster").size(); ics++){
        std::vector<int> tmp_towerID = m_HFclus->getHalfClusterCol("CousinCluster")[ics]->getTowerID()[0];

        if( map_2DCluster.find( tmp_towerID )==map_2DCluster.end() )
          tmp_delClus.push_back( m_HFclus->getHalfClusterCol("CousinCluster")[ics] );
      }
      for(int ics=0; ics<tmp_delClus.size(); ics++) m_HFclus->deleteCousinCluster( tmp_delClus[ics] );
    }

    //Convert to const
    std::vector<const PandoraPlus::CaloHalfCluster*> const_HFClusU; const_HFClusU.clear();
    std::vector<const PandoraPlus::CaloHalfCluster*> const_HFClusV; const_HFClusV.clear();
    for(int ics=0; ics<m_HFClusUInTower.size(); ics++){ m_HFClusUInTower[ics]->getLinkedMCPfromUnit(); const_HFClusU.push_back(m_HFClusUInTower[ics]); }
    for(int ics=0; ics<m_HFClusVInTower.size(); ics++){ m_HFClusVInTower[ics]->getLinkedMCPfromUnit(); const_HFClusV.push_back(m_HFClusVInTower[ics]); }

    std::shared_ptr<PandoraPlus::Calo3DCluster> m_tower = std::make_shared<PandoraPlus::Calo3DCluster>();
    m_tower->addTowerID( m_towerID );
    for(int i2d=0; i2d<map_2DCluster[m_towerID].size(); i2d++) m_tower->addUnit(map_2DCluster[m_towerID][i2d]);
    m_tower->setHalfClusters( settings.map_stringPars["OutputClusName"]+"U", const_HFClusU,
                              settings.map_stringPars["OutputClusName"]+"V", const_HFClusV );
    m_towers.push_back(m_tower);
  }

  return StatusCode::SUCCESS;

}

#endif
