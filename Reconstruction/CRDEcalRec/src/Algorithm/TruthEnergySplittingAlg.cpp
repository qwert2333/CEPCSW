#ifndef _TRUTHENERGYSPLITTING_ALG_C
#define _TRUTHENERGYSPLITTING_ALG_C

#include "Algorithm/TruthEnergySplittingAlg.h"

StatusCode TruthEnergySplittingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_stringPars.find("ReadinAxisName")==settings.map_stringPars.end())  settings.map_stringPars["ReadinAxisName"] = "TruthAxis";
  if(settings.map_stringPars.find("OutputClusName")==settings.map_stringPars.end())  settings.map_stringPars["OutputClusName"] = "TruthESCluster";
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
  m_bkCol.Clear();

  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColU"].size(); ih++)
    p_HalfClusterU.push_back( m_datacol.map_HalfCluster["HalfClusterColU"][ih].get() );
  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColV"].size(); ih++)
    p_HalfClusterV.push_back( m_datacol.map_HalfCluster["HalfClusterColV"][ih].get() );

  return StatusCode::SUCCESS;
};

StatusCode TruthEnergySplittingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){
cout<<"  TruthEnergySplittingAlg: readin HFCluster size "<<p_HalfClusterU.size()<<", "<<p_HalfClusterV.size()<<endl;

  for(int ih=0; ih<p_HalfClusterU.size(); ih++){
    std::vector<const PandoraPlus::CaloHalfCluster*> m_axisUCol;
    if( settings.map_stringPars["ReadinAxisName"] == "AllAxis" ) m_axisUCol = p_HalfClusterU[ih]->getAllHalfClusterCol();
    else m_axisUCol = p_HalfClusterU[ih]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"]);    

    std::vector<const PandoraPlus::Calo1DCluster*> m_1dclusCol = p_HalfClusterU[ih]->getCluster();
cout<<"    In HFU #"<<ih<<": axis size "<<m_axisUCol.size()<<", 1DCluster size "<<m_1dclusCol.size()<<endl;
    
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

      if(m_newClus && m_newClus->getCluster().size()!=0){ 
        m_newClusUCol.push_back(m_newClus);
        m_datacol.map_HalfCluster["bkHalfCluster"].push_back(m_newClus);
      }
    }//End loop axis
  }//End loop HalfClusters.
cout<<"  Splitted HFClusterU size "<<m_newClusUCol.size()<<endl;

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
cout<<"  Splitted HFClusterU size after merging "<<m_newClusUCol.size()<<endl;

  for(int ih=0; ih<p_HalfClusterV.size(); ih++){
    std::vector<const PandoraPlus::CaloHalfCluster*> m_axisVCol;
    if( settings.map_stringPars["ReadinAxisName"] == "AllAxis" ) m_axisVCol = p_HalfClusterV[ih]->getAllHalfClusterCol();
    else m_axisVCol = p_HalfClusterV[ih]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"]);

    std::vector<const PandoraPlus::Calo1DCluster*> m_1dclusCol = p_HalfClusterV[ih]->getCluster();
cout<<"    In HFV #"<<ih<<": axis size "<<m_axisVCol.size()<<", 1DCluster size "<<m_1dclusCol.size()<<endl;

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

      if(m_newClus && m_newClus->getCluster().size()!=0){
        m_newClusVCol.push_back(m_newClus);
        m_datacol.map_HalfCluster["bkHalfCluster"].push_back(m_newClus);
      }
    }//End loop axis
  }//End loop HalfClusters.
cout<<"  Splitted HFClusterV size "<<m_newClusVCol.size()<<endl;

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
cout<<"  Splitted HFClusterU size after merging "<<m_newClusVCol.size()<<endl;

  m_datacol.map_HalfCluster[settings.map_stringPars["OutputClusName"]+"U"] = m_newClusUCol;
  m_datacol.map_HalfCluster[settings.map_stringPars["OutputClusName"]+"V"] = m_newClusVCol;  

  //Make tower
  m_towerCol.clear();
  HalfClusterToTowers(m_newClusUCol, m_newClusVCol, m_towerCol);
  m_datacol.map_CaloCluster[settings.map_stringPars["OutputTowerName"]] = m_towerCol;

/*
cout<<"  After splitting: tower size "<<m_towerCol.size()<<". Print Tower: "<<endl;
for(auto it : m_towerCol){
  std::vector<const CaloHalfCluster*> m_HFClusUInTower = it->getHalfClusterUCol(settings.map_stringPars["OutputClusName"]+"U");
  std::vector<const CaloHalfCluster*> m_HFClusVInTower = it->getHalfClusterVCol(settings.map_stringPars["OutputClusName"]+"V");

cout<<"Check tower ID: ";
for(int i=0; i<it->getTowerID().size(); i++) printf("[%d, %d, %d], ", it->getTowerID()[i][0], it->getTowerID()[i][1], it->getTowerID()[i][2]);
cout<<endl;
  printf("    In Tower [%d, %d, %d], ", it->getTowerID()[0][0], it->getTowerID()[0][1], it->getTowerID()[0][2] );
  printf("    HalfCluster size: (%d, %d) \n", m_HFClusUInTower.size(), m_HFClusVInTower.size() );
  cout<<"    Loop print HalfClusterU: "<<endl;
  for(int icl=0; icl<m_HFClusUInTower.size(); icl++){
    cout<<"      In HFClusU #"<<icl<<": shower size = "<<m_HFClusUInTower[icl]->getCluster().size()<<", En = "<<m_HFClusUInTower[icl]->getEnergy();
    printf(", Position (%.3f, %.3f, %.3f), address %p ",m_HFClusUInTower[icl]->getPos().x(), m_HFClusUInTower[icl]->getPos().y(), m_HFClusUInTower[icl]->getPos().z(), m_HFClusUInTower[icl]);
    printf(", cousin size %d, address: ", m_HFClusUInTower[icl]->getHalfClusterCol("CousinCluster").size());
    for(int ics=0; ics<m_HFClusUInTower[icl]->getHalfClusterCol("CousinCluster").size(); ics++) printf("%p, ", m_HFClusUInTower[icl]->getHalfClusterCol("CousinCluster")[ics]);
    printf(", track size %d, address: ", m_HFClusUInTower[icl]->getAssociatedTracks().size());
    for(int itrk=0; itrk<m_HFClusUInTower[icl]->getAssociatedTracks().size(); itrk++) printf("%p, ", m_HFClusUInTower[icl]->getAssociatedTracks()[itrk]);
    cout<<endl;
    for(auto ish : m_HFClusUInTower[icl]->getCluster()){
      printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, Address %p \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(),  ish->getNseeds(), ish );
    }
  }
  cout<<endl;

  cout<<"    Loop print HalfClusterV: "<<endl;
  for(int icl=0; icl<m_HFClusVInTower.size(); icl++){
    cout<<"      In HFClusV #"<<icl<<": shower size = "<<m_HFClusVInTower[icl]->getCluster().size()<<", En = "<<m_HFClusVInTower[icl]->getEnergy();
    printf(", Position (%.3f, %.3f, %.3f), address %p ",m_HFClusVInTower[icl]->getPos().x(), m_HFClusVInTower[icl]->getPos().y(), m_HFClusVInTower[icl]->getPos().z(), m_HFClusVInTower[icl]);
    printf(", cousin size %d, address: ", m_HFClusVInTower[icl]->getHalfClusterCol("CousinCluster").size());
    for(int ics=0; ics<m_HFClusVInTower[icl]->getHalfClusterCol("CousinCluster").size(); ics++) printf("%p, ", m_HFClusVInTower[icl]->getHalfClusterCol("CousinCluster")[ics]);
    printf(", track size %d, address: ", m_HFClusVInTower[icl]->getAssociatedTracks().size());
    for(int itrk=0; itrk<m_HFClusVInTower[icl]->getAssociatedTracks().size(); itrk++) printf("%p, ", m_HFClusVInTower[icl]->getAssociatedTracks()[itrk]);
    cout<<endl;
    for(auto ish : m_HFClusVInTower[icl]->getCluster()){
      printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, Address %p \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(), ish->getNseeds(), ish );
    }
  }
  cout<<endl;
}
*/

  m_datacol.map_1DCluster["bk1DCluster"].insert( m_datacol.map_1DCluster["bk1DCluster"].end(), m_bkCol.map_1DCluster["bk1DCluster"].begin(), m_bkCol.map_1DCluster["bk1DCluster"].end() );
  m_datacol.map_2DCluster["bk2DCluster"].insert( m_datacol.map_2DCluster["bk2DCluster"].end(), m_bkCol.map_2DCluster["bk2DCluster"].begin(), m_bkCol.map_2DCluster["bk2DCluster"].end() );
  m_datacol.map_HalfCluster["bkHalfCluster"].insert( m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_bkCol.map_HalfCluster["bkHalfCluster"].begin(), m_bkCol.map_HalfCluster["bkHalfCluster"].end() );


  return StatusCode::SUCCESS;
};

StatusCode TruthEnergySplittingAlg::ClearAlgorithm(){
  p_HalfClusterU.clear();
  p_HalfClusterV.clear();
  m_newClusUCol.clear();
  m_newClusVCol.clear();

  m_bkCol.Clear();
  return StatusCode::SUCCESS;
};


StatusCode TruthEnergySplittingAlg::HalfClusterToTowers( std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& m_halfClusU,
                                                         std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& m_halfClusV,
                                                         std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>& m_towers ){
  m_towers.clear();

  std::map<std::vector<int>, std::vector<const PandoraPlus::CaloUnit*>> map_barCol; map_barCol.clear();
  for(int il=0; il<m_halfClusU.size(); il++){
    for(int is=0; is<m_halfClusU[il]->getCluster().size(); is++){
      const PandoraPlus::Calo1DCluster* p_shower = m_halfClusU[il]->getCluster()[is];
      for(int ibar=0; ibar<p_shower->getBars().size(); ibar++){
        std::vector<int> barID(3);
        barID[0] = p_shower->getBars()[ibar]->getModule();
        barID[1] = p_shower->getBars()[ibar]->getPart();
        barID[2] = p_shower->getBars()[ibar]->getStave();
        map_barCol[barID].push_back(p_shower->getBars()[ibar]);
      }
    }
  }
  for(int il=0; il<m_halfClusV.size(); il++){
    for(int is=0; is<m_halfClusV[il]->getCluster().size(); is++){
      const PandoraPlus::Calo1DCluster* p_shower = m_halfClusV[il]->getCluster()[is];
      for(int ibar=0; ibar<p_shower->getBars().size(); ibar++){
        std::vector<int> barID(3);
        barID[0] = p_shower->getBars()[ibar]->getModule();
        barID[1] = p_shower->getBars()[ibar]->getPart();
        barID[2] = p_shower->getBars()[ibar]->getStave();
        map_barCol[barID].push_back(p_shower->getBars()[ibar]);
      }
    }
  }

  //Re-build the objects in tower
  for(auto& itower: map_barCol){
    std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_newHFClusUCol;
    std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_newHFClusVCol;
    std::vector<std::shared_ptr<PandoraPlus::Calo2DCluster>> m_new2DClusCol;

    //Get map for MCP-bars
    std::map<edm4hep::MCParticle, std::vector<const PandoraPlus::CaloUnit*> > map_matchBar;
    for(int ibar=0; ibar<itower.second.size(); ibar++){
      edm4hep::MCParticle mcp = itower.second[ibar]->getLeadingMCP();
      map_matchBar[mcp].push_back(itower.second[ibar]);
    }

    //Build objects for MCP
    for(auto& iter: map_matchBar){
//cout<<"    Print bar info"<<endl;
//for(int ibar=0; ibar<iter.second.size(); ibar++){
//printf("    cellID (%d, %d, %d, %d, %d, %d), En %.3f \n", iter.second[ibar]->getModule(), iter.second[ibar]->getPart(), iter.second[ibar]->getStave(), iter.second[ibar]->getDlayer(), iter.second[ibar]->getSlayer(), iter.second[ibar]->getBar(), iter.second[ibar]->getEnergy());
//}
      std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>> m_new1DClusUCol;
      std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>> m_new1DClusVCol;

      //Ordered bars
      std::map<int, std::vector<const PandoraPlus::CaloUnit*> > m_orderedBars; m_orderedBars.clear();
      for(int ibar=0; ibar<iter.second.size(); ibar++)
        m_orderedBars[iter.second[ibar]->getDlayer()].push_back(iter.second[ibar]);
      //1D&2D cluster
      for(auto& iter_bar: m_orderedBars){
//cout<<"    In layer #"<<iter_bar.first<<": bar size "<<iter_bar.second.size()<<endl;

        std::shared_ptr<PandoraPlus::Calo1DCluster> tmp_new1dclusterU = std::make_shared<PandoraPlus::Calo1DCluster>();
        std::shared_ptr<PandoraPlus::Calo1DCluster> tmp_new1dclusterV = std::make_shared<PandoraPlus::Calo1DCluster>();
        for(int ibar=0; ibar<iter_bar.second.size(); ibar++){
          if(iter_bar.second[ibar]->getSlayer()==0) tmp_new1dclusterU->addUnit( iter_bar.second[ibar] );
          if(iter_bar.second[ibar]->getSlayer()==1) tmp_new1dclusterV->addUnit( iter_bar.second[ibar] );
        }
        if(tmp_new1dclusterU && tmp_new1dclusterU->getBars().size()>0){
          tmp_new1dclusterU->getLinkedMCPfromUnit();
          tmp_new1dclusterU->setSeed();
          tmp_new1dclusterU->setIDInfo();
          m_new1DClusUCol.push_back(tmp_new1dclusterU);
          m_bkCol.map_1DCluster["bk1DCluster"].push_back(tmp_new1dclusterU);
        }
        if(tmp_new1dclusterV && tmp_new1dclusterV->getBars().size()>0){
          tmp_new1dclusterV->getLinkedMCPfromUnit();
          tmp_new1dclusterV->setSeed();
          tmp_new1dclusterV->setIDInfo();
          m_new1DClusVCol.push_back(tmp_new1dclusterV);
          m_bkCol.map_1DCluster["bk1DCluster"].push_back(tmp_new1dclusterV);
        }
        if(tmp_new1dclusterU && tmp_new1dclusterV && tmp_new1dclusterU->getBars().size()>0 && tmp_new1dclusterV->getBars().size()>0){
          std::shared_ptr<PandoraPlus::Calo2DCluster> tmp_block = std::make_shared<PandoraPlus::Calo2DCluster>();
          for(int ibar=0; ibar<iter_bar.second.size(); ibar++) tmp_block->addBar(iter_bar.second[ibar]);
          tmp_block->addUnit(tmp_new1dclusterU.get());
          tmp_block->addUnit(tmp_new1dclusterV.get());
          tmp_block->setTowerID( itower.first );
          m_new2DClusCol.push_back(tmp_block);
          m_bkCol.map_2DCluster["bk2DCluster"].push_back(tmp_block);
        }
      }
//cout<<"    1DClusU size "<<m_new1DClusUCol.size()<<", 1DClusV size "<<m_new1DClusVCol.size();
//cout<<", 2DClus size "<<m_new2DClusCol.size()<<endl;
      //Half cluster
      std::shared_ptr<PandoraPlus::CaloHalfCluster> m_newHFClusterU = std::make_shared<PandoraPlus::CaloHalfCluster>();
      std::shared_ptr<PandoraPlus::CaloHalfCluster> m_newHFClusterV = std::make_shared<PandoraPlus::CaloHalfCluster>();
      for(int i1d=0; i1d<m_new1DClusUCol.size(); i1d++)
        m_newHFClusterU->addUnit(m_new1DClusUCol[i1d].get());
      for(int i1d=0; i1d<m_new1DClusVCol.size(); i1d++)
        m_newHFClusterV->addUnit(m_new1DClusVCol[i1d].get());

      m_newHFClusterU->getLinkedMCPfromUnit();
      m_newHFClusterV->getLinkedMCPfromUnit();
      m_newHFClusUCol.push_back(m_newHFClusterU);
      m_newHFClusVCol.push_back(m_newHFClusterV);
      m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(m_newHFClusterU);
      m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(m_newHFClusterV);
    }

    //Form a tower
    std::vector<const PandoraPlus::CaloHalfCluster*> const_HFClusU; const_HFClusU.clear();
    std::vector<const PandoraPlus::CaloHalfCluster*> const_HFClusV; const_HFClusV.clear();
    for(int ics=0; ics<m_newHFClusUCol.size(); ics++){ const_HFClusU.push_back(m_newHFClusUCol[ics].get()); }
    for(int ics=0; ics<m_newHFClusVCol.size(); ics++){ const_HFClusV.push_back(m_newHFClusVCol[ics].get()); }

    std::shared_ptr<PandoraPlus::Calo3DCluster> m_tower = std::make_shared<PandoraPlus::Calo3DCluster>();
    m_tower->addTowerID(itower.first);
    for(int i2d=0; i2d<m_new2DClusCol.size(); i2d++) m_tower->addUnit(m_new2DClusCol[i2d].get());
    m_tower->setHalfClusters( settings.map_stringPars["OutputECALHalfClusters"]+"U", const_HFClusU,
                              settings.map_stringPars["OutputECALHalfClusters"]+"V", const_HFClusV );
    m_towers.push_back(m_tower);
  }

  return StatusCode::SUCCESS;
}


#endif
