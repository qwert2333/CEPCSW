#ifndef _TRUTHENERGYSPLITTING_ALG_C
#define _TRUTHENERGYSPLITTING_ALG_C

#include "Algorithm/TruthEnergySplittingAlg.h"

StatusCode TruthEnergySplittingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_stringPars.find("ReadinAxisName")==settings.map_stringPars.end())  settings.map_stringPars["ReadinAxisName"] = "TruthMergedAxis";
  if(settings.map_stringPars.find("OutputClusName")==settings.map_stringPars.end())  settings.map_stringPars["OutputClusName"] = "TruthESCluster";
  if(settings.map_stringPars.find("OutputTowerName")==settings.map_stringPars.end()) settings.map_stringPars["OutputTowerName"] = "TruthESTower";

  if(settings.map_floatPars.find("Eth_HFClus")==settings.map_floatPars.end())        settings.map_floatPars["Eth_HFClus"] = 0.05;
  if(settings.map_intPars.find("th_Nhit")==settings.map_intPars.end())               settings.map_intPars["th_Nhit"] = 2;
  if(settings.map_boolPars.find("CompactHFCluster")==settings.map_boolPars.end())    settings.map_boolPars["CompactHFCluster"] = 1;

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
//cout<<"  TruthEnergySplittingAlg: readin HFCluster size "<<p_HalfClusterU.size()<<", "<<p_HalfClusterV.size()<<endl;

  for(int ih=0; ih<p_HalfClusterU.size(); ih++){
    std::vector<const PandoraPlus::CaloHalfCluster*> m_axisUCol;
    if( settings.map_stringPars["ReadinAxisName"] == "AllAxis" ) m_axisUCol = p_HalfClusterU[ih]->getAllHalfClusterCol();
    else m_axisUCol = p_HalfClusterU[ih]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"]);    

    std::vector<const PandoraPlus::Calo1DCluster*> m_1dclusCol = p_HalfClusterU[ih]->getCluster();
//cout<<"    In HFU #"<<ih<<": axis size "<<m_axisUCol.size()<<", 1DCluster size "<<m_1dclusCol.size()<<endl;
    
    for(int iax=0; iax<m_axisUCol.size(); iax++){
      edm4hep::MCParticle truthMCP_axis = m_axisUCol[iax]->getLeadingMCP();
      std::shared_ptr<PandoraPlus::CaloHalfCluster> m_newClus = std::make_shared<PandoraPlus::CaloHalfCluster>();
//cout<<"      Axis #"<<iax<<": track size "<<m_axisUCol[iax]->getAssociatedTracks().size()<<", truth MC pid "<<truthMCP_axis.getPDG()<<endl;

      for(int ish=0; ish<m_1dclusCol.size(); ish++){
        std::shared_ptr<PandoraPlus::Calo1DCluster> m_shower = std::make_shared<PandoraPlus::Calo1DCluster>(); 
//cout<<"        1D cluster #"<<ish<<": bar size "<<m_1dclusCol[ish]->getBars().size();
        for(int ibar=0; ibar<m_1dclusCol[ish]->getBars().size(); ibar++){
          auto truthMap = m_1dclusCol[ish]->getBars()[ibar]->getLinkedMCP();
          for(auto& iter: truthMap){
            if(!(iter.first==truthMCP_axis)) continue;
            if(iter.second<0.05) continue;

            //This bar has contribution from MCP. Split it. 
            auto m_bar = m_1dclusCol[ish]->getBars()[ibar]->Clone();
            m_bar->setQ(m_bar->getQ1()*iter.second, m_bar->getQ2()*iter.second);
            m_bar->setLinkedMCP( std::vector<std::pair<edm4hep::MCParticle, float>>{std::make_pair(iter.first, 1)} );

            m_shower->addUnit(m_bar.get());
            m_datacol.map_BarCol["bkBars"].push_back(m_bar);
          }
        }

        if(m_shower && m_shower->getBars().size()>0){
//cout<<". Make a new shower. "<<endl;
          m_shower->getLinkedMCPfromUnit();
          m_shower->setSeed();
          m_shower->setIDInfo();
          m_newClus->addUnit(m_shower.get());
          m_datacol.map_1DCluster["bk1DCluster"].push_back(m_shower);
        }
      }//End loop 1DClusters in HFCluster

      if(m_newClus && m_newClus->getCluster().size()!=0){ 
        for(int itrk=0; itrk<m_axisUCol[iax]->getAssociatedTracks().size(); itrk++) m_newClus->addAssociatedTrack(m_axisUCol[iax]->getAssociatedTracks()[itrk]);
        m_newClus->getLinkedMCPfromUnit();
        m_newClus->fitAxis("");
        m_newClus->mergeClusterInLayer();
        m_newClusUCol.push_back(m_newClus);
        m_datacol.map_HalfCluster["bkHalfCluster"].push_back(m_newClus);
      }
    }//End loop axis
  }//End loop HalfClusters.
//cout<<"  Splitted HFClusterU size "<<m_newClusUCol.size()<<endl;

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
//cout<<"  Splitted HFClusterU size after merging "<<m_newClusUCol.size()<<endl;

  for(int ih=0; ih<p_HalfClusterV.size(); ih++){
    std::vector<const PandoraPlus::CaloHalfCluster*> m_axisVCol;
    if( settings.map_stringPars["ReadinAxisName"] == "AllAxis" ) m_axisVCol = p_HalfClusterV[ih]->getAllHalfClusterCol();
    else m_axisVCol = p_HalfClusterV[ih]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"]);

    std::vector<const PandoraPlus::Calo1DCluster*> m_1dclusCol = p_HalfClusterV[ih]->getCluster();
//cout<<"    In HFV #"<<ih<<": axis size "<<m_axisVCol.size()<<", 1DCluster size "<<m_1dclusCol.size()<<endl;

    for(int iax=0; iax<m_axisVCol.size(); iax++){
      edm4hep::MCParticle truthMCP_axis = m_axisVCol[iax]->getLeadingMCP();
      std::shared_ptr<PandoraPlus::CaloHalfCluster> m_newClus = std::make_shared<PandoraPlus::CaloHalfCluster>();
//cout<<"      Axis #"<<iax<<": track size "<<m_axisVCol[iax]->getAssociatedTracks().size()<<", truth MC pid "<<truthMCP_axis.getPDG()<<endl;

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
            m_bar->setLinkedMCP( std::vector<std::pair<edm4hep::MCParticle, float>>{std::make_pair(iter.first, 1)} );

            m_shower->addUnit(m_bar.get());
            m_datacol.map_BarCol["bkBars"].push_back(m_bar);
          }
        }
        if(m_shower && m_shower->getBars().size()>0){
          m_shower->getLinkedMCPfromUnit();
          m_shower->setSeed();
          m_shower->setIDInfo();
          m_newClus->addUnit(m_shower.get());
          m_datacol.map_1DCluster["bk1DCluster"].push_back(m_shower);
        }
      }//End loop 1DClusters in HFCluster

      if(m_newClus && m_newClus->getCluster().size()!=0){
        for(int itrk=0; itrk<m_axisVCol[iax]->getAssociatedTracks().size(); itrk++) m_newClus->addAssociatedTrack(m_axisVCol[iax]->getAssociatedTracks()[itrk]);
        m_newClus->getLinkedMCPfromUnit();
        m_newClus->fitAxis("");
        m_newClus->mergeClusterInLayer();
        m_newClusVCol.push_back(m_newClus);
        m_datacol.map_HalfCluster["bkHalfCluster"].push_back(m_newClus);
      }
    }//End loop axis
  }//End loop HalfClusters.
//cout<<"  Splitted HFClusterV size "<<m_newClusVCol.size()<<endl;

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
//cout<<"  Splitted HFClusterU size after merging "<<m_newClusVCol.size()<<endl;

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
  //std::vector<const CaloHalfCluster*> m_HFClusUInTower, m_HFClusVInTower; 
  //for(int i=0; i<m_newClusUCol.size(); i++) m_HFClusUInTower.push_back(m_newClusUCol[i].get());
  //for(int i=0; i<m_newClusVCol.size(); i++) m_HFClusVInTower.push_back(m_newClusVCol[i].get());


cout<<"Check tower ID: ";
for(int i=0; i<it->getTowerID().size(); i++) printf("[%d, %d, %d], ", it->getTowerID()[i][0], it->getTowerID()[i][1], it->getTowerID()[i][2]);
cout<<endl;
  printf("    In Tower [%d, %d, %d], ", it->getTowerID()[0][0], it->getTowerID()[0][1], it->getTowerID()[0][2] );
  printf("    HalfCluster size: (%d, %d) \n", m_HFClusUInTower.size(), m_HFClusVInTower.size() );
  cout<<"    Loop print HalfClusterU: "<<endl;
  for(int icl=0; icl<m_HFClusUInTower.size(); icl++){
    cout<<"      In HFClusU #"<<icl<<": shower size = "<<m_HFClusUInTower[icl]->getCluster().size()<<", En = "<<m_HFClusUInTower[icl]->getEnergy()<<", linked MC size "<<m_HFClusUInTower[icl]->getLinkedMCP().size();
    for(int imc=0; imc<m_HFClusUInTower[icl]->getLinkedMCP().size(); imc++) printf("[%d, %.3f], ", m_HFClusUInTower[icl]->getLinkedMCP()[imc].first.getPDG(), m_HFClusUInTower[icl]->getLinkedMCP()[imc].second);
    printf(", Position (%.3f, %.3f, %.3f), address %p ",m_HFClusUInTower[icl]->getPos().x(), m_HFClusUInTower[icl]->getPos().y(), m_HFClusUInTower[icl]->getPos().z(), m_HFClusUInTower[icl]);
    printf(", cousin size %d, address: ", m_HFClusUInTower[icl]->getHalfClusterCol("CousinCluster").size());
    for(int ics=0; ics<m_HFClusUInTower[icl]->getHalfClusterCol("CousinCluster").size(); ics++) printf("%p, ", m_HFClusUInTower[icl]->getHalfClusterCol("CousinCluster")[ics]);
    printf(", track size %d, address: ", m_HFClusUInTower[icl]->getAssociatedTracks().size());
    for(int itrk=0; itrk<m_HFClusUInTower[icl]->getAssociatedTracks().size(); itrk++) printf("%p, ", m_HFClusUInTower[icl]->getAssociatedTracks()[itrk]);
    cout<<endl;
    for(auto ish : m_HFClusUInTower[icl]->getCluster()){
      printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), leading MCP [%d, %.3f], Nbars %d, NSeed %d, Address %p \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getLeadingMCP().getPDG(), ish->getLeadingMCPweight(),  ish->getBars().size(),  ish->getNseeds(), ish );
    }
  }
  cout<<endl;

  cout<<"    Loop print HalfClusterV: "<<endl;
  for(int icl=0; icl<m_HFClusVInTower.size(); icl++){
    cout<<"      In HFClusV #"<<icl<<": shower size = "<<m_HFClusVInTower[icl]->getCluster().size()<<", En = "<<m_HFClusVInTower[icl]->getEnergy()<<", linked MC size "<<m_HFClusVInTower[icl]->getLinkedMCP().size();
    for(int imc=0; imc<m_HFClusVInTower[icl]->getLinkedMCP().size(); imc++) printf("[%d, %.3f], ", m_HFClusVInTower[icl]->getLinkedMCP()[imc].first.getPDG(), m_HFClusVInTower[icl]->getLinkedMCP()[imc].second);
    printf(", Position (%.3f, %.3f, %.3f), address %p ",m_HFClusVInTower[icl]->getPos().x(), m_HFClusVInTower[icl]->getPos().y(), m_HFClusVInTower[icl]->getPos().z(), m_HFClusVInTower[icl]);
    printf(", cousin size %d, address: ", m_HFClusVInTower[icl]->getHalfClusterCol("CousinCluster").size());
    for(int ics=0; ics<m_HFClusVInTower[icl]->getHalfClusterCol("CousinCluster").size(); ics++) printf("%p, ", m_HFClusVInTower[icl]->getHalfClusterCol("CousinCluster")[ics]);
    printf(", track size %d, address: ", m_HFClusVInTower[icl]->getAssociatedTracks().size());
    for(int itrk=0; itrk<m_HFClusVInTower[icl]->getAssociatedTracks().size(); itrk++) printf("%p, ", m_HFClusVInTower[icl]->getAssociatedTracks()[itrk]);
    cout<<endl;
    for(auto ish : m_HFClusVInTower[icl]->getCluster()){
      printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), leading MCP [%d, %.3f], Nbars %d, NSeed %d, Address %p \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getLeadingMCP().getPDG(), ish->getLeadingMCPweight(), ish->getBars().size(), ish->getNseeds(), ish );
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
                                                         std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>& m_towers )
{

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
        m_bkCol.map_HalfCluster["bkHalfCluster"].push_back( tmp_clus );
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
        m_bkCol.map_HalfCluster["bkHalfCluster"].push_back( tmp_clus );
      }
      p_shower = nullptr;

    }

    //Connect cousins
//cout<<"  LongiClusVMap size: "<<tmp_LongiClusMaps.size()<<endl;
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
      m_bkCol.map_2DCluster["bk2DCluster"].push_back( tmp_block );
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
//printf("  In tower: [%d, %d, %d] \n", m_towerID[0], m_towerID[1], m_towerID[2]);
    //Check cousin clusters:
    std::vector<PandoraPlus::CaloHalfCluster*> m_HFClusUInTower = map_HalfClusterU[m_towerID];
    for(auto &m_HFclus : m_HFClusUInTower){
      std::vector<const CaloHalfCluster*> tmp_delClus; tmp_delClus.clear();
//printf("    Check the cousin of HFClus %p: cousin size %d \n",m_HFclus, m_HFclus->getHalfClusterCol("CousinCluster").size());
      for(int ics=0; ics<m_HFclus->getHalfClusterCol("CousinCluster").size(); ics++){
        std::vector<int> tmp_towerID = m_HFclus->getHalfClusterCol("CousinCluster")[ics]->getTowerID()[0];
//printf("      Cousin #%d: address %p, it's in tower [%d, %d, %d]. \n", ics,  m_HFclus->getHalfClusterCol("CousinCluster")[ics],  tmp_towerID[0], tmp_towerID[1], tmp_towerID[2]);
//if(m_HFclus->getHalfClusterCol("CousinCluster")[ics]->getTowerID().size()!=1) cout<<"ERROR: cousin cluster covers >1 towers. Check here! "<<endl;

        if( map_2DCluster.find( tmp_towerID )==map_2DCluster.end() )
          tmp_delClus.push_back( m_HFclus->getHalfClusterCol("CousinCluster")[ics] );
      }
//cout<<"Need to delete "<<tmp_delClus.size()<<" cousins"<<endl;
      for(int ics=0; ics<tmp_delClus.size(); ics++) m_HFclus->deleteCousinCluster( tmp_delClus[ics] );
    }

    std::vector<PandoraPlus::CaloHalfCluster*> m_HFClusVInTower = map_HalfClusterV[m_towerID];
    for(auto &m_HFclus : m_HFClusVInTower){
      std::vector<const CaloHalfCluster*> tmp_delClus; tmp_delClus.clear();
//printf("    Check the cousin of HFClus %p: cousin size %d \n",m_HFclus, m_HFclus->getHalfClusterCol("CousinCluster").size());
      for(int ics=0; ics<m_HFclus->getHalfClusterCol("CousinCluster").size(); ics++){
        std::vector<int> tmp_towerID = m_HFclus->getHalfClusterCol("CousinCluster")[ics]->getTowerID()[0];
//printf("      Cousin #%d: address %p, it's in tower [%d, %d, %d]. \n", ics,  m_HFclus->getHalfClusterCol("CousinCluster")[ics],  tmp_towerID[0], tmp_towerID[1], tmp_towerID[2]);
//if(m_HFclus->getHalfClusterCol("CousinCluster")[ics]->getTowerID().size()!=1) cout<<"ERROR: cousin cluster covers >1 towers. Check here! "<<endl;

        if( map_2DCluster.find( tmp_towerID )==map_2DCluster.end() )
          tmp_delClus.push_back( m_HFclus->getHalfClusterCol("CousinCluster")[ics] );
      }
//cout<<"Need to delete "<<tmp_delClus.size()<<" cousins"<<endl;
      for(int ics=0; ics<tmp_delClus.size(); ics++) m_HFclus->deleteCousinCluster( tmp_delClus[ics] );
    }

    //Convert to const
    std::vector<const PandoraPlus::CaloHalfCluster*> const_HFClusU; const_HFClusU.clear();
    std::vector<const PandoraPlus::CaloHalfCluster*> const_HFClusV; const_HFClusV.clear();
    for(int ics=0; ics<m_HFClusUInTower.size(); ics++){ m_HFClusUInTower[ics]->getLinkedMCPfromUnit(); const_HFClusU.push_back(m_HFClusUInTower[ics]); }
    for(int ics=0; ics<m_HFClusVInTower.size(); ics++){ m_HFClusVInTower[ics]->getLinkedMCPfromUnit(); const_HFClusV.push_back(m_HFClusVInTower[ics]); }

    std::shared_ptr<PandoraPlus::Calo3DCluster> m_tower = std::make_shared<PandoraPlus::Calo3DCluster>();
//printf("Creating tower: [%d, %d, %d] \n", m_towerID[0], m_towerID[1], m_towerID[2]);
    m_tower->addTowerID( m_towerID );
    for(int i2d=0; i2d<map_2DCluster[m_towerID].size(); i2d++) m_tower->addUnit(map_2DCluster[m_towerID][i2d]);
    m_tower->setHalfClusters( settings.map_stringPars["OutputClusName"]+"U", const_HFClusU,
                              settings.map_stringPars["OutputClusName"]+"V", const_HFClusV );
    m_towers.push_back(m_tower);
  }

  return StatusCode::SUCCESS;
}




/*
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
    m_tower->setHalfClusters( settings.map_stringPars["OutputClusName"]+"U", const_HFClusU,
                              settings.map_stringPars["OutputClusName"]+"V", const_HFClusV );
    m_towers.push_back(m_tower);
  }

  return StatusCode::SUCCESS;
}
*/

#endif
