#ifndef _TRUTHCLUS_ALG_C
#define _TRUTHCLUS_ALG_C

#include "Algorithm/TruthClusteringAlg.h"

StatusCode TruthClusteringAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  if(settings.map_boolPars.find("UseSplit")==settings.map_boolPars.end()) 
    settings.map_boolPars["UseSplit"] = 1;
  if(settings.map_boolPars.find("DoECALClustering")==settings.map_boolPars.end()) 
    settings.map_boolPars["DoECALClustering"] = 1;
  if(settings.map_boolPars.find("DoHCALClustering")==settings.map_boolPars.end()) 
    settings.map_boolPars["DoHCALClustering"] = 1;

  if(settings.map_stringPars.find("InputECALBars")==settings.map_stringPars.end()) 
    settings.map_stringPars["InputECALBars"] = "BarCol";
  if(settings.map_stringPars.find("InputHCALHits")==settings.map_stringPars.end()) 
    settings.map_stringPars["InputHCALHits"] = "HCALBarrel";
  if(settings.map_stringPars.find("OutputECAL1DClusters")==settings.map_stringPars.end()) 
    settings.map_stringPars["OutputECAL1DClusters"] = "TruthCluster1DCol";
  if(settings.map_stringPars.find("OutputECALHalfClusters")==settings.map_stringPars.end()) 
    settings.map_stringPars["OutputECALHalfClusters"] = "TruthESCluster";
  if(settings.map_stringPars.find("OutputECALTower")==settings.map_stringPars.end()) 
    settings.map_stringPars["OutputECALTower"] = "TruthESTower";
  if(settings.map_stringPars.find("OutputHCALClusters")==settings.map_stringPars.end()) 
    settings.map_stringPars["OutputHCALClusters"] = "TruthHcalCluster";

  return StatusCode::SUCCESS;
};

StatusCode TruthClusteringAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  m_TrackCol.clear();
  m_bars.clear();
  m_1dclusterUCol.clear();
  m_1dclusterVCol.clear();
  m_halfclusterU.clear();
  m_halfclusterV.clear();
  m_towers.clear();

  m_hits.clear();
  m_clusters.clear();
  m_bkCol.Clear();

  m_TrackCol = m_datacol.TrackCol; 
  m_bars = m_datacol.map_BarCol[settings.map_stringPars["InputECALBars"]];
  m_hits = m_datacol.map_CaloHit[settings.map_stringPars["InputHCALHits"]];

  return StatusCode::SUCCESS;
};

StatusCode TruthClusteringAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  //ECAL clustering
  if(settings.map_boolPars["DoECALClustering"]){
cout<<"Input bar size: "<<m_bars.size()<<endl;

    std::map<edm4hep::MCParticle, std::vector<std::shared_ptr<PandoraPlus::CaloUnit>> > map_barCol; map_barCol.clear();
    for(int ibar=0; ibar<m_bars.size(); ibar++){
      if(settings.map_boolPars["UseSplit"]){
        std::vector< std::pair<edm4hep::MCParticle, float> > m_linkVec = m_bars[ibar]->getLinkedMCP();
//printf("  Bar #%d: link size %d: ", ibar, m_linkVec.size());
//for(int il=0; il<m_linkVec.size(); il++) printf("(%d, %.2f), ", m_linkVec[il].first.getPDG(), m_linkVec[il].second);
//cout<<endl;

        if(m_linkVec.size()==0){
          std::cout<<"ERROR: No truth info in EcalBar #"<<ibar<<std::endl;
          continue;
        }
        if(m_linkVec.size()==1){
          edm4hep::MCParticle mcp = m_linkVec[0].first;
          map_barCol[mcp].push_back(m_bars[ibar]);
          continue;
        }
        for(int ilink=0; ilink<m_linkVec.size(); ilink++){
          auto tmp_newbar = m_bars[ibar]->Clone();
          tmp_newbar->setQ(tmp_newbar->getQ1()*m_linkVec[ilink].second, tmp_newbar->getQ2()*m_linkVec[ilink].second);
          map_barCol[m_linkVec[ilink].first].push_back(tmp_newbar);
          tmp_newbar->setLinkedMCP( std::vector<std::pair<edm4hep::MCParticle, float>>{std::make_pair(m_linkVec[ilink].first, 1.)} );
          m_datacol.map_BarCol["bkBars"].push_back(tmp_newbar); 
        }
      }
      else{
        edm4hep::MCParticle mcp = m_bars[ibar]->getLeadingMCP();
        map_barCol[mcp].push_back(m_bars[ibar]);
      }
    }
cout<<"truth map size "<<map_barCol.size()<<endl;
for(auto& iter: map_barCol){ 
  double totE = 0.;
  for(int ibar=0; ibar<iter.second.size(); ibar++)
    totE += iter.second[ibar]->getEnergy();
  printf("  Truth MCP pdgid %d, bar size %d, totE %.3f \n", iter.first.getPDG(), iter.second.size(), totE);
}   


    for(auto& iter: map_barCol){
      std::map<int, std::vector<std::shared_ptr<PandoraPlus::CaloUnit>> > m_orderedBars; m_orderedBars.clear();
      for(int ibar=0; ibar<iter.second.size(); ibar++)
        m_orderedBars[iter.second[ibar]->getDlayer()].push_back(iter.second[ibar]);
//cout<<"  MCP "<<iter.first.getPDG()<<" covers "<<m_orderedBars.size()<<"layers "<<endl;

      std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>> m_new1DClusUCol; m_new1DClusUCol.clear();
      std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>> m_new1DClusVCol; m_new1DClusVCol.clear();
   
      for(auto& iter_bar: m_orderedBars){
        std::shared_ptr<PandoraPlus::Calo1DCluster> tmp_new1dclusterU = std::make_shared<PandoraPlus::Calo1DCluster>();
        std::shared_ptr<PandoraPlus::Calo1DCluster> tmp_new1dclusterV = std::make_shared<PandoraPlus::Calo1DCluster>();
        for(int ibar=0; ibar<iter_bar.second.size(); ibar++){ 
          if(iter_bar.second[ibar]->getSlayer()==0) tmp_new1dclusterU->addUnit( iter_bar.second[ibar].get() );
          if(iter_bar.second[ibar]->getSlayer()==1) tmp_new1dclusterV->addUnit( iter_bar.second[ibar].get() );
        }
        if(tmp_new1dclusterU  && tmp_new1dclusterU->getBars().size()>0){
          tmp_new1dclusterU->getLinkedMCPfromUnit();
          tmp_new1dclusterU->setSeed();
          tmp_new1dclusterU->setIDInfo();
          m_new1DClusUCol.push_back(tmp_new1dclusterU);
          m_datacol.map_1DCluster["bk1DCluster"].push_back(tmp_new1dclusterU);
        }
        if(tmp_new1dclusterV  && tmp_new1dclusterV->getBars().size()>0){
          tmp_new1dclusterV->getLinkedMCPfromUnit();
          tmp_new1dclusterV->setSeed();
          tmp_new1dclusterV->setIDInfo();
          m_new1DClusVCol.push_back(tmp_new1dclusterV);
          m_datacol.map_1DCluster["bk1DCluster"].push_back(tmp_new1dclusterV);
        }
      }
//cout<<"  1D cluster size "<<m_new1DClusUCol.size()<<", "<<m_new1DClusVCol.size()<<endl;
//for(int icl=0; icl<m_new1DClusUCol.size(); icl++){
//  printf("    Dlayer %d, Slayer %d, En %.3f, bar size %d, seed size %d, MCP link size %d \n", m_new1DClusUCol[icl]->getDlayer(), m_new1DClusUCol[icl]->getSlayer(), m_new1DClusUCol[icl]->getEnergy(), m_new1DClusUCol[icl]->getBars().size(), m_new1DClusUCol[icl]->getNseeds(), m_new1DClusUCol[icl]->getLinkedMCP().size());
//}
//for(int icl=0; icl<m_new1DClusVCol.size(); icl++){
//  printf("    Dlayer %d, Slayer %d, En %.3f, bar size %d, seed size %d, MCP link size %d \n", m_new1DClusVCol[icl]->getDlayer(), m_new1DClusVCol[icl]->getSlayer(), m_new1DClusVCol[icl]->getEnergy(), m_new1DClusVCol[icl]->getBars().size(), m_new1DClusVCol[icl]->getNseeds(), m_new1DClusVCol[icl]->getLinkedMCP().size());
//}

      m_1dclusterUCol.insert(m_1dclusterUCol.end(), m_new1DClusUCol.begin(), m_new1DClusUCol.end());
      m_1dclusterVCol.insert(m_1dclusterVCol.end(), m_new1DClusVCol.begin(), m_new1DClusVCol.end());

      std::shared_ptr<PandoraPlus::CaloHalfCluster> m_newHFClusterU = std::make_shared<PandoraPlus::CaloHalfCluster>();
      std::shared_ptr<PandoraPlus::CaloHalfCluster> m_newHFClusterV = std::make_shared<PandoraPlus::CaloHalfCluster>();
      for(int i1d=0; i1d<m_new1DClusUCol.size(); i1d++)
        m_newHFClusterU->addUnit(m_new1DClusUCol[i1d].get());
      for(int i1d=0; i1d<m_new1DClusVCol.size(); i1d++)
        m_newHFClusterV->addUnit(m_new1DClusVCol[i1d].get());
      
      m_newHFClusterU->getLinkedMCPfromUnit();
      m_newHFClusterV->getLinkedMCPfromUnit();
      m_halfclusterU.push_back(m_newHFClusterU);
      m_halfclusterV.push_back(m_newHFClusterV);
      m_datacol.map_HalfCluster["bkHalfCluster"].push_back(m_newHFClusterU);
      m_datacol.map_HalfCluster["bkHalfCluster"].push_back(m_newHFClusterV);
    }

    //Track match
    for(int itrk=0; itrk<m_TrackCol.size(); itrk++){
      edm4hep::MCParticle mcp_trk = m_TrackCol[itrk]->getLeadingMCP();
      for(int ihf=0; ihf<m_halfclusterU.size(); ihf++){
        if(m_halfclusterU[ihf]->getLeadingMCP()==mcp_trk){
          m_halfclusterU[ihf]->addAssociatedTrack(m_TrackCol[itrk].get());
          m_TrackCol[itrk]->addAssociatedHalfClusterU(m_halfclusterU[ihf].get());
        }
      }
      for(int ihf=0; ihf<m_halfclusterV.size(); ihf++){
        if(m_halfclusterV[ihf]->getLeadingMCP()==mcp_trk){
          m_halfclusterV[ihf]->addAssociatedTrack(m_TrackCol[itrk].get());
          m_TrackCol[itrk]->addAssociatedHalfClusterV(m_halfclusterV[ihf].get());
        }
      }
    }

/*
cout<<"Check tower ID: ";
cout<<endl;
  printf("    HalfCluster size: (%d, %d) \n", m_halfclusterU.size(), m_halfclusterV.size() );
  cout<<"    Loop print HalfClusterU: "<<endl;
  for(int icl=0; icl<m_halfclusterU.size(); icl++){
    cout<<"      In HFClusU #"<<icl<<": shower size = "<<m_halfclusterU[icl]->getCluster().size()<<", En = "<<m_halfclusterU[icl]->getEnergy();
    printf(", Position (%.3f, %.3f, %.3f), address %p ",m_halfclusterU[icl]->getPos().x(), m_halfclusterU[icl]->getPos().y(), m_halfclusterU[icl]->getPos().z(), m_halfclusterU[icl]);
    printf(", cousin size %d, address: ", m_halfclusterU[icl]->getHalfClusterCol("CousinCluster").size());
    for(int ics=0; ics<m_halfclusterU[icl]->getHalfClusterCol("CousinCluster").size(); ics++) printf("%p, ", m_halfclusterU[icl]->getHalfClusterCol("CousinCluster")[ics]);
    printf(", track size %d, address: ", m_halfclusterU[icl]->getAssociatedTracks().size());
    for(int itrk=0; itrk<m_halfclusterU[icl]->getAssociatedTracks().size(); itrk++) printf("%p, ", m_halfclusterU[icl]->getAssociatedTracks()[itrk]);
    printf(", MCP link size %d, pid, Pz and weight: ", m_halfclusterU[icl]->getLinkedMCP().size());
    for(int imc=0; imc<m_halfclusterU[icl]->getLinkedMCP().size(); imc++) printf("(%d, %.3f, %.3f), ", m_halfclusterU[icl]->getLinkedMCP()[imc].first.getPDG(), m_halfclusterU[icl]->getLinkedMCP()[imc].first.getMomentum().z, m_halfclusterU[icl]->getLinkedMCP()[imc].second);
    cout<<endl;
    for(auto ish : m_halfclusterU[icl]->getCluster()){
      printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, Address %p \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(),  ish->getNseeds(), ish );
    }
  }
  cout<<endl;

  cout<<"    Loop print HalfClusterV: "<<endl;
  for(int icl=0; icl<m_halfclusterV.size(); icl++){
    cout<<"      In HFClusV #"<<icl<<": shower size = "<<m_halfclusterV[icl]->getCluster().size()<<", En = "<<m_halfclusterV[icl]->getEnergy();
    printf(", Position (%.3f, %.3f, %.3f), address %p ",m_halfclusterV[icl]->getPos().x(), m_halfclusterV[icl]->getPos().y(), m_halfclusterV[icl]->getPos().z(), m_halfclusterV[icl]);
    printf(", cousin size %d, address: ", m_halfclusterV[icl]->getHalfClusterCol("CousinCluster").size());
    for(int ics=0; ics<m_halfclusterV[icl]->getHalfClusterCol("CousinCluster").size(); ics++) printf("%p, ", m_halfclusterV[icl]->getHalfClusterCol("CousinCluster")[ics]);
    printf(", track size %d, address: ", m_halfclusterV[icl]->getAssociatedTracks().size());
    for(int itrk=0; itrk<m_halfclusterV[icl]->getAssociatedTracks().size(); itrk++) printf("%p, ", m_halfclusterV[icl]->getAssociatedTracks()[itrk]);
    printf(", MCP link size %d, pid, Pz and weight: ", m_halfclusterV[icl]->getLinkedMCP().size());
    for(int imc=0; imc<m_halfclusterV[icl]->getLinkedMCP().size(); imc++) printf("(%d, %.3f, %.3f), ", m_halfclusterV[icl]->getLinkedMCP()[imc].first.getPDG(), m_halfclusterV[icl]->getLinkedMCP()[imc].first.getMomentum().z, m_halfclusterV[icl]->getLinkedMCP()[imc].second);
    cout<<endl;
    for(auto ish : m_halfclusterV[icl]->getCluster()){
      printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, Address %p \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(), ish->getNseeds(), ish );
    }
  }
  cout<<endl;
*/

    //Create tower
    HalfClusterToTowers(m_halfclusterU, m_halfclusterV, m_towers);
cout<<"  Ecal halfcluster size "<<m_halfclusterU.size()<<", "<<m_halfclusterV.size()<<", tower size "<<m_towers.size()<<endl;

    m_datacol.map_1DCluster[settings.map_stringPars["OutputECAL1DClusters"]+"U"] = m_1dclusterUCol;
    m_datacol.map_1DCluster[settings.map_stringPars["OutputECAL1DClusters"]+"V"] = m_1dclusterVCol;
    m_datacol.map_HalfCluster[settings.map_stringPars["OutputECALHalfClusters"]+"U"] = m_halfclusterU;
    m_datacol.map_HalfCluster[settings.map_stringPars["OutputECALHalfClusters"]+"V"] = m_halfclusterV;  
    m_datacol.map_CaloCluster[settings.map_stringPars["OutputECALTower"]] = m_towers;
  }




  //HCAL clustering
  if(settings.map_boolPars["DoHCALClustering"]){
cout<<"Input HCAL hit size "<<m_hits.size()<<endl;
    std::map<edm4hep::MCParticle, std::vector<std::shared_ptr<PandoraPlus::CaloHit>> > map_hitCol; map_hitCol.clear();
    for(int ihit=0; ihit<m_hits.size(); ihit++){
      if(settings.map_boolPars["UseSplit"]){
        std::vector< std::pair<edm4hep::MCParticle, float> > m_linkVec = m_hits[ihit]->getLinkedMCP();
//printf("  Hit #%d: link size %d: ", ihit, m_linkVec.size());
//for(int il=0; il<m_linkVec.size(); il++) printf("(%d, %.2f), ", m_linkVec[il].first.getPDG(), m_linkVec[il].second);
//cout<<endl;

        if(m_linkVec.size()==0){
          std::cout<<"ERROR: No truth info in hit #"<<ihit<<std::endl;
          continue;
        }
        if(m_linkVec.size()==1){
          edm4hep::MCParticle mcp = m_linkVec[0].first;
          map_hitCol[mcp].push_back(m_hits[ihit]);
          continue;
        }
   
        for(int ilink=0; ilink<m_linkVec.size(); ilink++){
          auto tmp_newhit = m_hits[ihit]->Clone();
          tmp_newhit->setEnergy(tmp_newhit->getEnergy()*m_linkVec[ilink].second);
          tmp_newhit->setLinkedMCP(std::vector<std::pair<edm4hep::MCParticle, float>>{std::make_pair(m_linkVec[ilink].first, 1.)});
          map_hitCol[m_linkVec[ilink].first].push_back(tmp_newhit);
          m_datacol.map_CaloHit["bkHit"].push_back(tmp_newhit);
        }
      }
   
      else{
        edm4hep::MCParticle mcp = m_hits[ihit]->getLeadingMCP();
        map_hitCol[mcp].push_back(m_hits[ihit]);
      }
   
    }//End loop hits

cout<<"truth map size "<<map_hitCol.size()<<endl;
for(auto& iter: map_hitCol){
  double totE = 0.;
  for(int ibar=0; ibar<iter.second.size(); ibar++)
    totE += iter.second[ibar]->getEnergy();
  printf("  Truth MCP pdgid %d, bar size %d, totE %.3f \n", iter.first.getPDG(), iter.second.size(), totE);
} 
  
    for(auto& iter: map_hitCol){
      std::shared_ptr<PandoraPlus::Calo3DCluster> m_newCluster = std::make_shared<PandoraPlus::Calo3DCluster>();
      for(int ihit=0; ihit<iter.second.size(); ihit++) m_newCluster->addHit( iter.second[ihit].get() );    
      m_newCluster->getLinkedMCPfromHit();
      m_clusters.push_back(m_newCluster);
      m_datacol.map_CaloCluster["bk3DCluster"].push_back(m_newCluster);
    }

    for(int icl=0; icl<m_clusters.size(); icl++){
      for(int itrk=0; itrk<m_TrackCol.size(); itrk++){
        edm4hep::MCParticle mcp_trk = m_TrackCol[itrk]->getLeadingMCP();
        if(m_clusters[icl]->getLeadingMCP()==mcp_trk){
          m_clusters[icl]->addAssociatedTrack(m_TrackCol[itrk].get());
          break;
        }
      }
    }
cout<<"  HCAL cluster size "<<m_clusters.size()<<endl;
for(int ic=0; ic<m_clusters.size(); ic++){


}
    m_datacol.map_CaloCluster[settings.map_stringPars["OutputHCALClusters"]] = m_clusters;
  }

  m_datacol.map_1DCluster["bk1DCluster"].insert( m_datacol.map_1DCluster["bk1DCluster"].end(), m_bkCol.map_1DCluster["bk1DCluster"].begin(), m_bkCol.map_1DCluster["bk1DCluster"].end() );
  m_datacol.map_2DCluster["bk2DCluster"].insert( m_datacol.map_2DCluster["bk2DCluster"].end(), m_bkCol.map_2DCluster["bk2DCluster"].begin(), m_bkCol.map_2DCluster["bk2DCluster"].end() );
  m_datacol.map_HalfCluster["bkHalfCluster"].insert( m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_bkCol.map_HalfCluster["bkHalfCluster"].begin(), m_bkCol.map_HalfCluster["bkHalfCluster"].end() );


  return StatusCode::SUCCESS;
};

StatusCode TruthClusteringAlg::ClearAlgorithm(){
  m_bars.clear();
  m_1dclusterUCol.clear();
  m_halfclusterU.clear();
  m_halfclusterV.clear();

  m_hits.clear();
  m_clusters.clear();
  m_towers.clear();

  m_bkCol.Clear();
  return StatusCode::SUCCESS;
};

StatusCode TruthClusteringAlg::HalfClusterToTowers( std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& m_halfClusU,
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

//cout<<"  tower size: "<<map_barCol.size()<<endl;
//cout<<"  Print tower info "<<endl;
//for(auto& iter: map_barCol){
//printf("    TowerID [%d, %d, %d], bar size %d \n", iter.first[0], iter.first[1], iter.first[2], iter.second.size());
//}

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
//printf("    TowerID [%d, %d, %d], MCP size %d, bar size %d \n", itower.first[0], itower.first[1], itower.first[2], map_matchBar.size(), itower.second.size());

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

//cout<<"    Layer size "<<m_orderedBars.size()<<endl;

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
//printf("    In tower [%d, %d, %d] HFClus size (%d, %d) \n", itower.first[0], itower.first[1], itower.first[2], m_newHFClusUCol.size(), m_newHFClusVCol.size());
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

/*
cout<<"  After splitting: tower size "<<m_towers.size()<<". Print Tower: "<<endl;
for(auto it : m_towers){
  std::vector<const CaloHalfCluster*> m_HFClusUInTower = it->getHalfClusterUCol(settings.map_stringPars["OutputECALHalfClusters"]+"U");
  std::vector<const CaloHalfCluster*> m_HFClusVInTower = it->getHalfClusterVCol(settings.map_stringPars["OutputECALHalfClusters"]+"V");

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
    printf(", MCP link size %d, pid, Pz and weight: ", m_HFClusUInTower[icl]->getLinkedMCP().size());
    for(int imc=0; imc<m_HFClusUInTower[icl]->getLinkedMCP().size(); imc++) printf("(%d, %.3f, %.3f), ", m_HFClusUInTower[icl]->getLinkedMCP()[imc].first.getPDG(), m_HFClusUInTower[icl]->getLinkedMCP()[imc].first.getMomentum().z, m_HFClusUInTower[icl]->getLinkedMCP()[imc].second);
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
    printf(", MCP link size %d, pid, Pz and weight: ", m_HFClusVInTower[icl]->getLinkedMCP().size());
    for(int imc=0; imc<m_HFClusVInTower[icl]->getLinkedMCP().size(); imc++) printf("(%d, %.3f, %.3f), ", m_HFClusVInTower[icl]->getLinkedMCP()[imc].first.getPDG(), m_HFClusVInTower[icl]->getLinkedMCP()[imc].first.getMomentum().z, m_HFClusVInTower[icl]->getLinkedMCP()[imc].second);
    cout<<endl;
    for(auto ish : m_HFClusVInTower[icl]->getCluster()){
      printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, Address %p \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(), ish->getNseeds(), ish );
    }
  }
  cout<<endl;
}
*/


  return StatusCode::SUCCESS;
}


#endif
