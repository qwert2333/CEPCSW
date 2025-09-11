#ifndef _ENERGYSPLITTING_ALG_C
#define _ENERGYSPLITTING_ALG_C

#include "Algorithm/EnergySplittingAlg.h"

StatusCode EnergySplittingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  if(settings.map_floatPars.find("th_split")==settings.map_floatPars.end()) settings.map_floatPars["th_split"] = -1;
  //if(settings.map_floatPars.find("Eth_Seed")==settings.map_floatPars.end()) settings.map_floatPars["Eth_Seed"] = 0.005;
  //if(settings.map_floatPars.find("Eth_ShowerAbs")==settings.map_floatPars.end()) settings.map_floatPars["Eth_ShowerAbs"] = 0.005;
  if(settings.map_floatPars.find("Eth_unit")==settings.map_floatPars.end())  settings.map_floatPars["Eth_unit"] = 0.001;
  if(settings.map_floatPars.find("Eth_HFClus")==settings.map_floatPars.end())        settings.map_floatPars["Eth_HFClus"] = 0.05;
  if(settings.map_intPars.find("th_Nhit")==settings.map_intPars.end()) settings.map_intPars["th_Nhit"] = 2;
  if(settings.map_boolPars.find("CompactHFCluster")==settings.map_boolPars.end()) settings.map_boolPars["CompactHFCluster"] = true;
  if(settings.map_stringPars.find("ReadinAxisName")==settings.map_stringPars.end()) settings.map_stringPars["ReadinAxisName"] = "MergedAxis";
  if(settings.map_stringPars.find("OutputClusName")==settings.map_stringPars.end()) settings.map_stringPars["OutputClusName"] = "ESHalfCluster";
  if(settings.map_stringPars.find("OutputTowerName")==settings.map_stringPars.end()) settings.map_stringPars["OutputTowerName"] = "ESTower";

  return StatusCode::SUCCESS;
};


StatusCode EnergySplittingAlg::Initialize( CyberDataCol& m_datacol ){
  m_axisUCol.clear();
  m_axisVCol.clear();
  m_newClusUCol.clear();
  m_newClusVCol.clear();
  m_1dShowerUCol.clear();
  m_1dShowerVCol.clear();
  m_towerCol.clear(); 
  m_bkCol.Clear();

  p_HalfClusterU.clear();
  p_HalfClusterV.clear(); 

  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColU"].size(); ih++)
    p_HalfClusterU.push_back( m_datacol.map_HalfCluster["HalfClusterColU"][ih].get() );
  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColV"].size(); ih++)
    p_HalfClusterV.push_back( m_datacol.map_HalfCluster["HalfClusterColV"][ih].get() );
  for(int ih=0; ih<m_datacol.map_HalfCluster["emptyHalfClusterU"].size(); ih++){
    p_emptyHalfClusterU.push_back( m_datacol.map_HalfCluster["emptyHalfClusterU"][ih].get() );
  }
  for(int ih=0; ih<m_datacol.map_HalfCluster["emptyHalfClusterV"].size(); ih++)
    p_emptyHalfClusterV.push_back( m_datacol.map_HalfCluster["emptyHalfClusterV"][ih].get() );

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::RunAlgorithm( CyberDataCol& m_datacol ){

  //Input: HalfCluster (with 1D clusters and HoughAxis)
  //Output: Create Towers for the matching. 

  if( p_HalfClusterU.size() + p_HalfClusterV.size()==0 ) {
    if( p_emptyHalfClusterU.size() + p_emptyHalfClusterV.size() == 0 ){
      std::cout<<"EnergySplittingAlg: No HalfCluster input"<<std::endl;
      return StatusCode::SUCCESS;
    }
    m_towerCol.clear();
    std::vector<Cyber::CaloHalfCluster*> tmp_newClusUCol, tmp_newClusVCol;
    HalfClusterToTowers( tmp_newClusUCol, tmp_newClusVCol, p_emptyHalfClusterU, p_emptyHalfClusterV, m_towerCol );

    m_datacol.map_CaloCluster[settings.map_stringPars["OutputTowerName"]] = m_towerCol;

    m_datacol.map_BarCol["bkBar"].insert( m_datacol.map_BarCol["bkBar"].end(), m_bkCol.map_BarCol["bkBar"].begin(), m_bkCol.map_BarCol["bkBar"].end() );
    m_datacol.map_1DCluster["bk1DCluster"].insert( m_datacol.map_1DCluster["bk1DCluster"].end(), m_bkCol.map_1DCluster["bk1DCluster"].begin(), m_bkCol.map_1DCluster["bk1DCluster"].end() );
    m_datacol.map_2DCluster["bk2DCluster"].insert( m_datacol.map_2DCluster["bk2DCluster"].end(), m_bkCol.map_2DCluster["bk2DCluster"].begin(), m_bkCol.map_2DCluster["bk2DCluster"].end() );
    m_datacol.map_HalfCluster["bkHalfCluster"].insert( m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_bkCol.map_HalfCluster["bkHalfCluster"].begin(), m_bkCol.map_HalfCluster["bkHalfCluster"].end() );

    return StatusCode::SUCCESS;
  }

/*
float tmp_totE_U = 0;
for(int ih=0; ih<p_HalfClusterU.size(); ih++) tmp_totE_U += p_HalfClusterU[ih]->getEnergy();
float tmp_totE_V = 0;
for(int ih=0; ih<p_HalfClusterV.size(); ih++) tmp_totE_V += p_HalfClusterV[ih]->getEnergy();
cout<<"Total energy in HFCluster: "<<tmp_totE_U<<'\t'<<tmp_totE_V<<endl;
cout<<"Print readin axesU"<<endl;
for(int ih=0; ih<p_HalfClusterU.size(); ih++){
  cout<<"  HalfCluster #"<<ih<<": cluster En "<<p_HalfClusterU[ih]->getEnergy()<<endl;
  auto tmp_axisUCol = p_HalfClusterU[ih]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"]);
  for(int icl=0; icl<tmp_axisUCol.size(); icl++){
    cout<<"  Axis #"<<icl<<": cluster En "<<tmp_axisUCol[icl]->getEnergy()<<", type "<<tmp_axisUCol[icl]->getType()<<endl;
    for(auto ish : tmp_axisUCol[icl]->getCluster()){
      printf("    Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Address \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy());
    }
  }
}
cout<<"Print readin axesV"<<endl;
for(int ih=0; ih<p_HalfClusterV.size(); ih++){
  cout<<"  HalfCluster #"<<ih<<": cluster En "<<p_HalfClusterV[ih]->getEnergy()<<endl;
  auto tmp_axisVCol = p_HalfClusterV[ih]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"]);
  for(int icl=0; icl<tmp_axisVCol.size(); icl++){
    cout<<"  Axis #"<<icl<<": cluster En "<<tmp_axisVCol[icl]->getEnergy()<<", type "<<tmp_axisVCol[icl]->getType()<<endl;
    for(auto ish : tmp_axisVCol[icl]->getCluster()){
      printf("    Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Address \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy() );
    }
  }
}
*/

cout<<"Readin empty half cluster size "<<p_emptyHalfClusterU.size()<<", "<<p_emptyHalfClusterV.size()<<endl;
/*
tmp_totE_U = 0;
for(int ih=0; ih<p_emptyHalfClusterU.size(); ih++) tmp_totE_U += p_emptyHalfClusterU[ih]->getEnergy();
tmp_totE_V = 0;
for(int ih=0; ih<p_emptyHalfClusterV.size(); ih++) tmp_totE_V += p_emptyHalfClusterV[ih]->getEnergy();
cout<<"Empty half cluster energy "<<tmp_totE_U<<", "<<tmp_totE_V<<endl;
*/

  //Start in U direction. 
  m_newClusUCol.clear();
  for(int ih=0; ih<p_HalfClusterU.size(); ih++){

    //Get all axis in U cluster.
    m_axisUCol.clear();
    if( settings.map_stringPars["ReadinAxisName"] == "AllAxis" ) m_axisUCol = p_HalfClusterU[ih]->getAllHalfClusterCol(); 
    else m_axisUCol = p_HalfClusterU[ih]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"]);
   
    //  Loop for 1DClusters.
    m_1dShowerUCol.clear();
    std::vector<const Cyber::Calo1DCluster*> m_1dclusCol = p_HalfClusterU[ih]->getCluster();
    for(int icl=0; icl<m_1dclusCol.size(); icl++){
      std::vector<const Cyber::CaloUnit*> m_bars = m_1dclusCol[icl]->getBars();

      //Find the seed with axis in 1DCluster: 
      for(auto iaxis: m_axisUCol){
        for(auto iseed: iaxis->getCluster() ){

          if(iseed->getBars().size()!=1) {
            std::cout<<"WARNING: Axis has more than one bars! Check!"<<'\t'<<iseed->getBars().size()<<std::endl;
            continue;
          }

          bool findSeed = false;
          if(find( m_bars.begin(), m_bars.end(), iseed->getBars()[0]) != m_bars.end()) findSeed = true;
          for(int ib=0; ib<m_bars.size(); ib++){
            if(findSeed) break;
            if(m_bars[ib]->getcellID()==iseed->getBars()[0]->getcellID()) {findSeed=true; break; }
          }
          if(findSeed) const_cast<Cyber::Calo1DCluster*>(m_1dclusCol[icl])->addSeed( iseed->getBars()[0] );

        }
      }
      //if(m_1dclusCol[icl]->getNseeds()==0) const_cast<Cyber::Calo1DCluster*>(m_1dclusCol[icl])->setSeed();   

      //Split cluster to showers
      std::vector<std::shared_ptr<Cyber::Calo1DCluster>> tmp_showers; tmp_showers.clear();
      ClusterSplitting( m_1dclusCol[icl], tmp_showers );
      //if(tmp_showers.size()==0) continue;
      m_1dShowerUCol.insert( m_1dShowerUCol.end(), tmp_showers.begin(), tmp_showers.end() );
      tmp_showers.clear();
    }

//double tmp_En = 0.;
//for(int i1d=0; i1d<m_1dShowerUCol.size(); i1d++) tmp_En += m_1dShowerUCol[i1d]->getEnergy();
//cout<<"  HalfClusterU #"<<ih<<": total En after energy splitting "<<tmp_En<<endl;

//cout<<"  After energy splitting: 1D cluster size "<<m_1dShowerUCol.size()<<", print check"<<endl;
//for(int i1d=0; i1d<m_1dShowerUCol.size(); i1d++){
//  auto ish = m_1dShowerUCol[i1d];
//  printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, MC map size %d \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(),  ish->getNseeds(), ish->getLinkedMCP().size() );
//  //printf(", cover tower size: %d: ", ish->getTowerID().size());
//  //for(int atw=0; atw<ish->getTowerID().size(); atw++) printf("[%d, %d, %d], ", ish->getTowerID()[atw][0], ish->getTowerID()[atw][1], ish->getTowerID()[atw][2] );
//  //cout<<endl;
//}
//cout<<endl;

    //Clean showers without seed.
    for(int ic=0; ic<m_1dShowerUCol.size(); ic++){
      if(m_1dShowerUCol[ic]->getNseeds()==0){
        if(MergeToClosestCluster( m_1dShowerUCol[ic].get(), m_1dShowerUCol )==StatusCode::SUCCESS){
          m_1dShowerUCol.erase(m_1dShowerUCol.begin()+ic);
          ic--;
        }
        else
          m_1dShowerUCol[ic]->setSeed();
      }
    }
    for(int i1d=0; i1d<m_1dShowerUCol.size(); i1d++) m_1dShowerUCol[i1d]->getLinkedMCPfromUnit();
    m_datacol.map_1DCluster["bk1DCluster"].insert( m_datacol.map_1DCluster["bk1DCluster"].end(), m_1dShowerUCol.begin(), m_1dShowerUCol.end() );

//tmp_En = 0.;
//for(int i1d=0; i1d<m_1dShowerUCol.size(); i1d++) tmp_En += m_1dShowerUCol[i1d]->getEnergy();
//cout<<"  HalfClusterU #"<<ih<<": total En after shower cleaning "<<tmp_En<<endl;
//cout<<"  After shower cleaning: 1D cluster size "<<m_1dShowerUCol.size()<<", print check"<<endl;
//for(int i1d=0; i1d<m_1dShowerUCol.size(); i1d++){
//  auto ish = m_1dShowerUCol[i1d];
//  printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, MC map size %d \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(),  ish->getNseeds(), ish->getLinkedMCP().size() );
//  //printf(", cover tower size: %d: ", ish->getTowerID().size());
//  //for(int atw=0; atw<ish->getTowerID().size(); atw++) printf("[%d, %d, %d], ", ish->getTowerID()[atw][0], ish->getTowerID()[atw][1], ish->getTowerID()[atw][2] );
//  //cout<<endl;
//}
//cout<<endl;

    //  Longitudinal linking: update clusters' energy.
    std::vector<std::shared_ptr<Cyber::CaloHalfCluster>> tmp_newClus; tmp_newClus.clear();
    LongitudinalLinking(m_1dShowerUCol, m_axisUCol, tmp_newClus);
    m_newClusUCol.insert(m_newClusUCol.end(), tmp_newClus.begin(), tmp_newClus.end());

//tmp_En = 0.;
//for(int acl=0; acl<tmp_newClus.size(); acl++) tmp_En += tmp_newClus[acl]->getEnergy();
//cout<<"  HalfClusterU #"<<ih<<": total En after LongiLinking "<<tmp_En<<endl;
//cout<<"  After Longitudinally linking: HalfCluster U size "<<tmp_newClus.size()<<", print check"<<endl;
//for(int icl=0; icl<tmp_newClus.size(); icl++){
//  cout<<"      In HFClusU #"<<icl<<": shower size = "<<tmp_newClus[icl]->getCluster().size()<<endl;
//  for(auto ish : tmp_newClus[icl]->getCluster()){
//    printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, MC map size %d \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(),  ish->getNseeds(), ish->getLinkedMCP().size() );
//    //printf(", cover tower size: %d: ", ish->getTowerID().size());
//    //for(int atw=0; atw<ish->getTowerID().size(); atw++) printf("[%d, %d, %d], ", ish->getTowerID()[atw][0], ish->getTowerID()[atw][1], ish->getTowerID()[atw][2] );
//    //cout<<endl;
//  }
//}
//cout<<endl;

    tmp_newClus.clear();
  }

  //Start in V direction
  m_newClusVCol.clear();
  for(int ih=0; ih<p_HalfClusterV.size(); ih++){
    m_axisVCol.clear();
//cout<<"Readin axis in V: "<<settings.map_stringPars["ReadinAxisName"]<<endl;
    if( settings.map_stringPars["ReadinAxisName"] == "AllAxis" ) m_axisVCol = p_HalfClusterV[ih]->getAllHalfClusterCol(); 
    else m_axisVCol = p_HalfClusterV[ih]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"]);

    m_1dShowerVCol.clear();
    p_HalfClusterV[ih]->sortBarShowersByLayer();
    std::vector<const Cyber::Calo1DCluster*> m_1dclusCol = p_HalfClusterV[ih]->getCluster();

//cout<<"1D Cluster size: "<<m_1dclusCol.size()<<", axis size: "<<m_axisVCol.size()<<endl;
    for(int icl=0; icl<m_1dclusCol.size(); icl++){
      std::vector<const Cyber::CaloUnit*> m_bars = m_1dclusCol[icl]->getBars();
//cout<<"  In 1DCluster #"<<icl<<": layer = "<<m_1dclusCol[icl]->getDlayer();   
    
      //Find the seed with axis in 1DCluster:
      for(auto iaxis: m_axisVCol){
        for(auto iseed: iaxis->getCluster() ){
          if(iseed->getBars().size()!=1) {
            std::cout<<"WARNING: Axis has more than one bars! Check!"<<'\t'<<iseed->getBars().size()<<std::endl;
            continue;
          }
          if(find( m_bars.begin(), m_bars.end(), iseed->getBars()[0]) != m_bars.end())
            const_cast<Cyber::Calo1DCluster*>(m_1dclusCol[icl])->addSeed( iseed->getBars()[0] );
        }
      }
//cout<<", matched seed size "<<m_1dclusCol[icl]->getNseeds()<<endl;
      //if(m_1dclusCol[icl]->getNseeds()==0) const_cast<Cyber::Calo1DCluster*>(m_1dclusCol[icl])->setSeed();   

      //Split cluster to showers
      std::vector<std::shared_ptr<Cyber::Calo1DCluster>> tmp_showers; tmp_showers.clear();
      ClusterSplitting( m_1dclusCol[icl], tmp_showers );
//cout<<"  In 1DCluster #"<<icl<<": splitted into "<<tmp_showers.size ()<<" showers: ";
//for(int is=0; is<tmp_showers.size(); is++){
//  printf(" (%.3f, %.3f, %.3f, %.3f) \t", tmp_showers[is].get()->getPos().x(), tmp_showers[is].get()->getPos().y(), tmp_showers[is].get()->getPos().z(), tmp_showers[is].get()->getEnergy());
//}
//cout<<endl;
//cout<<endl;
      //if(tmp_showers.size()==0) continue;
      m_1dShowerVCol.insert( m_1dShowerVCol.end(), tmp_showers.begin(), tmp_showers.end() );
      tmp_showers.clear();
    }

//double tmp_En = 0.;
//for(int i1d=0; i1d<m_1dShowerVCol.size(); i1d++) tmp_En += m_1dShowerVCol[i1d]->getEnergy();
//cout<<"  HalfClusterV #"<<ih<<": total En after energy splitting "<<tmp_En<<endl;

    //Clean showers without seed.
    for(int ic=0; ic<m_1dShowerVCol.size(); ic++){
      if(m_1dShowerVCol[ic]->getNseeds()==0){
        if(MergeToClosestCluster( m_1dShowerVCol[ic].get(), m_1dShowerVCol )==StatusCode::SUCCESS){
          m_1dShowerVCol.erase(m_1dShowerVCol.begin()+ic);
          ic--;
        }
        else
          m_1dShowerVCol[ic]->setSeed();
      }
    }

//cout<<"After splitting: print all 1D showers"<<endl;
//for(int is=0; is<m_1dShowerVCol.size(); is++){
//printf("  Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, address %p \n", m_1dShowerVCol[is].get()->getDlayer(), m_1dShowerVCol[is].get()->getPos().x(), m_1dShowerVCol[is].get()->getPos().y(), m_1dShowerVCol[is].get()->getPos().z(), m_1dShowerVCol[is].get()->getEnergy(), m_1dShowerVCol[is].get()->getBars().size(), m_1dShowerVCol[is].get()->getNseeds(), m_1dShowerVCol[is].get() );
//}
//cout<<endl;
//tmp_En = 0.;
//for(int i1d=0; i1d<m_1dShowerVCol.size(); i1d++) tmp_En += m_1dShowerVCol[i1d]->getEnergy();
//cout<<"  HalfClusterV #"<<ih<<": total En after shower cleaning "<<tmp_En<<endl;

    for(int i1d=0; i1d<m_1dShowerVCol.size(); i1d++) m_1dShowerVCol[i1d]->getLinkedMCPfromUnit();
    m_datacol.map_1DCluster["bk1DCluster"].insert( m_datacol.map_1DCluster["bk1DCluster"].end(), m_1dShowerVCol.begin(), m_1dShowerVCol.end() );


    //  Longitudinal linking: update clusters' energy.
    std::vector<std::shared_ptr<CaloHalfCluster>> tmp_newClus; tmp_newClus.clear();
    LongitudinalLinking(m_1dShowerVCol, m_axisVCol, tmp_newClus);
    m_newClusVCol.insert(m_newClusVCol.end(), tmp_newClus.begin(), tmp_newClus.end());
//double tmp_En = 0.;
//for(int acl=0; acl<tmp_newClus.size(); acl++) tmp_En += tmp_newClus[acl]->getEnergy();
//cout<<"  HalfClusterV #"<<ih<<": total En after LongiLinking "<<tmp_En<<endl;
    tmp_newClus.clear();
  }


  //Assign MCtruth info to the HalfCluster. 
  for(int ic=0; ic<m_newClusUCol.size(); ic++) m_newClusUCol[ic].get()->getLinkedMCPfromUnit();
  for(int ic=0; ic<m_newClusVCol.size(); ic++) m_newClusVCol[ic].get()->getLinkedMCPfromUnit();

  m_datacol.map_HalfCluster[settings.map_stringPars["OutputClusName"]+"U"] = m_newClusUCol; 
  m_datacol.map_HalfCluster[settings.map_stringPars["OutputClusName"]+"V"] = m_newClusVCol; 


cout<<"  After splitting: HalfCluster U size "<<m_newClusUCol.size()<<", print check"<<endl;
/*
for(int icl=0; icl<m_newClusUCol.size(); icl++){
  cout<<"      In HFClusU #"<<icl<<": shower size = "<<m_newClusUCol[icl]->getCluster().size()<<endl;
  for(auto ish : m_newClusUCol[icl]->getCluster()){
    printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, MC map size %d \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(),  ish->getNseeds(), ish->getLinkedMCP().size() );
    //printf(", cover tower size: %d: ", ish->getTowerID().size());
    //for(int atw=0; atw<ish->getTowerID().size(); atw++) printf("[%d, %d, %d], ", ish->getTowerID()[atw][0], ish->getTowerID()[atw][1], ish->getTowerID()[atw][2] );
    //cout<<endl;
  }
}
cout<<endl;
*/
cout<<"  After splitting: HalfCluster V size "<<m_newClusVCol.size()<<", print check"<<endl;
/*
for(int icl=0; icl<m_newClusVCol.size(); icl++){
  cout<<"      In HFClusV #"<<icl<<": shower size = "<<m_newClusVCol[icl]->getCluster().size()<<endl;
  for(auto ish : m_newClusVCol[icl]->getCluster()){
    printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, MC map size %d \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(),  ish->getNseeds(), ish->getLinkedMCP().size() );
    //printf(", cover tower size: %d: ", ish->getTowerID().size());
    //for(int atw=0; atw<ish->getTowerID().size(); atw++) printf("[%d, %d, %d], ", ish->getTowerID()[atw][0], ish->getTowerID()[atw][1], ish->getTowerID()[atw][2] );
    //cout<<endl;
  }
}
*/
/*
float totE_U = 0.;
float totE_V = 0.;
for(auto iter : m_newClusUCol) totE_U += iter->getEnergy();
for(auto iter : m_newClusVCol) totE_V += iter->getEnergy();
printf("    Before split to tower: HalfCluster size: (%d, %d), energy (%.3f, %.3f) \n", m_newClusUCol.size(), m_newClusVCol.size(), totE_U, totE_V );
cout<<"    Loop print HalfClusterU: "<<endl;
for(int icl=0; icl<m_newClusUCol.size(); icl++){
  cout<<"      In HFClusU #"<<icl<<": shower size = "<<m_newClusUCol[icl]->getCluster().size()<<", En = "<<m_newClusUCol[icl]->getEnergy()<<", type "<<m_newClusUCol[icl]->getType();
  printf(", Position (%.3f, %.3f, %.3f), address ",m_newClusUCol[icl]->getPos().x(), m_newClusUCol[icl]->getPos().y(), m_newClusUCol[icl]->getPos().z());
  printf(", cousin size %d, address: ", m_newClusUCol[icl]->getHalfClusterCol("CousinCluster").size());
  //for(int ics=0; ics<m_newClusUCol[icl]->getHalfClusterCol("CousinCluster").size(); ics++) printf("%p, ", m_newClusUCol[icl]->getHalfClusterCol("CousinCluster")[ics]);
  printf(", track size %d, address: ", m_newClusUCol[icl]->getAssociatedTracks().size());
  //for(int itrk=0; itrk<m_newClusUCol[icl]->getAssociatedTracks().size(); itrk++) printf("%p, ", m_newClusUCol[icl]->getAssociatedTracks()[itrk]);
  cout<<endl;
  //for(auto ish : m_newClusUCol[icl]->getCluster()){
  //  printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, seedID [%d, %d, %d, %d], Address %p \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(),  ish->getNseeds(), ish->getSeeds()[0]->getModule(), ish->getSeeds()[0]->getPart(), ish->getSeeds()[0]->getStave(), ish->getSeeds()[0]->getBar(), ish );
  //}
}
cout<<endl;

cout<<"    Loop print HalfClusterV: "<<endl;
for(int icl=0; icl<m_newClusVCol.size(); icl++){
  cout<<"      In HFClusV #"<<icl<<": shower size = "<<m_newClusVCol[icl]->getCluster().size()<<", En = "<<m_newClusVCol[icl]->getEnergy()<<", type "<<m_newClusVCol[icl]->getType();
  printf(", Position (%.3f, %.3f, %.3f), address ",m_newClusVCol[icl]->getPos().x(), m_newClusVCol[icl]->getPos().y(), m_newClusVCol[icl]->getPos().z());
  printf(", cousin size %d, address: ", m_newClusVCol[icl]->getHalfClusterCol("CousinCluster").size());
  //for(int ics=0; ics<m_newClusVCol[icl]->getHalfClusterCol("CousinCluster").size(); ics++) printf("%p, ", m_newClusVCol[icl]->getHalfClusterCol("CousinCluster")[ics]);
  printf(", track size %d, address: ", m_newClusVCol[icl]->getAssociatedTracks().size());
  //for(int itrk=0; itrk<m_newClusVCol[icl]->getAssociatedTracks().size(); itrk++) printf("%p, ", m_newClusVCol[icl]->getAssociatedTracks()[itrk]);
  cout<<endl;
  //for(auto ish : m_newClusVCol[icl]->getCluster()){
  //  printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, seedID [%d, %d, %d, %d], Address %p \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(), ish->getNseeds(), ish->getSeeds()[0]->getModule(), ish->getSeeds()[0]->getPart(), ish->getSeeds()[0]->getStave(), ish->getSeeds()[0]->getBar(), ish );
  //}
}
cout<<endl;
*/


  //Make tower 
  m_towerCol.clear(); 
  std::vector<Cyber::CaloHalfCluster*> tmp_newClusUCol, tmp_newClusVCol;
  tmp_newClusUCol.clear(); tmp_newClusVCol.clear();
  for(int icl=0; icl<m_newClusUCol.size(); icl++) tmp_newClusUCol.push_back( m_newClusUCol[icl].get() );
  for(int icl=0; icl<m_newClusVCol.size(); icl++) tmp_newClusVCol.push_back( m_newClusVCol[icl].get() );
  HalfClusterToTowers( tmp_newClusUCol, tmp_newClusVCol, p_emptyHalfClusterU, p_emptyHalfClusterV, m_towerCol );

  m_datacol.map_CaloCluster[settings.map_stringPars["OutputTowerName"]] = m_towerCol;

  m_datacol.map_BarCol["bkBar"].insert( m_datacol.map_BarCol["bkBar"].end(), m_bkCol.map_BarCol["bkBar"].begin(), m_bkCol.map_BarCol["bkBar"].end() );
  m_datacol.map_1DCluster["bk1DCluster"].insert( m_datacol.map_1DCluster["bk1DCluster"].end(), m_bkCol.map_1DCluster["bk1DCluster"].begin(), m_bkCol.map_1DCluster["bk1DCluster"].end() );
  m_datacol.map_2DCluster["bk2DCluster"].insert( m_datacol.map_2DCluster["bk2DCluster"].end(), m_bkCol.map_2DCluster["bk2DCluster"].begin(), m_bkCol.map_2DCluster["bk2DCluster"].end() );
  m_datacol.map_HalfCluster["bkHalfCluster"].insert( m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_bkCol.map_HalfCluster["bkHalfCluster"].begin(), m_bkCol.map_HalfCluster["bkHalfCluster"].end() );

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::ClearAlgorithm(){
  p_HalfClusterU.clear();
  p_HalfClusterV.clear();
  p_emptyHalfClusterU.clear();
  p_emptyHalfClusterV.clear();

  m_axisUCol.clear();
  m_axisVCol.clear();

  m_newClusUCol.clear();
  m_newClusVCol.clear();  
  m_1dShowerUCol.clear();
  m_1dShowerVCol.clear();
  m_towerCol.clear(); 
  m_bkCol.Clear();

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::LongitudinalLinking( std::vector<std::shared_ptr<Cyber::Calo1DCluster>>& m_showers,
                                                    std::vector<const Cyber::CaloHalfCluster*>& m_oldClusCol,
                                                    std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& m_newClusCol )

{
  if(m_showers.size()==0 || m_oldClusCol.size()==0) return StatusCode::SUCCESS;
  m_newClusCol.clear();

//double tmp_totE = 0;
//for(int i=0; i<m_showers.size(); i++) tmp_totE += m_showers[i]->getEnergy();
//cout<<"  In LongitudinalLinking: input 1DShower total energy : "<<tmp_totE<<". old axis size "<<m_oldClusCol.size()<<endl;
//cout<<"  Print all 1DShower (address)"<<endl;
//for(int i=0; i<m_showers.size(); i++){
//  printf("    Shower #%d, pos+E [%.3f, %.3f, %.3f, %.3f], Nbars %d, address \n", i, m_showers[i]->getPos().x(), m_showers[i]->getPos().y(), m_showers[i]->getPos().z(), m_showers[i]->getBars().size(), m_showers[i]->getEnergy() );
//}

  //Update old Hough clusters
  for(int ic=0; ic<m_oldClusCol.size(); ic++){
    std::shared_ptr<Cyber::CaloHalfCluster> m_newClus = std::make_shared<Cyber::CaloHalfCluster>();

    for(int is=0; is<m_oldClusCol[ic]->getCluster().size(); is++){
      const Cyber::Calo1DCluster* m_shower = m_oldClusCol[ic]->getCluster()[is];

      const Cyber::Calo1DCluster* m_selshower = NULL;
      bool fl_foundshower = false; 
      for(int js=0; js<m_showers.size(); js++){
        bool fl_inTower = false;
        for(int itw=0; itw<m_showers[js].get()->getTowerID().size() && !fl_inTower; itw++){
        for(int jtw=0; jtw<m_shower->getTowerID().size(); jtw++){
          if(m_showers[js].get()->getTowerID()[itw]==m_shower->getTowerID()[jtw]) {fl_inTower=true; break;}
        }}

        if( fl_inTower &&
            m_showers[js].get()->getDlayer() == m_shower->getDlayer() && 
            m_showers[js].get()->getSeeds().size()==1 && 
            (m_showers[js].get()->getSeeds()[0]->getPosition()-m_shower->getPos()).Mag()<Cyber::CaloUnit::barsize ) 
        {m_selshower = m_showers[js].get(); fl_foundshower=true; break; }
      }
      if(fl_foundshower && m_selshower!=NULL) m_newClus->addUnit( m_selshower );
    }

    for(int itrk=0; itrk<m_oldClusCol[ic]->getAssociatedTracks().size(); itrk++) m_newClus->addAssociatedTrack( m_oldClusCol[ic]->getAssociatedTracks()[itrk] );
    m_newClus->setType( m_oldClusCol[ic]->getType() );
    m_newClus->setHoughPars(m_oldClusCol[ic]->getHoughAlpha(), m_oldClusCol[ic]->getHoughRho() );
    m_newClus->setIntercept(m_oldClusCol[ic]->getHoughIntercept());
    m_newClus->fitAxis("");
    m_newClus->addHalfCluster(settings.map_stringPars["ReadinAxisName"], m_oldClusCol[ic]);
    m_newClusCol.push_back( m_newClus );
  }

//tmp_totE = 0.;
//for(int i=0; i<m_newClusCol.size(); i++) tmp_totE += m_newClusCol[i]->getEnergy();
//cout<<"  In Longi-Linking: raw new HFCluster size : "<<m_newClusCol.size()<<", totE "<<tmp_totE<<endl;

//cout<<"  Print 1DShower in new HFCluster"<<endl;
//for(int ih=0; ih<m_newClusCol.size(); ih++){
//  cout<<"    In new HFCluster #"<<ih<<endl;
//  for(int i=0; i<m_newClusCol[ih]->getCluster().size(); i++){
//    printf("      Shower #%d, pos+E [%.3f, %.3f, %.3f, %.3f], address \n", i, m_newClusCol[ih]->getCluster()[i]->getPos().x(), m_newClusCol[ih]->getCluster()[i]->getPos().y(), m_newClusCol[ih]->getCluster()[i]->getPos().z(), m_newClusCol[ih]->getCluster()[i]->getEnergy() );
//  }
//}


  std::vector<const Cyber::Calo1DCluster*> m_leftshowers; m_leftshowers.clear(); 
  for(int ish=0; ish<m_showers.size(); ish++){
    bool fl_inclus = false; 
    for(int icl=0; icl<m_newClusCol.size(); icl++){
      std::vector<const Calo1DCluster*> p_showerCol = m_newClusCol[icl]->getCluster();
      if( find(p_showerCol.begin(), p_showerCol.end(), m_showers[ish].get())!=p_showerCol.end() ){ fl_inclus=true; break; }
    }

    if(!fl_inclus) m_leftshowers.push_back( m_showers[ish].get() );
  }

//tmp_totE = 0.;
//for(int i=0; i<m_leftshowers.size(); i++) tmp_totE += m_leftshowers[i]->getEnergy();
//cout<<"  In Longi-Linking: left 1D showers size: "<<m_leftshowers.size()<<", totE "<<tmp_totE<<endl;
//cout<<"  Print left 1DShower (address)"<<endl;
//for(int i=0; i<m_leftshowers.size(); i++){
//  printf("    Shower #%d, pos+E [%.3f, %.3f, %.3f, %.3f], address %p \n", i, m_leftshowers[i]->getPos().x(), m_leftshowers[i]->getPos().y(), m_leftshowers[i]->getPos().z(), m_leftshowers[i]->getEnergy(), m_leftshowers[i] );
//}


  //Merge showers into closest Half cluster
  std::sort( m_leftshowers.begin(), m_leftshowers.end(), compLayer );
  for(int is=0; is<m_leftshowers.size(); is++) 
    MergeToClosestCluster( m_leftshowers[is], m_newClusCol );

//tmp_totE = 0.;
//for(int i=0; i<m_newClusCol.size(); i++) tmp_totE += m_newClusCol[i]->getEnergy();
//cout<<"  In Longi-Linking: after merge left showers: HFCluster size "<<m_newClusCol.size()<<", totE "<<tmp_totE<<endl;

 //Split the double-used 1DClusters  
  SplitOverlapCluster(m_newClusCol);
  for(int is=0; is<m_newClusCol.size(); is++) m_newClusCol[is]->sortBarShowersByLayer();

  //Merge showers in the same layer
  if(settings.map_boolPars["CompactHFCluster"]){
    for(int icl=0; icl<m_newClusCol.size(); icl++) 
      m_newClusCol[icl].get()->mergeClusterInLayer();
  }

//cout<<"  Print 1DShower in new HFCluster"<<endl;
//for(int ih=0; ih<m_newClusCol.size(); ih++){
//  cout<<"    In new HFCluster #"<<ih<<endl;
//  for(int i=0; i<m_newClusCol[ih]->getCluster().size(); i++){
//    printf("      Shower #%d, pos+E [%.3f, %.3f, %.3f, %.3f], address \n", i, m_newClusCol[ih]->getCluster()[i]->getPos().x(), m_newClusCol[ih]->getCluster()[i]->getPos().y(), m_newClusCol[ih]->getCluster()[i]->getPos().z(), m_newClusCol[ih]->getCluster()[i]->getEnergy() );
//  }
//}

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::HalfClusterToTowers( std::vector<Cyber::CaloHalfCluster*>& m_halfClusU,
                                                    std::vector<Cyber::CaloHalfCluster*>& m_halfClusV,
                                                    std::vector<Cyber::CaloHalfCluster*>& m_emptyClusU,
                                                    std::vector<Cyber::CaloHalfCluster*>& m_emptyClusV,
                                                    std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_towers )
{
cout<<"  HalfClusterToTowers: Input halfcluster size "<<m_halfClusU.size()<<", "<<m_halfClusV.size()<<endl;
cout<<"    Input empty cluster size "<<m_emptyClusU.size()<<", "<<m_emptyClusV.size()<<endl;

  m_towers.clear(); 

  std::map<std::vector<int>, std::vector<const Cyber::Calo2DCluster*> > map_2DCluster;
  std::map<std::vector<int>, std::vector<Cyber::CaloHalfCluster*> > map_HalfClusterU; 
  std::map<std::vector<int>, std::vector<Cyber::CaloHalfCluster*> > map_HalfClusterV; 

  std::map<std::vector<int>, std::vector<Cyber::CaloHalfCluster*> > map_emptyHalfClusterU;
  std::map<std::vector<int>, std::vector<Cyber::CaloHalfCluster*> > map_emptyHalfClusterV;  

  //Assign the empty clusters into tower
  for(int il=0; il<m_emptyClusU.size(); il++){
cout<<"  Empty clusterU #"<<il<<": towerID size "<<m_emptyClusU[il]->getTowerID().size()<<endl;
    //map_emptyHalfClusterU[m_emptyClusU[il]->getTowerID()[0]].push_back(m_emptyClusU[il]);
    map_HalfClusterU[m_emptyClusU[il]->getTowerID()[0]].push_back(m_emptyClusU[il]);
  }
  for(int il=0; il<m_emptyClusV.size(); il++){
cout<<"  Empty clusterV #"<<il<<": towerID size "<<m_emptyClusV[il]->getTowerID().size()<<endl;
    //map_emptyHalfClusterV[m_emptyClusV[il]->getTowerID()[0]].push_back(m_emptyClusV[il]);
    map_HalfClusterV[m_emptyClusV[il]->getTowerID()[0]].push_back(m_emptyClusV[il]);
  }


  //Split CaloHalfClusterU
  for(int il=0; il<m_halfClusU.size(); il++){
    if(m_halfClusU[il]->getCluster().size()==0) {std::cout<<"WARNING: Have an empty CaloHalfCluster! Skip it! "<<std::endl; continue;}

    //HalfCluster does not cover tower: 
    if( m_halfClusU[il]->getTowerID().size()==1 && m_halfClusU[il]->getCluster().size()>=settings.map_intPars["th_Nhit"]){
      std::vector<int> cl_towerID =  m_halfClusU[il]->getTowerID()[0]; 
      if(settings.map_boolPars["CompactHFCluster"]) m_halfClusU[il]->mergeClusterInLayer();
      map_HalfClusterU[cl_towerID].push_back(m_halfClusU[il]);
      continue; 
    }

    //CaloHalfCluster covers towers: Loop check showers. 
    std::map<std::vector<int>, Cyber::CaloHalfCluster* > tmp_LongiClusMaps; tmp_LongiClusMaps.clear();
    for(int is=0; is<m_halfClusU[il]->getCluster().size(); is++){
      const Cyber::Calo1DCluster* p_shower = m_halfClusU[il]->getCluster()[is];

      std::map<std::vector<int>, Cyber::Calo1DCluster* > tmp_1DClusMaps; tmp_1DClusMaps.clear();
      for(int ib=0; ib<p_shower->getBars().size(); ib++){
        std::vector<int> towerID(4);
        towerID[0] = p_shower->getBars()[ib]->getSystem();
        towerID[1] = p_shower->getBars()[ib]->getModule();
        towerID[2] = p_shower->getBars()[ib]->getStave();
        towerID[3] = p_shower->getBars()[ib]->getPart();        

        if(tmp_1DClusMaps.find(towerID)!=tmp_1DClusMaps.end()){
          tmp_1DClusMaps[towerID]->addUnit(p_shower->getBars()[ib]);
          tmp_1DClusMaps[towerID]->setIDInfo();
        }
        else{
          std::shared_ptr<Cyber::Calo1DCluster> tmp_clus = std::make_shared<Cyber::Calo1DCluster>();
          tmp_clus->addUnit( p_shower->getBars()[ib] );
          tmp_clus->setIDInfo();
          tmp_1DClusMaps[towerID] = tmp_clus.get();
          m_bkCol.map_1DCluster["bk1DCluster"].push_back( tmp_clus );
        }
      }

      if(tmp_1DClusMaps.size()>1){
        for(auto &iter : tmp_1DClusMaps){
          for(auto &iter1 : tmp_1DClusMaps){
            if(iter!= iter1) iter.second->addCousinCluster(iter1.second);
          }
        }
      }

      for(auto &iter: tmp_1DClusMaps){
        std::vector<int> towerID = iter.first;
        iter.second->setSeed();
        iter.second->setIDInfo();
        iter.second->getLinkedMCPfromUnit();

        if( tmp_LongiClusMaps.find( towerID )!=tmp_LongiClusMaps.end() ){
          tmp_LongiClusMaps[towerID]->addUnit( iter.second );
          tmp_LongiClusMaps[towerID]->setTowerID( towerID );
        }
        else{
          std::shared_ptr<Cyber::CaloHalfCluster> tmp_clus = std::make_shared<Cyber::CaloHalfCluster>();
          tmp_clus->addUnit( iter.second );
          tmp_clus->setTowerID( towerID );
          tmp_LongiClusMaps[towerID] = tmp_clus.get();
          m_bkCol.map_HalfCluster["bkHalfCluster"].push_back( tmp_clus );
        }        
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
      iter.second->setType(m_halfClusU[il]->getType());
      if(m_halfClusU[il]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"]).size()>0 ) 
        iter.second->addHalfCluster(settings.map_stringPars["ReadinAxisName"], m_halfClusU[il]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"])[0]);

      if(iter.second->getEnergy()<settings.map_floatPars["Eth_HFClus"]) continue;
      iter.second->addHalfCluster("ParentCluster", m_halfClusU[il]);
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
      map_HalfClusterV[cl_towerID].push_back(m_halfClusV[il]);
      continue;
    }

    //CaloHalfCluster covers towers: Loop check showers.
    std::map<std::vector<int>, Cyber::CaloHalfCluster* > tmp_LongiClusMaps; tmp_LongiClusMaps.clear();
    for(int is=0; is<m_halfClusV[il]->getCluster().size(); is++){
      const Cyber::Calo1DCluster* p_shower = m_halfClusV[il]->getCluster()[is];

      std::map<std::vector<int>, Cyber::Calo1DCluster* > tmp_1DClusMaps; tmp_1DClusMaps.clear();
      for(int ib=0; ib<p_shower->getBars().size(); ib++){
        std::vector<int> towerID(4);
        towerID[0] = p_shower->getBars()[ib]->getSystem();
        towerID[1] = p_shower->getBars()[ib]->getModule();
        towerID[2] = p_shower->getBars()[ib]->getStave();
        towerID[3] = p_shower->getBars()[ib]->getPart();

        if(tmp_1DClusMaps.find(towerID)!=tmp_1DClusMaps.end()){
          tmp_1DClusMaps[towerID]->addUnit(p_shower->getBars()[ib]);
          tmp_1DClusMaps[towerID]->setIDInfo();
        }
        else{
          std::shared_ptr<Cyber::Calo1DCluster> tmp_clus = std::make_shared<Cyber::Calo1DCluster>();
          tmp_clus->addUnit( p_shower->getBars()[ib] );
          tmp_clus->setIDInfo();
          tmp_1DClusMaps[towerID] = tmp_clus.get();
          m_bkCol.map_1DCluster["bk1DCluster"].push_back( tmp_clus );
        }
      }

      if(tmp_1DClusMaps.size()>1){
        for(auto &iter : tmp_1DClusMaps){
          for(auto &iter1 : tmp_1DClusMaps){
            if(iter!= iter1) iter.second->addCousinCluster(iter1.second);
          }
        }
      }

      for(auto &iter: tmp_1DClusMaps){
        std::vector<int> towerID = iter.first;
        iter.second->setSeed();
        iter.second->setIDInfo();
        iter.second->getLinkedMCPfromUnit();

        if( tmp_LongiClusMaps.find( towerID )!=tmp_LongiClusMaps.end() ){
          tmp_LongiClusMaps[towerID]->addUnit( iter.second );
          tmp_LongiClusMaps[towerID]->setTowerID( towerID );
        }
        else{
          std::shared_ptr<Cyber::CaloHalfCluster> tmp_clus = std::make_shared<Cyber::CaloHalfCluster>();
          tmp_clus->addUnit( iter.second );
          tmp_clus->setTowerID( towerID );
          tmp_LongiClusMaps[towerID] = tmp_clus.get();
          m_bkCol.map_HalfCluster["bkHalfCluster"].push_back( tmp_clus );
        }
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
      iter.second->setType(m_halfClusV[il]->getType());
      if(m_halfClusV[il]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"]).size()>0 )
        iter.second->addHalfCluster(settings.map_stringPars["ReadinAxisName"], m_halfClusV[il]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"])[0]);

      if(iter.second->getEnergy()<settings.map_floatPars["Eth_HFClus"]) continue;
      iter.second->addHalfCluster("ParentCluster", m_halfClusV[il]);
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

    std::vector<Cyber::CaloHalfCluster*> p_halfClusU = iterU.second; 
    std::vector<Cyber::CaloHalfCluster*> p_halfClusV = map_HalfClusterV[iterU.first]; 

    //Get ordered showers for looping in layers. 
    std::map<int, std::vector<const Cyber::Calo1DCluster*>> m_orderedShowerU; m_orderedShowerU.clear();
    std::map<int, std::vector<const Cyber::Calo1DCluster*>> m_orderedShowerV; m_orderedShowerV.clear();

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
    std::vector<const Cyber::Calo2DCluster*> m_blocks; m_blocks.clear();
    for(auto &iter1 : m_orderedShowerU){
      if( m_orderedShowerV.find( iter1.first )==m_orderedShowerV.end() ) continue; 
      std::shared_ptr<Cyber::Calo2DCluster> tmp_block = std::make_shared<Cyber::Calo2DCluster>();
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
printf("  In tower: [%d, %d, %d, %d] \n", m_towerID[0], m_towerID[1], m_towerID[2], m_towerID[3]);
    //Check cousin clusters: 
    std::vector<Cyber::CaloHalfCluster*> m_HFClusUInTower = map_HalfClusterU[m_towerID];
    for(auto &m_HFclus : m_HFClusUInTower){
      std::vector<const CaloHalfCluster*> tmp_delClus; tmp_delClus.clear();
      for(int ics=0; ics<m_HFclus->getHalfClusterCol("CousinCluster").size(); ics++){
        std::vector<int> tmp_towerID = m_HFclus->getHalfClusterCol("CousinCluster")[ics]->getTowerID()[0];

        if( map_2DCluster.find( tmp_towerID )==map_2DCluster.end() )
          tmp_delClus.push_back( m_HFclus->getHalfClusterCol("CousinCluster")[ics] );
      }
      for(int ics=0; ics<tmp_delClus.size(); ics++) m_HFclus->deleteCousinCluster( tmp_delClus[ics] );
    }

    std::vector<Cyber::CaloHalfCluster*> m_HFClusVInTower = map_HalfClusterV[m_towerID];
    for(auto &m_HFclus : m_HFClusVInTower){
      std::vector<const CaloHalfCluster*> tmp_delClus; tmp_delClus.clear();
      for(int ics=0; ics<m_HFclus->getHalfClusterCol("CousinCluster").size(); ics++){
        std::vector<int> tmp_towerID = m_HFclus->getHalfClusterCol("CousinCluster")[ics]->getTowerID()[0];

        if( map_2DCluster.find( tmp_towerID )==map_2DCluster.end() )
          tmp_delClus.push_back( m_HFclus->getHalfClusterCol("CousinCluster")[ics] );
      }
      for(int ics=0; ics<tmp_delClus.size(); ics++) m_HFclus->deleteCousinCluster( tmp_delClus[ics] );
    }

    //std::vector<const Cyber::CaloHalfCluster*> const_emptyClusU; const_emptyClusU.clear();
    //std::vector<const Cyber::CaloHalfCluster*> const_emptyClusV; const_emptyClusV.clear();
    //if(map_emptyHalfClusterU.find(m_towerID)!=map_emptyHalfClusterU.end()){
    //  for(int icl=0; icl<map_emptyHalfClusterU[m_towerID].size(); icl++){
    //    map_emptyHalfClusterU[m_towerID][icl]->getLinkedMCPfromUnit();
    //    const_emptyClusU.push_back(map_emptyHalfClusterU[m_towerID][icl]);
    //  }
    //  map_emptyHalfClusterU.erase(m_towerID);
    //}
    //if(map_emptyHalfClusterV.find(m_towerID)!=map_emptyHalfClusterV.end()){
    //  for(int icl=0; icl<map_emptyHalfClusterV[m_towerID].size(); icl++){
    //    map_emptyHalfClusterV[m_towerID][icl]->getLinkedMCPfromUnit();
    //    const_emptyClusV.push_back(map_emptyHalfClusterV[m_towerID][icl]);
    //  }
    //  map_emptyHalfClusterV.erase(m_towerID);
    //}
//cout<<"  Found empty half cluster size: "<<const_emptyClusU.size()<<", "<<const_emptyClusV.size()<<endl;

    //Convert to const
    std::vector<const Cyber::CaloHalfCluster*> const_HFClusU; const_HFClusU.clear();
    std::vector<const Cyber::CaloHalfCluster*> const_HFClusV; const_HFClusV.clear();
    for(int ics=0; ics<m_HFClusUInTower.size(); ics++){ m_HFClusUInTower[ics]->getLinkedMCPfromUnit(); const_HFClusU.push_back(m_HFClusUInTower[ics]); }
    for(int ics=0; ics<m_HFClusVInTower.size(); ics++){ m_HFClusVInTower[ics]->getLinkedMCPfromUnit(); const_HFClusV.push_back(m_HFClusVInTower[ics]); }

    std::shared_ptr<Cyber::Calo3DCluster> m_tower = std::make_shared<Cyber::Calo3DCluster>();
//printf("  Creating tower: [%d, %d, %d] \n", m_towerID[0], m_towerID[1], m_towerID[2]);
    m_tower->addTowerID( m_towerID );
    for(int i2d=0; i2d<map_2DCluster[m_towerID].size(); i2d++) m_tower->addUnit(map_2DCluster[m_towerID][i2d]);
    m_tower->setHalfClusters( settings.map_stringPars["OutputClusName"]+"U", const_HFClusU, 
                              settings.map_stringPars["OutputClusName"]+"V", const_HFClusV );
    //m_tower->setHalfClusters( "emptyHalfClusterU", const_emptyClusU,
    //                          "emptyHalfClusterV", const_emptyClusV );
    m_towers.push_back(m_tower);
  }

  //if(map_emptyHalfClusterU.size()>0 && map_emptyHalfClusterV.size()>0){
  //  for(auto iter: map_emptyHalfClusterU){
  //    std::vector<int> m_towerID = iter.first;

  //    if( map_emptyHalfClusterV.find(m_towerID)==map_emptyHalfClusterV.end() ){
  //      iter.second.clear();
  //      //map_emptyHalfClusterU.erase(m_towerID);
  //      continue;
  //    }

  //    std::vector<const Cyber::CaloHalfCluster*> const_emptyClusU; const_emptyClusU.clear();
  //    std::vector<const Cyber::CaloHalfCluster*> const_emptyClusV; const_emptyClusV.clear();
  //    for(int icl=0; icl<map_emptyHalfClusterU[m_towerID].size(); icl++){
  //      map_emptyHalfClusterU[m_towerID][icl]->getLinkedMCPfromUnit();
  //      const_emptyClusU.push_back(map_emptyHalfClusterU[m_towerID][icl]);
  //    }
  //    //map_emptyHalfClusterU.erase(m_towerID);
  //    for(int icl=0; icl<map_emptyHalfClusterV[m_towerID].size(); icl++){
  //      map_emptyHalfClusterV[m_towerID][icl]->getLinkedMCPfromUnit();
  //      const_emptyClusV.push_back(map_emptyHalfClusterV[m_towerID][icl]);
  //    }
  //    //map_emptyHalfClusterV.erase(m_towerID);

  //    std::shared_ptr<Cyber::Calo3DCluster> m_tower = std::make_shared<Cyber::Calo3DCluster>();
  //    m_tower->addTowerID( m_towerID );
  //    m_tower->setHalfClusters( "emptyHalfClusterU", const_emptyClusU,
  //                              "emptyHalfClusterV", const_emptyClusV );
  //    m_towers.push_back(m_tower);
  //  }
  //}

/*
cout<<"  After splitting: tower size "<<m_towers.size()<<". Print Tower: "<<endl;
for(auto it : m_towers){
  std::vector<const CaloHalfCluster*> m_HFClusUInTower = it->getHalfClusterUCol(settings.map_stringPars["OutputClusName"]+"U");
  std::vector<const CaloHalfCluster*> m_HFClusVInTower = it->getHalfClusterVCol(settings.map_stringPars["OutputClusName"]+"V");
  //std::vector<const CaloHalfCluster*> m_HFClusUInTower = it->getHalfClusterUCol("emptyHalfClusterU");
  //std::vector<const CaloHalfCluster*> m_HFClusVInTower = it->getHalfClusterVCol("emptyHalfClusterV");


cout<<"Check tower ID: ";
for(int i=0; i<it->getTowerID().size(); i++) printf("[%d, %d, %d], ", it->getTowerID()[i][0], it->getTowerID()[i][1], it->getTowerID()[i][2]);
cout<<endl;

  printf("    In Tower [%d, %d, %d], ", it->getTowerID()[0][0], it->getTowerID()[0][1], it->getTowerID()[0][2] );
  printf("    HalfCluster size: (%d, %d) \n", m_HFClusUInTower.size(), m_HFClusVInTower.size() );
  cout<<"    Loop print HalfClusterU: "<<endl;
  for(int icl=0; icl<m_HFClusUInTower.size(); icl++){
    cout<<"      In HFClusU #"<<icl<<": shower size = "<<m_HFClusUInTower[icl]->getCluster().size()<<", En = "<<m_HFClusUInTower[icl]->getEnergy()<<", type "<<m_HFClusUInTower[icl]->getType();
    printf(", Position (%.3f, %.3f, %.3f), ",m_HFClusUInTower[icl]->getPos().x(), m_HFClusUInTower[icl]->getPos().y(), m_HFClusUInTower[icl]->getPos().z());
    printf(", cousin size %d, address: ", m_HFClusUInTower[icl]->getHalfClusterCol("CousinCluster").size());
    //for(int ics=0; ics<m_HFClusUInTower[icl]->getHalfClusterCol("CousinCluster").size(); ics++) printf("%p, ", m_HFClusUInTower[icl]->getHalfClusterCol("CousinCluster")[ics]);
    printf(", track size %d, address: ", m_HFClusUInTower[icl]->getAssociatedTracks().size());
    //for(int itrk=0; itrk<m_HFClusUInTower[icl]->getAssociatedTracks().size(); itrk++) printf("%p, ", m_HFClusUInTower[icl]->getAssociatedTracks()[itrk]);
    cout<<endl;
    for(auto ish : m_HFClusUInTower[icl]->getCluster()){
      printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(),  ish->getNseeds() );
    }
  }
  cout<<endl;

  cout<<"    Loop print HalfClusterV: "<<endl;
  for(int icl=0; icl<m_HFClusVInTower.size(); icl++){
    cout<<"      In HFClusV #"<<icl<<": shower size = "<<m_HFClusVInTower[icl]->getCluster().size()<<", En = "<<m_HFClusVInTower[icl]->getEnergy()<<", type "<<m_HFClusVInTower[icl]->getType();
    printf(", Position (%.3f, %.3f, %.3f), ",m_HFClusVInTower[icl]->getPos().x(), m_HFClusVInTower[icl]->getPos().y(), m_HFClusVInTower[icl]->getPos().z());
    printf(", cousin size %d, address: ", m_HFClusVInTower[icl]->getHalfClusterCol("CousinCluster").size());
    //for(int ics=0; ics<m_HFClusVInTower[icl]->getHalfClusterCol("CousinCluster").size(); ics++) printf("%p, ", m_HFClusVInTower[icl]->getHalfClusterCol("CousinCluster")[ics]);
    printf(", track size %d, address: ", m_HFClusVInTower[icl]->getAssociatedTracks().size());
    //for(int itrk=0; itrk<m_HFClusVInTower[icl]->getAssociatedTracks().size(); itrk++) printf("%p, ", m_HFClusVInTower[icl]->getAssociatedTracks()[itrk]);
    cout<<endl;
    for(auto ish : m_HFClusVInTower[icl]->getCluster()){
      printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(), ish->getNseeds() );
    }
  }
  cout<<endl;
}
*/

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::ClusterSplitting(const Cyber::Calo1DCluster* m_cluster, std::vector<std::shared_ptr<Cyber::Calo1DCluster>>& outshCol ){

//cout<<"ClusterSplitting: input cluster seed size = "<<m_cluster->getSeeds().size()<<endl;
//cout<<"ClusterSplitting: input bar size = "<<m_cluster->getBars().size()<<endl;
//cout<<"Seed position and E: "<<endl;
//for(int a=0; a<m_cluster->getNseeds(); a++) printf("\t (%.3f, %.3f, %.3f, %.3f), barID %d \n", m_cluster->getSeeds()[a]->getPosition().x(), 
//                                                                                             m_cluster->getSeeds()[a]->getPosition().y(),
//                                                                                             m_cluster->getSeeds()[a]->getPosition().z(),
//                                                                                             m_cluster->getSeeds()[a]->getEnergy(), 
//                                                                                             m_cluster->getSeeds()[a]->getBar() );
//cout<<endl;
//double tmp_totE = 0.;
//for(int i=0; i<m_cluster->getBars().size(); i++) tmp_totE += m_cluster->getBars()[i]->getEnergy();
//cout<<"Total bar energy: "<<tmp_totE<<endl;

  //No seed in cluster: return origin cluster.
  if(m_cluster->getNseeds()==0) { 
    auto shower = m_cluster->Clone();
    outshCol.push_back(shower);
    //std::cout<<"WARNING: Still have no-seed cluster!!"<<std::endl; 
    return StatusCode::SUCCESS; 
  }

  //1 seed or second moment less than threshold: Not split. Turn cluster to shower and return
  else if(m_cluster->getNseeds()<2 || m_cluster->getScndMoment()<settings.map_floatPars["th_split"]){
    auto shower = m_cluster->Clone();
    outshCol.push_back(shower);
    return StatusCode::SUCCESS;
  }

  
  //Separated bars: 
  else if( m_cluster->getNseeds()>=m_cluster->getBars().size() ){
    for(int ish=0; ish<m_cluster->getNseeds(); ish++){
      std::shared_ptr<Cyber::Calo1DCluster> shower = std::make_shared<Cyber::Calo1DCluster>();
      shower->addUnit(m_cluster->getSeeds()[ish]);
      shower->addSeed(m_cluster->getSeeds()[ish]);
      shower->setIDInfo();

      outshCol.push_back(shower);
    }
    return StatusCode::SUCCESS;
  }
  
  //2 or more seeds, large second moment: Split
  int Nshower = m_cluster->getNseeds();
  int Nbars = m_cluster->getBars().size();
  double Eseed[Nshower] = {0};
  double weight[Nbars][Nshower] = {0};
  TVector3 SeedPos[Nshower];
  TVector3 SeedPos_Origin[Nshower];
  for(int is=0;is<Nshower;is++){ 
    SeedPos[is] = m_cluster->getSeeds()[is]->getPosition(); 
    SeedPos_Origin[is]=SeedPos[is]; 
    Eseed[is]=m_cluster->getSeeds()[is]->getEnergy(); 
  }
  //CalculateInitialEseed(m_cluster->getSeeds(), SeedPos, Eseed);
/*
  bool isConverge = false;
  bool isShifted = false; 
  int iter=0;
  TVector3 SeedPos_prev[Nshower];
  do{

    for(int ibar=0;ibar<m_cluster->getBars().size();ibar++){
      double Eexp[Nshower];
      double Eexp_tot=0;
      for(int is=0;is<Nshower;is++){ Eexp[is] = Eseed[is]*GetShowerProfile(m_cluster->getBars()[ibar]->getPosition(), SeedPos[is] ); Eexp_tot+= Eexp[is];}
      for(int is=0;is<Nshower;is++) weight[ibar][is] = Eexp[is]/Eexp_tot;
    }
    for(int is=0;is<Nshower;is++){
      SeedPos_prev[is]=SeedPos[is];

      double Etot=0;
      double Ebar[Nbars] = {0};
      double Emax = -99;
      for(int ib=0;ib<Nbars;ib++){
        Ebar[ib] = ( m_cluster->getBars()[ib]->getQ1()*weight[ib][is] + m_cluster->getBars()[ib]->getQ2()*weight[ib][is] )/2.; 
        Etot += Ebar[ib];
        if(Ebar[ib]>Emax) Emax = Ebar[ib];
      }

      TVector3 pos(0,0,0);
      for(int ib=0; ib<Nbars; ib++) pos += m_cluster->getBars()[ib]->getPosition() * (Ebar[ib]/Etot);
      SeedPos[is] = pos; 
      Eseed[is] = Emax;
    }
    isConverge=true;
    for(int is=0;is<Nshower;is++) if( (SeedPos_prev[is]-SeedPos[is]).Mag2()>2.89 ){ isConverge=false; break;}
 
    isShifted = false;
    for(int is=0;is<Nshower;is++) if( (SeedPos[is]-SeedPos_Origin[is]).Mag2()<15 ){ isShifted=true; break;}
    iter++;

  }
  while(iter<20 && !isConverge && !isShifted );
  if(iter>=20){
    std::cout<<"WARNING: Iteration time larger than 20! Might not converge!"<<std::endl;
    std::cout<<"  For Check: NBars: "<<m_cluster->getBars().size()<<"  Nseeds: "<<Nshower<<std::endl;
  }
*/


//Check: no iteration, just calculate the weight.
  for(int ibar=0;ibar<m_cluster->getBars().size();ibar++){
    double Eexp[Nshower];
    double Eexp_tot=0;
    for(int is=0;is<Nshower;is++){ Eexp[is] = Eseed[is]*GetShowerProfile(m_cluster->getBars()[ibar]->getPosition(), SeedPos[is] ); Eexp_tot+= Eexp[is];}
    for(int is=0;is<Nshower;is++){ 
      weight[ibar][is] = Eexp[is]/Eexp_tot;
      if(weight[ibar][is]<1e-5) weight[ibar][is] = 0.;
    }
  }



//cout<<"Print weight: "<<endl;
//for(int abar=0; abar<Nbars; abar++){
//for(int ash=0; ash<Nshower; ash++){
//  printf("%.3f \t",weight[abar][ash]);
//}
//cout<<endl;
//}

//cout<<"After weight calculation. "<<endl;
  //for(int is=0;is<Nshower;is++) SeedPos[is] = m_cluster->getSeeds()[is]->getPosition();
  for(int is=0;is<Nshower;is++){
//cout<<"  In shower #"<<is<<endl;
    std::vector<const Cyber::CaloUnit*> Bars; Bars.clear();
    int iseed=-1;
    int icount=0;
    double _Emax = -99;
    
    for(int ib=0;ib<Nbars;ib++){
      double barEn = m_cluster->getBars()[ib]->getEnergy() * weight[ib][is];
//cout<<"    bar#"<<ib<<" barID "<<m_cluster->getBars()[ib]->getBar()<<" En: "<<barEn<<endl;
      if(barEn<settings.map_floatPars["Eth_unit"]){ 
        continue;
      }
      auto bar = m_cluster->getBars()[ib]->Clone();
      bar->setQ( bar->getQ1()*weight[ib][is], bar->getQ2()*weight[ib][is]  );
      if( bar->getEnergy()>_Emax ) { _Emax=bar->getEnergy(); iseed=icount; }
      Bars.push_back(bar.get());
      icount++;
      m_bkCol.map_BarCol["bkBar"].push_back( bar );
//cout<<"      Found seed: iseed="<<iseed<<", icount="<<icount<<endl;
    }
    if(iseed<0) { std::cout<<"ERROR: Can not find seed(max energy bar) in this shower! Nbars = "<<Nbars<<", icount "<<icount<<", maxEn "<<_Emax<<". Please Check!"<<std::endl; continue; }
    //if( (Bars[iseed]->getPosition()-SeedPos[is]).Mag()>15 ) { std::cout<<"ERROR: MaxEnergy bar is too far with original seed! Please Check! iSeed = "<<iseed<<std::endl; }

    std::shared_ptr<Cyber::Calo1DCluster> shower = std::make_shared<Cyber::Calo1DCluster>();
    shower->setBars(Bars);
    shower->addSeed(Bars[iseed]);
    shower->setIDInfo(); 

    outshCol.push_back(shower);
  }
  
  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::SplitOverlapCluster( std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& m_HFClusCol ){

  if(m_HFClusCol.size()<2) return StatusCode::SUCCESS;

  for(int ic=0; ic<m_HFClusCol.size()-1; ic++){
    for(int jc=ic+1; jc<m_HFClusCol.size(); jc++){
      double E_ratio = m_HFClusCol[ic]->OverlapRatioE(m_HFClusCol[jc].get());
      if(E_ratio>0){

        std::vector<const Calo1DCluster*> m_ClusCol1 =  m_HFClusCol[ic]->getCluster();
        std::vector<const Calo1DCluster*> m_ClusCol2 =  m_HFClusCol[jc]->getCluster();
        std::vector<const Calo1DCluster*> m_overlapCl; m_overlapCl.clear();
        for(int i1d=0; i1d<m_ClusCol1.size(); i1d++){
          if(find(m_ClusCol2.begin(), m_ClusCol2.end(), m_ClusCol1[i1d])!=m_ClusCol2.end() ) m_overlapCl.push_back(m_ClusCol1[i1d]);
        }

        for(int io=0; io<m_overlapCl.size(); io++){

          //Split the cluster with +-1 layer cluster energy.
          int dlayer = m_overlapCl[io]->getDlayer();

          //  For  m_HFClusCol[ic]
          std::vector<const Calo1DCluster*> m_clus1_front = m_HFClusCol[ic]->getClusterInLayer(dlayer-1);
          std::vector<const Calo1DCluster*> m_clus1_back  = m_HFClusCol[ic]->getClusterInLayer(dlayer+1);
          std::vector<const Calo1DCluster*> m_clus2_front = m_HFClusCol[jc]->getClusterInLayer(dlayer-1);
          std::vector<const Calo1DCluster*> m_clus2_back  = m_HFClusCol[jc]->getClusterInLayer(dlayer+1);

          if( m_clus1_front.size()==0 && m_clus1_back.size()==0 && m_clus2_front.size()==0 && m_clus2_back.size()==0 ){
            //An isolated hit in both: delete in m_HFClusCol[ic]. TODO: maybe need to optimize based on profile?
            m_HFClusCol[ic]->deleteUnit(m_overlapCl[io]);
            continue;
          }
          else if( m_clus1_front.size()==0 && m_clus1_back.size()==0 ){
            //Only isolated in m_HFClusCol[ic]: delete in m_HFClusCol[ic].
            m_HFClusCol[ic]->deleteUnit(m_overlapCl[io]);
            continue;
          }
          else if( m_clus2_front.size()==0 && m_clus2_back.size()==0 ){
            m_HFClusCol[jc]->deleteUnit(m_overlapCl[io]);
            continue;
          }
          else{
            //Use the average energy
            double aveEn1 = 0.;
            for(int i1d=0; i1d<m_clus1_front.size(); i1d++) aveEn1 += m_clus1_front[i1d]->getEnergy();
            for(int i1d=0; i1d<m_clus1_back.size(); i1d++) aveEn1 += m_clus1_back[i1d]->getEnergy();
            aveEn1 = aveEn1/(m_clus1_front.size()+m_clus1_back.size());

            double aveEn2 = 0.;
            for(int i1d=0; i1d<m_clus2_front.size(); i1d++) aveEn2 += m_clus2_front[i1d]->getEnergy();
            for(int i1d=0; i1d<m_clus2_back.size(); i1d++) aveEn2 += m_clus2_back[i1d]->getEnergy();
            aveEn2 = aveEn2/(m_clus2_front.size()+m_clus2_back.size());

            std::shared_ptr<Cyber::Calo1DCluster> m_splitCl1 = std::make_shared<Cyber::Calo1DCluster>();
            std::shared_ptr<Cyber::Calo1DCluster> m_splitCl2 = std::make_shared<Cyber::Calo1DCluster>();

            for(int ib=0; ib<m_overlapCl[io]->getBars().size(); ib++){
              auto bar1 = m_overlapCl[io]->getBars()[ib]->Clone();
              bar1->setQ( bar1->getQ1()*aveEn1/(aveEn1+aveEn2), bar1->getQ2()*aveEn1/(aveEn1+aveEn2) );
              m_splitCl1->addUnit(bar1.get());

              auto bar2 = m_overlapCl[io]->getBars()[ib]->Clone();
              bar2->setQ( bar2->getQ1()*aveEn2/(aveEn1+aveEn2), bar2->getQ2()*aveEn2/(aveEn1+aveEn2) );
              m_splitCl2->addUnit(bar2.get());

              m_bkCol.map_BarCol["bkBar"].push_back(bar1);
              m_bkCol.map_BarCol["bkBar"].push_back(bar2);
            }
            m_splitCl1->setSeed();
            m_splitCl1->setIDInfo();
            m_splitCl1->getLinkedMCPfromUnit();
            m_splitCl2->setSeed();
            m_splitCl2->setIDInfo();
            m_splitCl2->getLinkedMCPfromUnit();

            m_HFClusCol[ic]->deleteUnit(m_overlapCl[io]);
            m_HFClusCol[ic]->addUnit(m_splitCl1.get());

            m_HFClusCol[jc]->deleteUnit(m_overlapCl[io]);
            m_HFClusCol[jc]->addUnit(m_splitCl2.get());

            m_bkCol.map_1DCluster["bk1DCluster"].push_back(m_splitCl1);
            m_bkCol.map_1DCluster["bk1DCluster"].push_back(m_splitCl2);
          }


        }//end loop m_HFClusCol
      }

    }
  }
  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::MergeToClosestCluster( Cyber::Calo1DCluster* iclus, std::vector<std::shared_ptr<Cyber::Calo1DCluster>>& clusvec ){

  int cLedge = iclus->getLeftEdge();
  int cRedge = iclus->getRightEdge();

//printf("ClusterMerging: input cluster system %d, layer %d, energy %.4f, Nbars %d, edge [%d, %d]\n",iclus->getSystem(), iclus->getDlayer(), iclus->getEnergy(), iclus->getBars().size(), cLedge, cRedge);

  //Find the closest cluster with iclus.
  int minD = 99;
  int index = -1;
  for(int icl=0; icl<clusvec.size(); icl++){
    if(clusvec[icl].get()->getNseeds()==0) continue;
    if( iclus->getTowerID() != clusvec[icl].get()->getTowerID() ) continue;
    if( iclus->getDlayer() != clusvec[icl].get()->getDlayer() ) continue;

    int iLedge = clusvec[icl].get()->getLeftEdge();
    int iRedge = clusvec[icl].get()->getRightEdge();

    //int dis = (cLedge-iRedge>0 ? cLedge-iRedge : iLedge-cRedge );
    int dis = min( abs(cLedge-iRedge), abs(iLedge-cRedge) );

//printf("  Loop in clusterCol #%d: Nbars %d, range [%d, %d], distance %d, energy ratio %.3f \n", icl, clusvec[icl].get()->getBars().size(), iLedge, iRedge, dis, iclus->getEnergy()/clusvec[icl]->getEnergy());
    if(dis>10) continue; //Don't merge to a too far cluster.
    if(dis<minD){ minD = dis; index=icl; }
  }
//cout<<"Selected closest cluster #"<<index<<", minD "<<minD<<endl;

  if(index<0) return StatusCode::FAILURE;

  //Merge to the selected cluster
  for(int icl=0; icl<iclus->getBars().size(); icl++)  clusvec[index].get()->addUnit(iclus->getBars()[icl]);


  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::MergeToClosestCluster( const Cyber::Calo1DCluster* m_shower, 
                                                      std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& m_clusters )
{
  if(m_clusters.size()==0) return StatusCode::SUCCESS;


  int minLayer = 99;
  int maxLayer = -99;
  for(int ic=0; ic<m_clusters.size(); ic++){
    if(minLayer>m_clusters[ic].get()->getBeginningDlayer())  minLayer = m_clusters[ic].get()->getBeginningDlayer();
    if(maxLayer<m_clusters[ic].get()->getEndDlayer())        maxLayer = m_clusters[ic].get()->getEndDlayer();
  }
  if(minLayer==99 || minLayer<0 || maxLayer<0) return StatusCode::SUCCESS;

  int dlayer = m_shower->getDlayer();
  TVector3 sh_pos = m_shower->getPos();

//cout<<"  Cluster range: ("<<minLayer<<", "<<maxLayer<<") "<<endl;
//cout<<"  Merging shower into cluster: Input shower ";
//printf(" (%.3f, %.3f, %.3f, %.3f), towerID[%d, %d, %d], layer #%d \n", sh_pos.X(), sh_pos.Y(), sh_pos.Z(), m_shower->getEnergy(), m_shower->getTowerID()[0][0], m_shower->getTowerID()[0][1], m_shower->getTowerID()[0][2], dlayer);

  double minR = 999;
  int index_cluster = -1;
  double m_distance = 0;
  for(int ic=0; ic<m_clusters.size(); ic++ ){
//printf("    Cluster #%d: pos (%.2f, %.2f, %.2f) \n", ic, m_clusters[ic]->getPos().x(), m_clusters[ic]->getPos().y(), m_clusters[ic]->getPos().z());
    if(dlayer<=minLayer)
      m_distance = (sh_pos-m_clusters[ic].get()->getClusterInLayer(m_clusters[ic].get()->getBeginningDlayer())[0]->getPos()).Mag();     
    else if(dlayer>=maxLayer)
      m_distance = (sh_pos-m_clusters[ic].get()->getClusterInLayer(m_clusters[ic].get()->getEndDlayer())[0]->getPos()).Mag();    
    else
      m_distance = (sh_pos-m_clusters[ic].get()->getPos()).Mag();

//cout<<"    Cluster #"<<ic<<": distance with shower = "<<m_distance<<endl;
    
    if( m_distance<minR )  { minR=m_distance; index_cluster=ic; }
  }

//printf("  minR = %.3f, in #cl %d, cluster size %d \n", minR, index_cluster, m_clusters.size());

  if(index_cluster>=0) m_clusters[index_cluster].get()->addUnit( m_shower );

  return StatusCode::SUCCESS;
}


void EnergySplittingAlg::CalculateInitialEseed( const std::vector<const Cyber::CaloUnit*>& Seeds, const TVector3* pos, double* Eseed){
//Calculate Eseed by solving a linear function:
// [ f(11) .. f(1mu)  .. ]   [E_seed 1 ]   [E_bar 1]
// [  ..   ..   ..    .. ] * [...      ] = [...    ]
// [ f(i1) .. f(imu)  .. ]   [E_seed mu]   [E_bar i]
// [ f(N1) ..   ..  f(NN)]   [E_seed N ]   [E_bar N]

  const int Nele = Seeds.size();
  std::vector<double> Eratio;
  std::vector<double> vec_Etot; //bar energy

  TVector vecE(Nele);
  TMatrix matrixR(Nele, Nele);

  for(int i=0;i<Nele;i++){ //Loop bar
    vecE[i] = Seeds[i]->getEnergy();
    for(int j=0;j<Nele;j++) matrixR[i][j] = GetShowerProfile(Seeds[i]->getPosition(), pos[j]);
  }


  matrixR.Invert();
  TVector sol = matrixR*vecE;

  for(int i=0;i<Nele;i++) Eseed[i] = sol[i];

//cout<<"Print initial Eseed: "<<endl;
//for(int i=0;i<Nele;i++) cout<<Eseed[i]<<'\t';
//cout<<endl;

}


double EnergySplittingAlg::GetShowerProfile(const TVector3& p_bar, const TVector3& p_seed ){
  TVector3 rpos = p_bar-p_seed;
  double dis = sqrt(rpos.Mag2());
  double a1, a2, b1, b2, Rm;
  a1=0.037; a2=0.265; b1=0.101; b2=0.437; Rm=1.868;

  return a1*exp(-b1*dis/Rm)+a2*exp(-b2*dis/Rm);
}
#endif
