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


StatusCode EnergySplittingAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  m_axisUCol.clear();
  m_axisVCol.clear();
  m_newClusUCol.clear();
  m_newClusVCol.clear();
  m_1dShowerUCol.clear();
  m_1dShowerVCol.clear();
  m_towerCol.clear(); 

  p_HalfClusterU.clear();
  p_HalfClusterV.clear(); 

  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColU"].size(); ih++)
    p_HalfClusterU.push_back( m_datacol.map_HalfCluster["HalfClusterColU"][ih].get() );
  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColV"].size(); ih++)
    p_HalfClusterV.push_back( m_datacol.map_HalfCluster["HalfClusterColV"][ih].get() );

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  //Input: HalfCluster (with 1D clusters and HoughAxis)
  //Output: Create Towers for the matching. 

  if( p_HalfClusterU.size() + p_HalfClusterV.size()==0 ) {
    std::cout<<"EnergySplittingAlg: No HalfCluster input"<<std::endl;
    return StatusCode::SUCCESS;
  }

  //Start in U direction. 
  m_newClusUCol.clear();
  for(int ih=0; ih<p_HalfClusterU.size(); ih++){

    //Get all axis in U cluster.
    m_axisUCol.clear();
    if( settings.map_stringPars["ReadinAxisName"] == "AllAxis" ) m_axisUCol = p_HalfClusterU[ih]->getAllHalfClusterCol(); 
    else m_axisUCol = p_HalfClusterU[ih]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"]);
   
    //  Loop for 1DClusters.
    m_1dShowerUCol.clear();
    std::vector<const PandoraPlus::Calo1DCluster*> m_1dclusCol = p_HalfClusterU[ih]->getCluster();
    for(int icl=0; icl<m_1dclusCol.size(); icl++){
      std::vector<const PandoraPlus::CaloUnit*> m_bars = m_1dclusCol[icl]->getBars();

      //Find the seed with axis in 1DCluster: 
      for(auto iaxis: m_axisUCol){
        for(auto iseed: iaxis->getCluster() ){
          if(iseed->getBars().size()!=1) {
            std::cout<<"WARNING: Axis has more than one bars! Check!"<<'\t'<<iseed->getBars().size()<<std::endl;
            continue;
          }
          if(find( m_bars.begin(), m_bars.end(), iseed->getBars()[0]) != m_bars.end()) 
            const_cast<PandoraPlus::Calo1DCluster*>(m_1dclusCol[icl])->addSeed( iseed->getBars()[0] );
        }
      }
      //if(m_1dclusCol[icl]->getNseeds()==0) const_cast<PandoraPlus::Calo1DCluster*>(m_1dclusCol[icl])->setSeed();   

      //Split cluster to showers
      std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>> tmp_showers; tmp_showers.clear();
      ClusterSplitting( m_1dclusCol[icl], tmp_showers, m_datacol.map_BarCol["bkBars"] );
      //if(tmp_showers.size()==0) continue;
      m_1dShowerUCol.insert( m_1dShowerUCol.end(), tmp_showers.begin(), tmp_showers.end() );
      tmp_showers.clear();
    }

    //Clean showers without seed.
    for(int ic=0; ic<m_1dShowerUCol.size(); ic++){
      if(m_1dShowerUCol[ic]->getNseeds()==0){
        MergeToClosestCluster( m_1dShowerUCol[ic].get(), m_1dShowerUCol );
        m_1dShowerUCol.erase(m_1dShowerUCol.begin()+ic);
        ic--;
      }
    }
    m_datacol.map_1DCluster["bk1DCluster"].insert( m_datacol.map_1DCluster["bk1DCluster"].end(), m_1dShowerUCol.begin(), m_1dShowerUCol.end() );

    //  Longitudinal linking: update clusters' energy.
    std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> tmp_newClus; tmp_newClus.clear();
    LongitudinalLinking(m_1dShowerUCol, m_axisUCol, tmp_newClus);
    m_newClusUCol.insert(m_newClusUCol.end(), tmp_newClus.begin(), tmp_newClus.end());
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
    std::vector<const PandoraPlus::Calo1DCluster*> m_1dclusCol = p_HalfClusterV[ih]->getCluster();

//cout<<"1D Cluster size: "<<m_1dclusCol.size()<<", axis size: "<<m_axisVCol.size()<<endl;
    for(int icl=0; icl<m_1dclusCol.size(); icl++){
      std::vector<const PandoraPlus::CaloUnit*> m_bars = m_1dclusCol[icl]->getBars();
//cout<<"  In 1DCluster #"<<icl<<": layer = "<<m_1dclusCol[icl]->getDlayer();   
    
      //Find the seed with axis in 1DCluster:
      for(auto iaxis: m_axisVCol){
        for(auto iseed: iaxis->getCluster() ){
          if(iseed->getBars().size()!=1) {
            std::cout<<"WARNING: Axis has more than one bars! Check!"<<'\t'<<iseed->getBars().size()<<std::endl;
            continue;
          }
          if(find( m_bars.begin(), m_bars.end(), iseed->getBars()[0]) != m_bars.end())
            const_cast<PandoraPlus::Calo1DCluster*>(m_1dclusCol[icl])->addSeed( iseed->getBars()[0] );
        }
      }
//cout<<", matched seed size "<<m_1dclusCol[icl]->getNseeds()<<endl;
      //if(m_1dclusCol[icl]->getNseeds()==0) const_cast<PandoraPlus::Calo1DCluster*>(m_1dclusCol[icl])->setSeed();   

      //Split cluster to showers
      std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>> tmp_showers; tmp_showers.clear();
      ClusterSplitting( m_1dclusCol[icl], tmp_showers, m_datacol.map_BarCol["bkBars"] );
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

    //Clean showers without seed.
    for(int ic=0; ic<m_1dShowerVCol.size(); ic++){
      if(m_1dShowerVCol[ic]->getNseeds()==0){
        MergeToClosestCluster( m_1dShowerVCol[ic].get(), m_1dShowerVCol );
        //delete m_1dShowerVCol[ic]; m_1dShowerVCol[ic]=NULL;
        m_1dShowerVCol.erase(m_1dShowerVCol.begin()+ic);
        ic--;
      }
    }

//cout<<"After splitting: print all 1D showers"<<endl;
//for(int is=0; is<m_1dShowerVCol.size(); is++){
//printf("  Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, address %p \n", m_1dShowerVCol[is].get()->getDlayer(), m_1dShowerVCol[is].get()->getPos().x(), m_1dShowerVCol[is].get()->getPos().y(), m_1dShowerVCol[is].get()->getPos().z(), m_1dShowerVCol[is].get()->getEnergy(), m_1dShowerVCol[is].get()->getBars().size(), m_1dShowerVCol[is].get()->getNseeds(), m_1dShowerVCol[is].get() );
//}
//cout<<endl;


    m_datacol.map_1DCluster["bk1DCluster"].insert( m_datacol.map_1DCluster["bk1DCluster"].end(), m_1dShowerVCol.begin(), m_1dShowerVCol.end() );


    //  Longitudinal linking: update clusters' energy.
    std::vector<std::shared_ptr<CaloHalfCluster>> tmp_newClus; tmp_newClus.clear();
    LongitudinalLinking(m_1dShowerVCol, m_axisVCol, tmp_newClus);
    m_newClusVCol.insert(m_newClusVCol.end(), tmp_newClus.begin(), tmp_newClus.end());
    tmp_newClus.clear();
  }


  //Assign MCtruth info to the HalfCluster. 
  for(int ic=0; ic<m_newClusUCol.size(); ic++) m_newClusUCol[ic].get()->getLinkedMCPfromUnit();
  for(int ic=0; ic<m_newClusVCol.size(); ic++) m_newClusVCol[ic].get()->getLinkedMCPfromUnit();

  m_datacol.map_HalfCluster[settings.map_stringPars["OutputClusName"]+"U"] = m_newClusUCol; 
  m_datacol.map_HalfCluster[settings.map_stringPars["OutputClusName"]+"V"] = m_newClusVCol; 



  //Make tower 
  m_towerCol.clear(); 
  std::vector<PandoraPlus::CaloHalfCluster*> tmp_newClusUCol, tmp_newClusVCol;
  tmp_newClusUCol.clear(); tmp_newClusVCol.clear();
  for(int icl=0; icl<m_newClusUCol.size(); icl++) tmp_newClusUCol.push_back( m_newClusUCol[icl].get() );
  for(int icl=0; icl<m_newClusVCol.size(); icl++) tmp_newClusVCol.push_back( m_newClusVCol[icl].get() );
  HalfClusterToTowers( tmp_newClusUCol, tmp_newClusVCol, m_towerCol, 
                       m_datacol.map_HalfCluster["bkHalfCluster"], 
                       m_datacol.map_1DCluster["bk1DCluster"], 
                       m_datacol.map_2DCluster["bk2DCluster"] );  

  m_datacol.map_CaloCluster[settings.map_stringPars["OutputTowerName"]] = m_towerCol;

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::ClearAlgorithm(){
  p_HalfClusterU.clear();
  p_HalfClusterV.clear();

  m_axisUCol.clear();
  m_axisVCol.clear();

  m_newClusUCol.clear();
  m_newClusVCol.clear();  
  m_1dShowerUCol.clear();
  m_1dShowerVCol.clear();
  m_towerCol.clear(); 

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::LongitudinalLinking( std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>>& m_showers,
                                                    std::vector<const PandoraPlus::CaloHalfCluster*>& m_oldClusCol,
                                                    std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& m_newClusCol )

{
  if(m_showers.size()==0 || m_oldClusCol.size()==0) return StatusCode::SUCCESS;
  m_newClusCol.clear();

  //Update old Hough clusters
  for(int ic=0; ic<m_oldClusCol.size(); ic++){
    std::shared_ptr<PandoraPlus::CaloHalfCluster> m_newClus = std::make_shared<PandoraPlus::CaloHalfCluster>();

    for(int is=0; is<m_oldClusCol[ic]->getCluster().size(); is++){
      const PandoraPlus::Calo1DCluster* m_shower = m_oldClusCol[ic]->getCluster()[is];

      const PandoraPlus::Calo1DCluster* m_selshower = NULL;
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
            (m_showers[js].get()->getSeeds()[0]->getPosition()-m_shower->getPos()).Mag()<10 ) 
        {m_selshower = m_showers[js].get(); fl_foundshower=true; break; }
      }
      if(fl_foundshower && m_selshower!=NULL) m_newClus->addUnit( m_selshower );
    }

    for(int itrk=0; itrk<m_oldClusCol[ic]->getAssociatedTracks().size(); itrk++) m_newClus->addAssociatedTrack( m_oldClusCol[ic]->getAssociatedTracks()[itrk] );
    m_newClus->setHoughPars(m_oldClusCol[ic]->getHoughAlpha(), m_oldClusCol[ic]->getHoughRho() );
    m_newClus->setIntercept(m_oldClusCol[ic]->getHoughIntercept());
    m_newClus->fitAxis("");
    m_newClusCol.push_back( m_newClus );
  }

/*
cout<<"  In Longi-Linking: raw new HFCluster is: "<<endl;
for(int icl=0; icl<m_newClusCol.size(); icl++){
  cout<<"    In HFClusU #"<<icl<<": shower size = "<<m_newClusCol[icl]->getCluster().size()<<endl;
  for(auto ish : m_newClusCol[icl]->getCluster()){
    printf("      Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, ", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(),  ish->getNseeds() );
    printf("Seed En %.3f, in tower [%d, %d, %d] \n",ish->getSeeds()[0]->getEnergy(), ish->getSeeds()[0]->getModule(), ish->getSeeds()[0]->getPart(), ish->getSeeds()[0]->getStave());
    printf("        Shower cover tower size: %d: ", ish->getTowerID().size());
    for(int atw=0; atw<ish->getTowerID().size(); atw++) printf("[%d, %d, %d], ", ish->getTowerID()[atw][0], ish->getTowerID()[atw][1], ish->getTowerID()[atw][2] );
    cout<<endl;
    cout<<endl;
  }
}
cout<<endl;
*/

  std::vector<const PandoraPlus::Calo1DCluster*> m_leftshowers; m_leftshowers.clear(); 
  for(int ish=0; ish<m_showers.size(); ish++){
    bool fl_inclus = false; 
    for(int icl=0; icl<m_newClusCol.size(); icl++){
      std::vector<const Calo1DCluster*> p_showerCol = m_newClusCol[icl]->getCluster();
      if( find(p_showerCol.begin(), p_showerCol.end(), m_showers[ish].get())!=p_showerCol.end() ){ fl_inclus=true; break; }
    }

    if(!fl_inclus) m_leftshowers.push_back( m_showers[ish].get() );
  }

/*
cout<<"  In Longi-Linking: left 1D showers size: "<<m_leftshowers.size()<<endl;
for(auto ish : m_leftshowers){
  printf("    Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, Address %p \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(), ish->getNseeds(), ish);
  for(int iseed=0; iseed<ish->getNseeds(); iseed++)
    printf("      Seed #%d: (%.3f, %.3f, %.3f, %.3f), towerID [%d, %d, %d], Address %p \n", iseed,
      ish->getSeeds()[iseed]->getPosition().x(),
      ish->getSeeds()[iseed]->getPosition().y(),
      ish->getSeeds()[iseed]->getPosition().z(),
      ish->getSeeds()[iseed]->getEnergy(),
      ish->getSeeds()[iseed]->getModule(),
      ish->getSeeds()[iseed]->getPart(),
      ish->getSeeds()[iseed]->getStave(),
      ish->getSeeds()[iseed] );
  cout<<endl;
}
*/
  //Merge showers into closest Half cluster
  std::sort( m_leftshowers.begin(), m_leftshowers.end(), compLayer );
  for(int is=0; is<m_leftshowers.size(); is++) 
    MergeToClosestCluster( m_leftshowers[is], m_newClusCol );

  //Merge showers in the same layer
  if(settings.map_boolPars["CompactHFCluster"]){
    for(int icl=0; icl<m_newClusCol.size(); icl++) 
      m_newClusCol[icl].get()->mergeClusterInLayer();
  }

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::HalfClusterToTowers( std::vector<PandoraPlus::CaloHalfCluster*>& m_halfClusU,
                                                    std::vector<PandoraPlus::CaloHalfCluster*>& m_halfClusV,
                                                    std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>& m_towers,
                                                    std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& bk_HFclus,
                                                    std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>>& bk_1Dclus,
                                                    std::vector<std::shared_ptr<PandoraPlus::Calo2DCluster>>& bk_2Dclus )
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
      map_HalfClusterU[cl_towerID].push_back(m_halfClusU[il]);
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
        bk_HFclus.push_back( tmp_clus );
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
        bk_HFclus.push_back( tmp_clus );
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
      bk_2Dclus.push_back( tmp_block );
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

/*
cout<<"  After splitting: tower size "<<m_towers.size()<<". Print Tower: "<<endl;
for(auto it : m_towers){
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

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::ClusterSplitting(const PandoraPlus::Calo1DCluster* m_cluster, std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>>& outshCol, std::vector<std::shared_ptr<PandoraPlus::CaloUnit>>& bk_bars ){

//cout<<"ClusterSplitting: input cluster seed size = "<<m_cluster->getSeeds().size()<<endl;
//cout<<"ClusterSplitting: input bar size = "<<m_cluster->getBars().size()<<endl;
//cout<<"Seed position and E: ";
//for(int a=0; a<m_cluster->getNseeds(); a++) printf(" (%.2f, %.2f, %.2f, %.2f) \t", m_cluster->getSeeds()[a]->getPosition().x(), 
//                                                                                   m_cluster->getSeeds()[a]->getPosition().y(),
//                                                                                   m_cluster->getSeeds()[a]->getPosition().z(),
//                                                                                   m_cluster->getSeeds()[a]->getEnergy() );
//cout<<endl;

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
      std::shared_ptr<PandoraPlus::Calo1DCluster> shower = std::make_shared<PandoraPlus::Calo1DCluster>();
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
  for(int is=0;is<Nshower;is++){ SeedPos[is] = m_cluster->getSeeds()[is]->getPosition(); SeedPos_Origin[is]=SeedPos[is]; Eseed[is]=m_cluster->getSeeds()[is]->getEnergy(); }
  //CalculateInitialEseed(m_cluster->getSeeds(), SeedPos, Eseed);

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
    std::vector<const PandoraPlus::CaloUnit*> Bars; Bars.clear();
    int iseed=-1;
    int icount=0;
    double _Emax = -99;
    
    for(int ib=0;ib<Nbars;ib++){
      double barEn = (m_cluster->getBars()[ib]->getQ1()*weight[ib][is] + m_cluster->getBars()[ib]->getQ2()*weight[ib][is])/2.;
//cout<<"    bar#"<<ib<<" En: "<<barEn;
      if(barEn<settings.map_floatPars["Eth_unit"]) continue;

      auto bar = m_cluster->getBars()[ib]->Clone();
      bar->setQ( bar->getQ1()*weight[ib][is], bar->getQ2()*weight[ib][is]  );
      if( bar->getEnergy()>_Emax ) { _Emax=bar->getEnergy(); iseed=icount; }
      Bars.push_back(bar.get());
      icount++;
      bk_bars.push_back( bar );
//cout<<", iseed="<<iseed<<", icount="<<icount<<endl;
    }
    if(iseed<0) { std::cout<<"ERROR: Can not find seed(max energy bar) in this shower! Nbars = "<<Nbars<<". Please Check!"<<std::endl; continue; }
    //if( (Bars[iseed]->getPosition()-SeedPos[is]).Mag()>15 ) { std::cout<<"ERROR: MaxEnergy bar is too far with original seed! Please Check! iSeed = "<<iseed<<std::endl; }

    std::shared_ptr<PandoraPlus::Calo1DCluster> shower = std::make_shared<PandoraPlus::Calo1DCluster>();
    shower->setBars(Bars);
    shower->addSeed(Bars[iseed]);
    shower->setIDInfo(); 


    outshCol.push_back(shower);
  }
  
  return StatusCode::SUCCESS;
}

StatusCode EnergySplittingAlg::MergeToClosestCluster( PandoraPlus::Calo1DCluster* iclus, std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>>& clusvec ){

  int cLedge = iclus->getLeftEdge();
  int cRedge = iclus->getRightEdge();

//printf("ClusterMerging: input cluster edge [%d, %d]\n", cLedge, cRedge);

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

//printf("  Loop in clusterCol #%d: range [%d, %d], distance %d, energy ratio %.3f \n", icl, iLedge, iRedge, dis, iclus->getEnergy()/clusvec[icl]->getEnergy());
    if(dis>10) continue; //Don't merge to a too far cluster.
    if(dis<minD){ minD = dis; index=icl; }
  }
//cout<<"Selected closest cluster #"<<index<<endl;

  if(index<0) return StatusCode::FAILURE;

  //Merge to the selected cluster
  for(int icl=0; icl<iclus->getBars().size(); icl++)  clusvec[index].get()->addUnit(iclus->getBars()[icl]);


  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::MergeToClosestCluster( const PandoraPlus::Calo1DCluster* m_shower, 
                                                      std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& m_clusters )
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
//printf(" (%.3f, %.3f, %.3f, %.3f), layer #%d \n", sh_pos.X(), sh_pos.Y(), sh_pos.Z(), m_shower->getEnergy(), dlayer);

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


void EnergySplittingAlg::CalculateInitialEseed( const std::vector<const PandoraPlus::CaloUnit*>& Seeds, const TVector3* pos, double* Eseed){
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
