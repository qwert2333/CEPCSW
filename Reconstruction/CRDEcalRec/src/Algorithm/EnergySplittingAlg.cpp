#ifndef _ENERGYSPLITTING_ALG_C
#define _ENERGYSPLITTING_ALG_C

#include "Algorithm/EnergySplittingAlg.h"

StatusCode EnergySplittingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  if(settings.map_floatPars.find("th_split")==settings.map_floatPars.end()) settings.map_floatPars["th_split"] = -1;
  if(settings.map_floatPars.find("Eth_Seed")==settings.map_floatPars.end()) settings.map_floatPars["Eth_Seed"] = 0.005;
  if(settings.map_floatPars.find("Eth_ShowerAbs")==settings.map_floatPars.end()) settings.map_floatPars["Eth_Shower"] = 0.005;
  if(settings.map_stringPars.find("OutputClusName")==settings.map_stringPars.end()) settings.map_stringPars["OutputClusName"] = "ESHalfCluster";
  if(settings.map_stringPars.find("OutputTowerName")==settings.map_stringPars.end()) settings.map_stringPars["OutputTowerName"] = "ESTower";

  return StatusCode::SUCCESS;
};


StatusCode EnergySplittingAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  m_axisUCol.clear();
  m_axisVCol.clear();
  m_newClusUCol.clear();
  m_newClusVCol.clear();
  m_1dclusCol.clear(); 
  m_1dShowerUCol.clear();
  m_1dShowerVCol.clear();
  m_towerCol.clear(); 

  p_HalfClustersU = &(m_datacol.map_HalfCluster["HalfClusterColU"]);
  p_HalfClustersV = &(m_datacol.map_HalfCluster["HalfClusterColV"]);


  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  //Input: HalfCluster (with 1D clusters and HoughAxis)
  //Output: Create Towers for the matching. 

  if( p_HalfClustersU->size() + p_HalfClustersV->size()==0 ) {
    std::cout<<"EnergySplittingAlg: No HalfCluster input"<<std::endl;
    return StatusCode::SUCCESS;
  }

  //Start in U direction. 
  m_newClusUCol.clear();
  for(int ih=0; ih<p_HalfClustersU->size(); ih++){
cout<<"    In HalfClusterU #"<<ih;
    //Get all axis in U cluster.
    m_axisUCol.clear();
    m_axisUCol = p_HalfClustersU->at(ih)->getAllHalfClusterCol(); 
cout<<", axis size "<<m_axisUCol.size();
   
    //  Loop for 1DClusters.
    m_1dShowerUCol.clear();
    m_1dclusCol.clear();
    m_1dclusCol = p_HalfClustersU->at(ih)->getCluster();

cout<<", 1DClusters: size = "<<m_1dclusCol.size();

cout<<"  Print HoughAxisUCol: "<<endl;
for(int icl=0; icl<m_axisUCol.size(); icl++){
cout<<"In HoughAxis #"<<icl<<endl;
for(int i=0; i<m_axisUCol[icl]->getCluster().size(); i++)
  cout<<'\t'<<m_axisUCol[icl]->getCluster()[i]->getDlayer()<<'\t'<<m_axisUCol[icl]->getCluster()[i]->getEnergy()<<endl;
}

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
      if(m_1dclusCol[icl]->getNseeds()==0) const_cast<PandoraPlus::Calo1DCluster*>(m_1dclusCol[icl])->setSeed();   

      //Split cluster to showers
      std::vector<const PandoraPlus::Calo1DCluster*> tmp_showers; tmp_showers.clear();
      ClusterSplitting( m_1dclusCol[icl], tmp_showers );
      //if(tmp_showers.size()==0) continue;
      m_1dShowerUCol.insert( m_1dShowerUCol.end(), tmp_showers.begin(), tmp_showers.end() );
      tmp_showers.clear();
    }

    //Clean showers without seed.
    std::vector<PandoraPlus::Calo1DCluster*> tmp_showers; tmp_showers.clear();
    for(int ic=0; ic<m_1dclusCol.size(); ic++) tmp_showers.push_back( const_cast<PandoraPlus::Calo1DCluster*>(m_1dShowerUCol[ic]) );
    m_1dShowerUCol.clear();
    for(int ic=0; ic<tmp_showers.size(); ic++){
      if(tmp_showers[ic]->getNseeds()==0){
        MergeToClosestCluster( tmp_showers[ic],tmp_showers );
        delete tmp_showers[ic]; tmp_showers[ic]=NULL;
        tmp_showers.erase(tmp_showers.begin()+ic);
        ic--;
      }
    }

    for(int as=0; as<tmp_showers.size(); as++){
      m_1dShowerUCol.push_back( tmp_showers[as] );
      m_datacol.bk_Cluster1DCol.push_back( tmp_showers[as] );
    }

    //  Longitudinal linking: update clusters' energy.
    std::vector<const CaloHalfCluster*> tmp_newClus; tmp_newClus.clear();
    LongitudinalLinking(m_1dShowerUCol, m_axisUCol, tmp_newClus);
    m_newClusUCol.insert(m_newClusUCol.end(), tmp_newClus.begin(), tmp_newClus.end());
    tmp_newClus.clear();
  }

  //Start in V direction
  m_newClusVCol.clear();
  for(int ih=0; ih<p_HalfClustersV->size(); ih++){
cout<<"    In HalfClusterV #"<<ih;
    m_axisVCol.clear();
    m_axisVCol = p_HalfClustersV->at(ih)->getAllHalfClusterCol();
cout<<", axis size "<<m_axisVCol.size();

    m_1dclusCol.clear();
    m_1dclusCol = p_HalfClustersV->at(ih)->getCluster();

cout<<", 1DClusters: size = "<<m_1dclusCol.size();

cout<<"  Print HoughAxisVCol: "<<endl;
for(int icl=0; icl<m_axisVCol.size(); icl++){
cout<<"In HoughAxis #"<<icl<<endl;
for(int i=0; i<m_axisVCol[icl]->getCluster().size(); i++)
  cout<<'\t'<<m_axisVCol[icl]->getCluster()[i]->getDlayer()<<'\t'<<m_axisVCol[icl]->getCluster()[i]->getEnergy()<<endl;
}

    for(int icl=0; icl<m_1dclusCol.size(); icl++){
      std::vector<const PandoraPlus::CaloUnit*> m_bars = m_1dclusCol[icl]->getBars();
    
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

cout<<"  Splitting in Cl #"<<icl<<", Dlayer = "<<m_1dclusCol[icl]->getDlayer();    
      //Split cluster to showers
      std::vector<const PandoraPlus::Calo1DCluster*> tmp_showers; tmp_showers.clear();
      ClusterSplitting( m_1dclusCol[icl], tmp_showers );
      //if(tmp_showers.size()==0) continue;
      m_1dShowerVCol.insert( m_1dShowerVCol.end(), tmp_showers.begin(), tmp_showers.end() );
cout<<", splitted shower size = "<<tmp_showers.size()<<endl;
      tmp_showers.clear();
    }
    //Clean showers without seed.
    std::vector<PandoraPlus::Calo1DCluster*> tmp_showers; tmp_showers.clear();
    for(int ic=0; ic<m_1dclusCol.size(); ic++) tmp_showers.push_back( const_cast<PandoraPlus::Calo1DCluster*>(m_1dShowerVCol[ic]) );
    m_1dShowerVCol.clear();
cout<<"  Clean showers without seed"<<endl;
    for(int ic=0; ic<tmp_showers.size(); ic++){
cout<<" shower #"<<ic<<", Dlayer = "<<tmp_showers[ic]->getDlayer()<<", Nseed = "<<tmp_showers[ic]->getNseeds()<<", En = "<<tmp_showers[ic]->getEnergy()<<endl;
      if(tmp_showers[ic]->getNseeds()==0){
        MergeToClosestCluster( tmp_showers[ic],tmp_showers );
        delete tmp_showers[ic]; tmp_showers[ic]=NULL;
        tmp_showers.erase(tmp_showers.begin()+ic);
        ic--;
      }
    }

cout<<"  Print Clenned showers: "<<endl;
for(int i=0; i<tmp_showers.size(); i++){
  printf("  Shower #%d: Dlayer %d, Nseed %d, En %.3f \n", i, tmp_showers[i]->getDlayer(), tmp_showers[i]->getNseeds(), tmp_showers[i]->getEnergy() );
}
    for(int as=0; as<tmp_showers.size(); as++){
      m_1dShowerVCol.push_back( tmp_showers[as] );
      m_datacol.bk_Cluster1DCol.push_back( tmp_showers[as] );
    }

    //  Longitudinal linking: update clusters' energy.
    std::vector<const CaloHalfCluster*> tmp_newClus; tmp_newClus.clear();
    LongitudinalLinking(m_1dShowerVCol, m_axisVCol, tmp_newClus);
    m_newClusVCol.insert(m_newClusVCol.end(), tmp_newClus.begin(), tmp_newClus.end());
    tmp_newClus.clear();
  }



  //Convert to non-const object
  std::vector<PandoraPlus::CaloHalfCluster*> tmp_nonconstClusU; tmp_nonconstClusU.clear();
  for(int icl=0; icl<m_newClusUCol.size(); icl++)
    tmp_nonconstClusU.push_back( const_cast<PandoraPlus::CaloHalfCluster*>(m_newClusUCol[icl]) );
  std::vector<PandoraPlus::CaloHalfCluster*> tmp_nonconstClusV; tmp_nonconstClusV.clear();
  for(int icl=0; icl<m_newClusVCol.size(); icl++)
    tmp_nonconstClusV.push_back( const_cast<PandoraPlus::CaloHalfCluster*>(m_newClusVCol[icl]) );


  m_datacol.map_HalfCluster[settings.map_stringPars["OutputClusName"]+"U"] = tmp_nonconstClusU; 
  m_datacol.map_HalfCluster[settings.map_stringPars["OutputClusName"]+"V"] = tmp_nonconstClusV; 

cout<<"Print HalfClusterU: "<<endl;
for(int i=0; i<m_newClusUCol.size(); i++){
  cout<<"  Total En = "<<m_newClusUCol[i]->getEnergy()<<endl;
  for(int is=0; is<m_newClusUCol[i]->getCluster().size(); is++)
    cout<<m_newClusUCol[i]->getCluster()[is]->getDlayer()<<'\t'<<m_newClusUCol[i]->getCluster()[is]->getEnergy()<<endl;
}
cout<<"Print HalfClusterV: "<<endl;
for(int i=0; i<m_newClusVCol.size(); i++){ 
  cout<<"  Total En = "<<m_newClusVCol[i]->getEnergy()<<endl;
  for(int is=0; is<m_newClusVCol[i]->getCluster().size(); is++)
    cout<<m_newClusVCol[i]->getCluster()[is]->getDlayer()<<'\t'<<m_newClusVCol[i]->getCluster()[is]->getEnergy()<<endl;
}


  tmp_nonconstClusU.clear(); tmp_nonconstClusV.clear();


  //Make tower 
  m_towerCol.clear(); 
  HalfClusterToTowers( m_newClusUCol, m_newClusVCol, m_towerCol);  

cout<<"  EnergySplittingAlg: Output tower size: "<<m_towerCol.size()<<", Print Energy: "<<endl;
for(int it=0; it<m_towerCol.size(); it++) cout<<m_towerCol[it]->getEnergy()<<endl;
cout<<endl;


  //  Convert to non-const
  std::vector<PandoraPlus::Calo3DCluster*> tmp_nonconstTower; tmp_nonconstTower.clear();
  for(int itw=0; itw<m_towerCol.size(); itw++)
    tmp_nonconstTower.push_back( const_cast<PandoraPlus::Calo3DCluster*>(m_towerCol[itw]) );
  
  m_datacol.map_CaloCluster[settings.map_stringPars["OutputTowerName"]] = tmp_nonconstTower;
  tmp_nonconstTower.clear();

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::ClearAlgorithm(){
  p_HalfClustersU = NULL;
  p_HalfClustersV = NULL;

  m_axisUCol.clear();
  m_axisVCol.clear();
  m_newClusUCol.clear();
  m_newClusVCol.clear();  
  m_1dclusCol.clear();
  m_1dShowerUCol.clear();
  m_1dShowerVCol.clear();
  m_towerCol.clear(); 

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::LongitudinalLinking( std::vector<const PandoraPlus::Calo1DCluster*>& m_showers, 
                                                    std::vector<const PandoraPlus::CaloHalfCluster*>& m_oldClusCol, 
                                                    std::vector<const PandoraPlus::CaloHalfCluster*>& m_newClusCol )
{
  if(m_showers.size()==0 || m_oldClusCol.size()==0) return StatusCode::SUCCESS;
  std::vector<PandoraPlus::CaloHalfCluster*> m_cluscol; m_cluscol.clear();

cout<<"  LongitudinalLinking: Print Input showers"<<endl;
for(int is=0; is<m_showers.size(); is++) 
  printf("    Dlayer: %d, Nseed: %d, Pos (%.2f, %.2f, %.2f), En %.3f \n", m_showers[is]->getDlayer(), m_showers[is]->getNseeds(), m_showers[is]->getPos().x(), m_showers[is]->getPos().y(), m_showers[is]->getPos().z(), m_showers[is]->getEnergy() );


  //Update old Hough clusters
  for(int ic=0; ic<m_oldClusCol.size(); ic++){
    PandoraPlus::CaloHalfCluster* m_newClus = new PandoraPlus::CaloHalfCluster();

    for(int is=0; is<m_oldClusCol[ic]->getCluster().size(); is++){
      const PandoraPlus::Calo1DCluster* m_shower = m_oldClusCol[ic]->getCluster()[is];

      const PandoraPlus::Calo1DCluster* m_selshower = NULL;
      bool fl_foundshower = false; 
      for(int js=0; js<m_showers.size(); js++){
        bool fl_inTower = false;
        for(int itw=0; itw<m_showers[js]->getTowerID().size() && !fl_inTower; itw++){
        for(int jtw=0; jtw<m_shower->getTowerID().size(); jtw++){
          if(m_showers[js]->getTowerID()[itw]==m_shower->getTowerID()[jtw]) {fl_inTower=true; break;}
        }}
cout<<"  Check axis finding: inTower = "<<fl_inTower<<endl;
printf("    Input shower: Dlayer %d, Nseed %d, Pos (%.2f, %.2f, %.2f), En %.3f \n", m_showers[js]->getDlayer(), m_showers[js]->getNseeds(), m_showers[js]->getPos().x(), m_showers[js]->getPos().y(), m_showers[js]->getPos().z(), m_showers[js]->getEnergy());
printf("    Axis shower: Dlayer %d, Dis %.2f \n", m_shower->getDlayer(), (m_showers[js]->getSeeds()[0]->getPosition()-m_shower->getPos()).Mag());

        if( fl_inTower &&
            m_showers[js]->getDlayer() == m_shower->getDlayer() && 
            m_showers[js]->getSeeds().size()==1 && 
            (m_showers[js]->getSeeds()[0]->getPosition()-m_shower->getPos()).Mag()<10 ) 
        {m_selshower = m_showers[js]; fl_foundshower=true; break; }
      }
cout<<"  Axis shower: Dlayer = "<<m_shower->getDlayer()<<", Find shower: "<<fl_foundshower<<endl;
      if(fl_foundshower && m_selshower!=NULL) m_newClus->addUnit( m_selshower );
    }

    m_newClus->setHoughPars(m_oldClusCol[ic]->getHoughAlpha(), m_oldClusCol[ic]->getHoughRho() );
    m_newClus->setIntercept(m_oldClusCol[ic]->getHoughIntercept());
    m_newClus->fitAxis("");
    m_cluscol.push_back( m_newClus );
  }

  std::vector<const PandoraPlus::Calo1DCluster*> m_leftshowers; m_leftshowers.clear(); 
  for(int ish=0; ish<m_showers.size(); ish++){
    bool fl_inclus = false; 
    for(int icl=0; icl<m_cluscol.size(); icl++){
      std::vector<const Calo1DCluster*> p_showerCol = m_cluscol[icl]->getCluster();
      if( find(p_showerCol.begin(), p_showerCol.end(), m_showers[ish])!=p_showerCol.end() ){ fl_inclus=true; break; }
    }

    if(!fl_inclus) m_leftshowers.push_back( m_showers[ish] );
  }


  //Merge showers into closest Hough cluster
  std::sort( m_leftshowers.begin(), m_leftshowers.end(), compLayer );
  for(int is=0; is<m_leftshowers.size(); is++) 
    MergeToClosestCluster( m_leftshowers[is], m_cluscol );


  //convert to const
  for(int icl=0; icl<m_cluscol.size(); icl++) m_newClusCol.push_back(m_cluscol[icl]);

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::HalfClusterToTowers( std::vector<const PandoraPlus::CaloHalfCluster*>& m_halfClusU,
                                                    std::vector<const PandoraPlus::CaloHalfCluster*>& m_halfClusV,
                                                    std::vector<const PandoraPlus::Calo3DCluster*>& m_towers ){

  m_towers.clear(); 

  std::map<std::vector<int>, std::vector<const PandoraPlus::Calo2DCluster*> > map_2DCluster;
  std::map<std::vector<int>, std::vector<const PandoraPlus::CaloHalfCluster*> > map_HalfClusterU; 
  std::map<std::vector<int>, std::vector<const PandoraPlus::CaloHalfCluster*> > map_HalfClusterV; 

  //Split CaloHalfClusterU
  for(int il=0; il<m_halfClusU.size(); il++){
    if(m_halfClusU[il]->getCluster().size()==0) {std::cout<<"WARNING: Have an empty CaloHalfCluster! Skip it! "<<std::endl; continue;}

cout<<"  HalfClusterU #"<<il<<": tower size = "<<m_halfClusU[il]->getTowerID().size()<<endl;

    //HalfCluster does not cover tower: 
    if( m_halfClusU[il]->getTowerID().size()==1 ){
      std::vector<int> cl_towerID =  m_halfClusU[il]->getTowerID()[0]; 
      map_HalfClusterU[cl_towerID].push_back(m_halfClusU[il]);
      continue; 
    }

    //CaloHalfCluster covers towers: Loop check showers. 
    std::map<std::vector<int>, PandoraPlus::CaloHalfCluster* > tmp_LongiClusMaps; tmp_LongiClusMaps.clear();
    for(int is=0; is<m_halfClusU[il]->getCluster().size(); is++){
      const PandoraPlus::Calo1DCluster* p_shower = m_halfClusU[il]->getCluster()[is];

      std::vector<std::vector<int>> towerIDs = p_shower->getTowerID();
      //Shower does not cover tower: 
      if(towerIDs.size()==1){
        if( tmp_LongiClusMaps.find( towerIDs[0] )!=tmp_LongiClusMaps.end() )
          tmp_LongiClusMaps[towerIDs[0]]->addUnit( p_shower );
        else{
          PandoraPlus::CaloHalfCluster* tmp_clus = new PandoraPlus::CaloHalfCluster();
          tmp_clus->addUnit( p_shower );
          tmp_LongiClusMaps[towerIDs[0]] = tmp_clus;
        }
        p_shower = nullptr;
        continue;
      }

      //Shower cover towers: split
      std::map< std::vector<int>, Calo1DCluster* > tmp_showerMap; tmp_showerMap.clear();      
      for(int ib=0; ib<p_shower->getBars().size(); ib++){
        std::vector<int> m_id(3);
        m_id[0] = p_shower->getBars()[ib]->getModule();
        m_id[1] = p_shower->getBars()[ib]->getPart();
        m_id[2] = p_shower->getBars()[ib]->getStave();

        if(tmp_showerMap.find(m_id)==tmp_showerMap.end()){
          Calo1DCluster* tmp_shower = new Calo1DCluster();
          tmp_shower->addUnit(p_shower->getBars()[ib]);
          tmp_shower->setIDInfo();
          tmp_showerMap[m_id] = tmp_shower;
        }
        else tmp_showerMap[m_id]->addUnit(p_shower->getBars()[ib]);
      }


      //Build connection between splitted showers. 
      for(auto &iter:tmp_showerMap){
        //iter.second->setSeed();
        for(auto &iter1:tmp_showerMap){
          if(iter==iter1) continue;
          iter.second->addCousinCluster( iter1.second );
        }
      }     


      for(auto &iter : tmp_showerMap){
        if(tmp_LongiClusMaps.find(iter.first) != tmp_LongiClusMaps.end())
          tmp_LongiClusMaps[iter.first]->addUnit(iter.second);
        else{
          PandoraPlus::CaloHalfCluster* tmp_clus = new PandoraPlus::CaloHalfCluster();
          tmp_clus->addUnit( iter.second );
          tmp_LongiClusMaps[iter.first] = tmp_clus;
        }
      }
      tmp_showerMap.clear();
      p_shower = nullptr;
    }

    //Connect cousins
    if(tmp_LongiClusMaps.size()>1){
      //std::map<std::vector<int>, PandoraPlus::CaloHalfCluster* >::iterator iter = tmp_LongiClusMaps.begin();
      //for(iter; iter!=tmp_LongiClusMaps.end(); iter++){
      //  std::map<std::vector<int>, PandoraPlus::CaloHalfCluster* >::iterator iter1 = tmp_LongiClusMaps.begin();
      //  for(iter1; iter1!=tmp_LongiClusMaps.end(); iter1++)
      //    if(iter!=iter1) iter->second->addCousinCluster(iter1->second);
      //}
      for(auto &iter : tmp_LongiClusMaps){
        for(auto &iter1 : tmp_LongiClusMaps){
          if(iter!= iter1) iter.second->addCousinCluster(iter1.second);
        }
      }
    }
    for(auto &iter : tmp_LongiClusMaps) map_HalfClusterU[iter.first].push_back(iter.second);
  }
cout<<"  After Splitting HalfClusterU: mapsize = "<<map_HalfClusterU.size()<<", Print En: "<<endl;
for(auto iter : map_HalfClusterU){ 
  printf("    Tower [%d, %d, %d], HalfClusterU size %d: ", iter.first[0], iter.first[1], iter.first[2], iter.second.size());
  for(int i=0; i<iter.second.size(); i++) cout<<iter.second[i]->getEnergy()<<'\t';
  cout<<endl;
}

  //Split CaloHalfClusterV
  for(int il=0; il<m_halfClusV.size(); il++){
    if(m_halfClusV[il]->getCluster().size()==0) {std::cout<<"WARNING: Have an empty CaloHalfCluster! Skip it! "<<std::endl; continue;}

    //HalfCluster does not cover tower:
    if( m_halfClusV[il]->getTowerID().size()==1 ){
      std::vector<int> cl_towerID =  m_halfClusV[il]->getTowerID()[0];
      map_HalfClusterV[cl_towerID].push_back(m_halfClusV[il]);
      continue;
    }

    //CaloHalfCluster covers towers: Loop check showers.
    std::map<std::vector<int>, PandoraPlus::CaloHalfCluster* > tmp_LongiClusMaps; tmp_LongiClusMaps.clear();
    for(int is=0; is<m_halfClusV[il]->getCluster().size(); is++){
      const PandoraPlus::Calo1DCluster* p_shower = m_halfClusV[il]->getCluster()[is];

      std::vector<std::vector<int>> towerIDs = p_shower->getTowerID();
      //Shower does not cover tower:
      if(towerIDs.size()==1){
        if( tmp_LongiClusMaps.find( towerIDs[0] )!=tmp_LongiClusMaps.end() )
          tmp_LongiClusMaps[towerIDs[0]]->addUnit( p_shower );
        else{
          PandoraPlus::CaloHalfCluster* tmp_clus = new PandoraPlus::CaloHalfCluster();
          tmp_clus->addUnit( p_shower );
          tmp_LongiClusMaps[towerIDs[0]] = tmp_clus;
        }
        p_shower = nullptr;
        continue;
      }

      //Shower cover towers: split
      std::map< std::vector<int>, Calo1DCluster* > tmp_showerMap; tmp_showerMap.clear();
      for(int ib=0; ib<p_shower->getBars().size(); ib++){
        std::vector<int> m_id(3);
        m_id[0] = p_shower->getBars()[ib]->getModule();
        m_id[1] = p_shower->getBars()[ib]->getPart();
        m_id[2] = p_shower->getBars()[ib]->getStave();

        if(tmp_showerMap.find(m_id)==tmp_showerMap.end()){
          Calo1DCluster* tmp_shower = new Calo1DCluster();
          tmp_shower->addUnit(p_shower->getBars()[ib]);
          tmp_shower->setIDInfo();
          tmp_showerMap[m_id] = tmp_shower;
        }
        else tmp_showerMap[m_id]->addUnit(p_shower->getBars()[ib]);
      }

      //Build connection between splitted showers.
      for(auto &iter:tmp_showerMap){
        //iter.second->setSeed();
        for(auto &iter1:tmp_showerMap){
          if(iter==iter1) continue;
          iter.second->addCousinCluster( iter1.second );
        }
      }

      for(auto &iter : tmp_showerMap){
        if(tmp_LongiClusMaps.find(iter.first) != tmp_LongiClusMaps.end())
          tmp_LongiClusMaps[iter.first]->addUnit(iter.second);
        else{
          PandoraPlus::CaloHalfCluster* tmp_clus = new PandoraPlus::CaloHalfCluster();
          tmp_clus->addUnit( iter.second );
          tmp_LongiClusMaps[iter.first] = tmp_clus;
        }
      }
      tmp_showerMap.clear();
      p_shower = nullptr;
    }

    //Connect cousins
    if(tmp_LongiClusMaps.size()>1){
      //std::map<std::vector<int>, PandoraPlus::CaloHalfCluster* >::iterator iter = tmp_LongiClusMaps.begin();
      //for(iter; iter!=tmp_LongiClusMaps.end(); iter++){
      //  std::map<std::vector<int>, PandoraPlus::CaloHalfCluster* >::iterator iter1 = tmp_LongiClusMaps.begin();
      //  for(iter1; iter1!=tmp_LongiClusMaps.end(); iter1++)
      //    if(iter!=iter1) iter->second->addCousinCluster(iter1->second);
      //}
      for(auto &iter : tmp_LongiClusMaps){
        for(auto &iter1 : tmp_LongiClusMaps){
          if(iter!= iter1) iter.second->addCousinCluster(iter1.second);
        }
      }
    }
    for(auto &iter : tmp_LongiClusMaps) map_HalfClusterV[iter.first].push_back(iter.second);
  }
cout<<"  After Splitting HalfClusterV: mapsize = "<<map_HalfClusterV.size()<<", Print En: "<<endl;
for(auto iter : map_HalfClusterV){
  printf("    Tower [%d, %d, %d], HalfClusterV size %d: ", iter.first[0], iter.first[1], iter.first[2], iter.second.size());
  for(int i=0; i<iter.second.size(); i++) cout<<iter.second[i]->getEnergy()<<'\t';
  cout<<endl;
}

  //Build 2DCluster
  for(auto &iterU : map_HalfClusterU){
    if( map_HalfClusterV.find(iterU.first)==map_HalfClusterV.end() ) continue; 

printf("  Building 2DCluster in Tower: [%d, %d, %d] \n", iterU.first[0], iterU.first[1], iterU.first[2] );

    std::vector<const PandoraPlus::CaloHalfCluster*> p_halfClusU = iterU.second; 
    std::vector<const PandoraPlus::CaloHalfCluster*> p_halfClusV = map_HalfClusterV[iterU.first]; 
cout<<"  Get HalfCluster: "<<p_halfClusU.size()<<", "<<p_halfClusV.size()<<endl;

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

cout<<"  Print ordered shower map in U: "<<endl;
for(auto iter : m_orderedShowerU) cout<<iter.first<<'\t'<<iter.second.size()<<endl;
cout<<"  Print ordered shower map in V: "<<endl;
for(auto iter : m_orderedShowerV) cout<<iter.first<<'\t'<<iter.second.size()<<endl;


    //Create super-layers (block)
    std::vector<const PandoraPlus::Calo2DCluster*> m_blocks; m_blocks.clear();
    for(auto &iter1 : m_orderedShowerU){
      if( m_orderedShowerV.find( iter1.first )==m_orderedShowerV.end() ) continue; 
cout<<"    Make Block in Layer #"<<iter1.first;
      PandoraPlus::Calo2DCluster* tmp_block = new PandoraPlus::Calo2DCluster();
      tmp_block->addTowerID( iterU.first );
      for(int is=0; is<iter1.second.size(); is++) tmp_block->addUnit( iter1.second.at(is) );
      for(int is=0; is<m_orderedShowerV[iter1.first].size(); is++) tmp_block->addUnit( m_orderedShowerV[iter1.first].at(is) );
cout<<", Block En = "<<tmp_block->getEnergy()<<endl;
      m_blocks.push_back( tmp_block );
    }
    map_2DCluster[iterU.first] = m_blocks;
  }
cout<<"  2D cluster map size: "<<map_2DCluster.size()<<", Print: "<<endl;
for(auto iter:map_2DCluster){
  printf("    Tower: [%d, %d, %d], 2DCluster size %d: \n", iter.first[0], iter.first[1], iter.first[2], iter.second.size() );
  for(int i=0; i<iter.second.size(); i++) cout<<iter.second[i]->getDlayer()<<'\t'<<iter.second[i]->getEnergy()<<endl;
  cout<<endl;
}

  //Form a tower:
  for(auto &iter : map_2DCluster){
    std::vector<int> m_towerID = iter.first;

    PandoraPlus::Calo3DCluster* m_tower = new PandoraPlus::Calo3DCluster();
    m_tower->addTowerID( m_towerID );
    for(int i2d=0; i2d<map_2DCluster[m_towerID].size(); i2d++) m_tower->addUnit(map_2DCluster[m_towerID][i2d]);
    m_tower->setHalfClusters( settings.map_stringPars["OutputClusName"]+"U", map_HalfClusterU[m_towerID], 
                              settings.map_stringPars["OutputClusName"]+"V", map_HalfClusterV[m_towerID] );
    m_towers.push_back(m_tower);
  }

//cout<<"  After splitting: tower size "<<m_towerCol.size()<<endl;

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::ClusterSplitting( const PandoraPlus::Calo1DCluster* m_cluster, std::vector<const PandoraPlus::Calo1DCluster*>& outshCol ){
//cout<<"ClusterSplitting: input cluster seed size = "<<m_cluster->getSeeds().size()<<endl;
//cout<<"Seed position: ";
//for(int a=0; a<m_cluster->getNseeds(); a++) printf(" (%.2f, %.2f, %.2f) \t", m_cluster->getSeeds()[a]->getPosition().x(), 
//                                                                             m_cluster->getSeeds()[a]->getPosition().y(),
//                                                                             m_cluster->getSeeds()[a]->getPosition().z() );
//cout<<endl;

  //No seed in cluster: return empty vector.
  if(m_cluster->getNseeds()==0) { 
    outshCol.push_back(m_cluster);
    //std::cout<<"WARNING: Still have no-seed cluster!!"<<std::endl; 
    return StatusCode::SUCCESS; 
  }

  //1 seed or second moment less than threshold: Not split. Turn cluster to shower and return
  else if(m_cluster->getNseeds()<2 || m_cluster->getScndMoment()<settings.map_floatPars["th_split"]){
    //PandoraPlus::Calo1DCluster* shower = new PandoraPlus::Calo1DCluster();
    //shower->setBars(m_cluster->getBars());
    //if(m_cluster->getNseeds()!=0) shower->addSeed(m_cluster->getSeeds()[0]);
    //else shower->addSeed(nullptr);
    //shower->setIDInfo();
    outshCol.push_back(m_cluster);
    return StatusCode::SUCCESS;
  }
  
  //2 or more seeds, large second moment: Split
  int Nshower = m_cluster->getNseeds();
  int Nbars = m_cluster->getBars().size();
  double Eseed[Nshower] = {0};
  double weight[Nbars][Nshower] = {0};
  TVector3 SeedPos[Nshower];
  TVector3 SeedPos_Origin[Nshower];
  for(int is=0;is<Nshower;is++){ SeedPos[is] = m_cluster->getSeeds()[is]->getPosition(); SeedPos_Origin[is]=SeedPos[is]; }
  CalculateInitialEseed(m_cluster->getSeeds(), SeedPos, Eseed);

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

//cout<<"Iter #"<<iter<<": Seed pos ";
//for(int a=0; a<Nshower; a++) printf(" (%.2f, %.2f, %.2f) \t", SeedPos[a].x(), SeedPos[a].y(), SeedPos[a].z());
//cout<<endl;
  }
  while(iter<20 && !isConverge && !isShifted );
  if(iter>=20){
    std::cout<<"WARNING: Iteration time larger than 20! Might not converge!"<<std::endl;
    std::cout<<"  For Check: NBars: "<<m_cluster->getBars().size()<<"  Nseeds: "<<Nshower<<std::endl;
  }

//cout<<"After weight calculation. "<<endl;
  for(int is=0;is<Nshower;is++) SeedPos[is] = m_cluster->getSeeds()[is]->getPosition();
  for(int is=0;is<Nshower;is++){
    PandoraPlus::Calo1DCluster* shower = new PandoraPlus::Calo1DCluster();
    std::vector<const PandoraPlus::CaloUnit*> Bars; Bars.clear();
    int iseed=-1;
    double _Emax = -99;
    
    for(int ib=0;ib<Nbars;ib++){
      PandoraPlus::CaloUnit* bar = m_cluster->getBars()[ib]->Clone();

      bar->setQ( bar->getQ1()*weight[ib][is], bar->getQ2()*weight[ib][is]  );
      if( bar->getEnergy()>_Emax ) { _Emax=bar->getEnergy(); iseed=ib; }
      Bars.push_back(bar);
    }
    if(iseed<0) { std::cout<<"ERROR: Can not find seed(max energy bar) in this shower! Please Check!"<<std::endl; iseed=0;}
    if( (Bars[iseed]->getPosition()-SeedPos[is]).Mag()>15 ) { std::cout<<"ERROR: MaxEnergy bar is too far with original seed! Please Check! iSeed = "<<iseed<<std::endl; }

    shower->setBars(Bars);
    shower->addSeed(Bars[iseed]);
    shower->setIDInfo(); 

/*
cout<<"  Print shower #"<<is<<endl;
printf("    Shower Nbars = %d, pos/E (%.2f, %.2f, %.2f, %.3f) \n", shower->getBars().size(), 
                                                                   shower->getPos().x(), shower->getPos().y(), shower->getPos().z(), shower->getEnergy() );
cout<<"    Bars: "<<endl;
for(int a=0; a<shower->getBars().size(); a++) 
  printf("      Bar pos/E (%.2f, %.2f, %.2f, %.3f), ID [%d, %d, %d] \n",  shower->getBars()[a]->getPosition().x(), 
                                                         shower->getBars()[a]->getPosition().y(),
                                                         shower->getBars()[a]->getPosition().z(),
                                                         shower->getBars()[a]->getEnergy(), 
                                                         shower->getBars()[a]->getModule(),
                                                         shower->getBars()[a]->getPart(),
                                                         shower->getBars()[a]->getStave() );
*/


    outshCol.push_back(shower);
  }
  
  return StatusCode::SUCCESS;
}

StatusCode EnergySplittingAlg::MergeToClosestCluster( PandoraPlus::Calo1DCluster* iclus, std::vector<PandoraPlus::Calo1DCluster*>& clusvec ){

  int cLedge = iclus->getLeftEdge();
  int cRedge = iclus->getRightEdge();

//printf("ClusterMerging: input cluster edge [%d, %d]\n", cLedge, cRedge);

  //Find the closest cluster with iclus.
  int minD = 99;
  int index = -1;
  for(int icl=0; icl<clusvec.size(); icl++){
    if(clusvec[icl]->getNseeds()==0) continue;
    if( iclus->getTowerID() != clusvec[icl]->getTowerID() ) continue;
    if( iclus->getDlayer() != clusvec[icl]->getDlayer() ) continue;

    int iLedge = clusvec[icl]->getLeftEdge();
    int iRedge = clusvec[icl]->getRightEdge();

    //int dis = (cLedge-iRedge>0 ? cLedge-iRedge : iLedge-cRedge );
    int dis = min( abs(cLedge-iRedge), abs(iLedge-cRedge) );

//printf("  Loop in clusterCol #%d: range [%d, %d], distance %d, energy ratio %.3f \n", icl, iLedge, iRedge, dis, iclus->getEnergy()/clusvec[icl]->getEnergy());
    if(dis>10) continue; //Don't merge to a too far cluster.
    if(dis<minD){ minD = dis; index=icl; }
  }
//cout<<"Selected closest cluster #"<<index<<endl;

  if(index<0) return StatusCode::FAILURE;

  //Merge to the selected cluster
  for(int icl=0; icl<iclus->getBars().size(); icl++)  clusvec[index]->addUnit(iclus->getBars()[icl]);


  return StatusCode::SUCCESS;
}

StatusCode EnergySplittingAlg::MergeToClosestCluster( const PandoraPlus::Calo1DCluster* m_shower, 
                                                      std::vector<PandoraPlus::CaloHalfCluster*>& m_clusters ){
  if(m_clusters.size()==0) return StatusCode::SUCCESS;


  int minLayer = 99;
  int maxLayer = -99;
  for(int ic=0; ic<m_clusters.size(); ic++){
    if(minLayer>m_clusters[ic]->getBeginningDlayer())  minLayer = m_clusters[ic]->getBeginningDlayer();
    if(maxLayer<m_clusters[ic]->getEndDlayer())        maxLayer = m_clusters[ic]->getEndDlayer();
  }
  if(minLayer==99 || maxLayer<0) return StatusCode::SUCCESS;

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
      m_distance = (sh_pos-m_clusters[ic]->getClusterInLayer(m_clusters[ic]->getBeginningDlayer())[0]->getPos()).Mag();     
    else if(dlayer>=maxLayer)
      m_distance = (sh_pos-m_clusters[ic]->getClusterInLayer(m_clusters[ic]->getEndDlayer())[0]->getPos()).Mag();    
    else
      m_distance = (sh_pos-m_clusters[ic]->getPos()).Mag();

//cout<<"    Cluster #"<<ic<<": distance with shower = "<<m_distance<<endl;
    
    if( m_distance<minR )  { minR=m_distance; index_cluster=ic; }
  }

//printf("  minR = %.3f, in #cl %d \n", minR, index_cluster);

  if(index_cluster>=0) m_clusters[index_cluster]->addUnit( m_shower );
  //delete m_shower; m_shower = NULL;

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

}


double EnergySplittingAlg::GetShowerProfile(const TVector3& p_bar, const TVector3& p_seed ){
  TVector3 rpos = p_bar-p_seed;
  double dis = sqrt(rpos.Mag2());
  double a1, a2, b1, b2, Rm;
  a1=0.037; a2=0.265; b1=0.101; b2=0.437; Rm=1.868;

  return a1*exp(-b1*dis/Rm)+a2*exp(-b2*dis/Rm);
}
#endif
