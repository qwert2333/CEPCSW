#ifndef _ENERGYSPLITTING_ALG_C
#define _ENERGYSPLITTING_ALG_C

#include "Algorithm/EnergySplittingAlg.h"

StatusCode EnergySplittingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  if(settings.map_floatPars.find("th_split")==settings.map_floatPars.end()) settings.map_floatPars["th_split"] = -1;
  if(settings.map_floatPars.find("Eth_Seed")==settings.map_floatPars.end()) settings.map_floatPars["Eth_Seed"] = 0.005;
  if(settings.map_floatPars.find("Eth_ShowerAbs")==settings.map_floatPars.end()) settings.map_floatPars["Eth_Shower"] = 0.005;
  if(settings.map_stringPars.find("OutputLongiClusName")==settings.map_stringPars.end()) settings.map_stringPars["OutputLongiClusName"] = "ESLongiCluster";

  return StatusCode::SUCCESS;
};


StatusCode EnergySplittingAlg::Initialize(){

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  //Input: HalfCluster (with 1D clusters and HoughAxis)
  //Output: Create Towers for the matching. 

  std::vector<PandoraPlus::CaloHalfCluster*>* p_HalfClusters = &(m_datacol.map_LongiCluster["HalfClusterCol"]);


      //Split cluster to showers
      std::vector<const PandoraPlus::Calo1DCluster*> m_barShowerUCol; m_barShowerUCol.clear();
      for(int ic=0; ic<m_barClusUCol.size(); ic++){
        std::vector<const PandoraPlus::Calo1DCluster*> m_showers; m_showers.clear();
        ClusterSplitting( m_barClusUCol[ic], m_showers );
        if(m_showers.size()==0) continue;
        m_barShowerUCol.insert( m_barShowerUCol.end(), m_showers.begin(), m_showers.end() );
      }
      std::vector<const PandoraPlus::Calo1DCluster*> m_barShowerVCol; m_barShowerVCol.clear();
      for(int ic=0; ic<m_barClusVCol.size(); ic++){
        std::vector<const PandoraPlus::Calo1DCluster*> m_showers; m_showers.clear();
        ClusterSplitting( m_barClusVCol[ic], m_showers );
        if(m_showers.size()==0) continue;
        m_barShowerVCol.insert( m_barShowerVCol.end(), m_showers.begin(), m_showers.end() );
      }

  std::vector<PandoraPlus::Calo3DCluster*>* p_3DClusters = &(m_datacol.Cluster3DCol);
  if( !p_3DClusters || p_3DClusters->size()==0){ std::cout<<"Warning: Empty or invalid input in EnergySplittingAlg! Please check previous algorithm!"<<std::endl; return StatusCode::SUCCESS; }

  for(int it=0; it<p_3DClusters->size(); it++){

    //Need to transfer const blocks to mutable. 
    std::vector<PandoraPlus::Calo2DCluster*> m_2DclusCol; m_2DclusCol.clear();
    for(int ib=0; ib<p_3DClusters->at(it)->getCluster().size(); ib++ )
      m_2DclusCol.push_back( const_cast<PandoraPlus::Calo2DCluster *>(p_3DClusters->at(it)->getCluster()[ib]) );

    //Get all LongiClusters from maps. 
    std::vector<const PandoraPlus::LongiCluster*> m_LongiClusUCol; m_LongiClusUCol.clear(); 
    std::vector<const PandoraPlus::LongiCluster*> m_LongiClusVCol; m_LongiClusVCol.clear(); 
    std::map<std::string, std::vector<const PandoraPlus::LongiCluster*> > tmp_mapU = p_3DClusters->at(it)->getLongiClusterUMap();
    std::map<std::string, std::vector<const PandoraPlus::LongiCluster*> > tmp_mapV = p_3DClusters->at(it)->getLongiClusterVMap();
    for(auto &iter : tmp_mapU){
      for(int icl=0; icl<iter.second.size(); icl++) m_LongiClusUCol.push_back(iter.second[icl]);
    }
    for(auto &iter : tmp_mapV){
      for(int icl=0; icl<iter.second.size(); icl++) m_LongiClusVCol.push_back(iter.second[icl]);
    }
    tmp_mapU.clear(); tmp_mapV.clear(); 
    if(m_LongiClusUCol.size()==0 || m_LongiClusVCol.size()==0) continue;


    //Get the layer range of LongiClusters
    int minLayerU = 99;
    int minLayerV = 99;
    int maxLayerU = -99;
    int maxLayerV = -99;
    for(int ic=0; ic<m_LongiClusUCol.size(); ic++){
      if(minLayerU>m_LongiClusUCol[ic]->getBeginningDlayer()) minLayerU = m_LongiClusUCol[ic]->getBeginningDlayer();
      if(maxLayerU<m_LongiClusUCol[ic]->getEndDlayer()) maxLayerU = m_LongiClusUCol[ic]->getEndDlayer();
    }
    for(int ic=0; ic<m_LongiClusVCol.size(); ic++){
      if(minLayerV>m_LongiClusVCol[ic]->getBeginningDlayer()) minLayerV = m_LongiClusVCol[ic]->getBeginningDlayer();
      if(maxLayerV<m_LongiClusVCol[ic]->getEndDlayer()) maxLayerV = m_LongiClusVCol[ic]->getEndDlayer();
    }
printf("  Layer range from LongiClusters: [%d, %d] / [%d, %d] \n", minLayerU, maxLayerU, minLayerV, maxLayerV);


    //Make clusters and shwoers in each layer:
    for(int ib=0; ib<m_2DclusCol.size(); ib++ ){

      std::vector<const PandoraPlus::CaloUnit*> m_barColU = m_2DclusCol[ib]->getBarUCol();
      std::vector<const PandoraPlus::CaloUnit*> m_barColV = m_2DclusCol[ib]->getBarVCol();



      //Neighbor clustering
      std::vector<PandoraPlus::Calo1DCluster*> m_barClusUCol; m_barClusUCol.clear();
      std::vector<PandoraPlus::Calo1DCluster*> m_barClusVCol; m_barClusVCol.clear();
      int dlayer = m_2DclusCol[ib]->getDlayer();

//double totE1 = 0;
//double totE2 = 0;
//for(int i=0; i<m_barColU.size(); i++) totE1 += m_barColU[i]->getEnergy();
//for(int i=0; i<m_barColV.size(); i++) totE2 += m_barColV[i]->getEnergy();
//printf("  Raw bar energy in Layer %d: U %.4f, V %.4f. \n", dlayer, totE1, totE2);

      if( dlayer>maxLayerU || dlayer<minLayerU ) Clustering( m_barColU, m_barClusUCol );
      else Clustering( m_barColU, m_barClusUCol, m_LongiClusUCol );

      if( dlayer>maxLayerV || dlayer<minLayerV ) Clustering( m_barColV, m_barClusVCol );
      else Clustering( m_barColV, m_barClusVCol, m_LongiClusVCol );

//printf("  Neighbor clustering in Layer %d: cluster energy U:", dlayer);
//for(int a=0; a<m_barClusUCol.size(); a++) cout<<m_barClusUCol[a]->getEnergy()<<", ";
//cout<<". V: ";
//for(int a=0; a<m_barClusVCol.size(); a++) cout<<m_barClusVCol[a]->getEnergy()<<", ";
//cout<<endl;

      m_datacol.bk_Cluster1DCol.insert( m_datacol.bk_Cluster1DCol.end(), m_barClusUCol.begin(), m_barClusUCol.end() );
      m_datacol.bk_Cluster1DCol.insert( m_datacol.bk_Cluster1DCol.end(), m_barClusVCol.begin(), m_barClusVCol.end() );


      //Split cluster to showers
      std::vector<const PandoraPlus::Calo1DCluster*> m_barShowerUCol; m_barShowerUCol.clear();
      for(int ic=0; ic<m_barClusUCol.size(); ic++){
        std::vector<const PandoraPlus::Calo1DCluster*> m_showers; m_showers.clear();
        ClusterSplitting( m_barClusUCol[ic], m_showers );
        if(m_showers.size()==0) continue;
        m_barShowerUCol.insert( m_barShowerUCol.end(), m_showers.begin(), m_showers.end() );
      }
      std::vector<const PandoraPlus::Calo1DCluster*> m_barShowerVCol; m_barShowerVCol.clear();
      for(int ic=0; ic<m_barClusVCol.size(); ic++){
        std::vector<const PandoraPlus::Calo1DCluster*> m_showers; m_showers.clear();
        ClusterSplitting( m_barClusVCol[ic], m_showers );
        if(m_showers.size()==0) continue;
        m_barShowerVCol.insert( m_barShowerVCol.end(), m_showers.begin(), m_showers.end() );
      }

//printf("  After splitting in Layer %d: shower energy U:", dlayer);
//for(int a=0; a<m_barShowerUCol.size(); a++) cout<<m_barShowerUCol[a]->getEnergy()<<", ";
//cout<<". V: ";
//for(int a=0; a<m_barShowerVCol.size(); a++) cout<<m_barShowerVCol[a]->getEnergy()<<", ";
//cout<<endl;
//cout<<endl;

      m_2DclusCol[ib]->setShowerUCol( m_barShowerUCol );
      m_2DclusCol[ib]->setShowerVCol( m_barShowerVCol );     

      for(int is=0; is<m_barShowerUCol.size(); is++) m_datacol.bk_Cluster1DCol.push_back( const_cast<PandoraPlus::Calo1DCluster*>(m_barShowerUCol[is]) );
      for(int is=0; is<m_barShowerVCol.size(); is++) m_datacol.bk_Cluster1DCol.push_back( const_cast<PandoraPlus::Calo1DCluster*>(m_barShowerVCol[is]) );

    }//End loop 2DCluster (Layer)


/*
cout<<"  After Clustering and cluster splitting in each layer: check BarShowers in Layers:"<<endl;
for(int a=0; a<m_2DclusCol.size(); a++){
  printf("    In 2DClus #%d: Layer %d, tower size: %d, tower ID ", a, m_2DclusCol[a]->getDlayer(), m_2DclusCol[a]->getTowerID().size() );
  for(int id=0; id<m_2DclusCol[a]->getTowerID().size(); id++) 
    printf("[%d, %d, %d], ", m_2DclusCol[a]->getTowerID()[id][0], m_2DclusCol[a]->getTowerID()[id][1], m_2DclusCol[a]->getTowerID()[id][2]);
  cout<<endl;
  cout<<"    Cluster U energy: ";
  for(int ic=0; ic<m_2DclusCol[a]->getClusterU().size(); ic++) cout<<m_2DclusCol[a]->getClusterU()[ic]->getEnergy()<<"  ";
  cout<<endl;
  cout<<"    Shower U energy:  ";
  for(int ic=0; ic<m_2DclusCol[a]->getShowerUCol().size(); ic++) cout<<m_2DclusCol[a]->getShowerUCol()[ic]->getEnergy()<<"  ";
  cout<<endl;
  cout<<"    Cluster V energy: ";
  for(int ic=0; ic<m_2DclusCol[a]->getClusterV().size(); ic++) cout<<m_2DclusCol[a]->getClusterV()[ic]->getEnergy()<<"  ";
  cout<<endl;
  cout<<"    Shower V energy:  ";
  for(int ic=0; ic<m_2DclusCol[a]->getShowerVCol().size(); ic++) cout<<m_2DclusCol[a]->getShowerVCol()[ic]->getEnergy()<<"  ";
  cout<<endl;
}
*/
    //Longitudinal linking: update clusters' energy.

    std::vector<const LongiCluster*> m_newLongiClusUCol; m_newLongiClusUCol.clear();
    std::vector<const LongiCluster*> m_newLongiClusVCol; m_newLongiClusVCol.clear();
    LongiClusterLinking(m_2DclusCol, m_LongiClusUCol, m_newLongiClusUCol);
    LongiClusterLinking(m_2DclusCol, m_LongiClusVCol, m_newLongiClusVCol);

/*
cout<<"  2DClusters in alg: size = "<<m_2DclusCol.size()<<endl;
for(int a=0; a<m_2DclusCol.size(); a++){
  printf("    In 2DClus #%d: Layer %d, tower size: %d, tower ID ", a, m_2DclusCol[a]->getDlayer(), m_2DclusCol[a]->getTowerID().size() );
  for(int id=0; id<m_2DclusCol[a]->getTowerID().size(); id++)
    printf("[%d, %d, %d], ", m_2DclusCol[a]->getTowerID()[id][0], m_2DclusCol[a]->getTowerID()[id][1], m_2DclusCol[a]->getTowerID()[id][2]);
  cout<<endl;
  printf("    BarShowerSize (%d, %d) \n", m_2DclusCol[a]->getShowerUCol().size(), m_2DclusCol[a]->getShowerVCol().size());
  cout<<"      ";
  for(int is=0; is<m_2DclusCol[a]->getShowerUCol().size(); is++)
    printf("[ (%.3f, %.3f, %.3f, %.3f), %d, %p ], ", m_2DclusCol[a]->getShowerUCol()[is]->getPos().x(), 
                                                     m_2DclusCol[a]->getShowerUCol()[is]->getPos().y(), 
                                                     m_2DclusCol[a]->getShowerUCol()[is]->getPos().z(),  
                                                     m_2DclusCol[a]->getShowerUCol()[is]->getEnergy(),  
                                                     m_2DclusCol[a]->getShowerUCol()[is]->getBars().size(),  
                                                     m_2DclusCol[a]->getShowerUCol()[is]  );
  cout<<endl;
  cout<<"      ";
  for(int is=0; is<m_2DclusCol[a]->getShowerVCol().size(); is++)
    printf("[ (%.3f, %.3f, %.3f, %.3f), %d, %p ], ", m_2DclusCol[a]->getShowerVCol()[is]->getPos().x(),
                                                     m_2DclusCol[a]->getShowerVCol()[is]->getPos().y(),
                                                     m_2DclusCol[a]->getShowerVCol()[is]->getPos().z(),
                                                     m_2DclusCol[a]->getShowerVCol()[is]->getEnergy(),
                                                     m_2DclusCol[a]->getShowerVCol()[is]->getBars().size(),
                                                     m_2DclusCol[a]->getShowerVCol()[is]  );
  cout<<endl;

}
*/
/*
cout<<"  Current 3DCluster E = "<<p_3DClusters->at(it)->getEnergy()<<endl;
cout<<"  2DClusters in 3DCluster: size "<<p_3DClusters->at(it)->getCluster().size()<<endl;
for(int a=0; a<p_3DClusters->at(it)->getCluster().size(); a++){
  const Calo2DCluster* p_clus = p_3DClusters->at(it)->getCluster()[a];
  printf("    In 2DClus #%d: Layer %d, Energy %.3f (%.3f), tower size: %d, tower ID ", a, p_clus->getDlayer(), p_clus->getEnergy(), p_clus->getEnergy1(), p_clus->getTowerID().size() );
  for(int id=0; id<p_clus->getTowerID().size(); id++)
    printf("[%d, %d, %d], ", p_clus->getTowerID()[id][0], p_clus->getTowerID()[id][1], p_clus->getTowerID()[id][2]);
  cout<<endl;
  printf("    BarShowerSize (%d, %d) \n", p_clus->getShowerUCol().size(), p_clus->getShowerVCol().size());
  cout<<"      ";
  for(int is=0; is<p_clus->getShowerUCol().size(); is++)
    printf("[ (%.3f, %.3f, %.3f, %.3f), %d, %d, %p ], ", p_clus->getShowerUCol()[is]->getPos().x(),
                                                     p_clus->getShowerUCol()[is]->getPos().y(),
                                                     p_clus->getShowerUCol()[is]->getPos().z(),
                                                     p_clus->getShowerUCol()[is]->getEnergy(),
                                                     p_clus->getShowerUCol()[is]->getBars().size(),
                                                     p_clus->getShowerUCol()[is]->getTowerID().size(),
                                                     p_clus->getShowerUCol()[is]  );
  cout<<endl;
  cout<<"      ";
  for(int is=0; is<p_clus->getShowerVCol().size(); is++)
    printf("[ (%.3f, %.3f, %.3f, %.3f), %d, %d, %p ], ", p_clus->getShowerVCol()[is]->getPos().x(),
                                                     p_clus->getShowerVCol()[is]->getPos().y(),
                                                     p_clus->getShowerVCol()[is]->getPos().z(),
                                                     p_clus->getShowerVCol()[is]->getEnergy(),
                                                     p_clus->getShowerVCol()[is]->getBars().size(),
                                                     p_clus->getShowerVCol()[is]->getTowerID().size(),
                                                     p_clus->getShowerVCol()[is]  );
  cout<<endl;

}
*/

/*
cout<<"  After LongiClusterLinking: check LongiClusterU "<<endl;
for(int a=0; a<m_newLongiClusUCol.size(); a++){
  printf("    LongiCluster #%d BarShower size: %d, address %p \n", a, m_newLongiClusUCol[a]->getBarShowers().size(), m_newLongiClusUCol[a]);
  for(int is=0; is<m_newLongiClusUCol[a]->getBarShowers().size(); is++)
    printf("      Layer %d, address %p \n", m_newLongiClusUCol[a]->getBarShowers()[is]->getDlayer(), m_newLongiClusUCol[a]->getBarShowers()[is]);
}
cout<<"  After LongiClusterLinking: check LongiClusterV "<<endl;
for(int a=0; a<m_newLongiClusVCol.size(); a++){
  printf("    LongiCluster #%d BarShower size: %d, address %p \n", a, m_newLongiClusVCol[a]->getBarShowers().size(), m_newLongiClusVCol[a]);
  for(int is=0; is<m_newLongiClusVCol[a]->getBarShowers().size(); is++)
    printf("      Layer %d, address %p \n", m_newLongiClusVCol[a]->getBarShowers()[is]->getDlayer(), m_newLongiClusVCol[a]->getBarShowers()[is]);
}
*/

    p_3DClusters->at(it)->setLongiClusters( "tmpCluster", m_newLongiClusUCol, 
                                            "tmpCluster", m_newLongiClusVCol );

    //Split 3DCluster to towers for later matching. 
    Split3DClusterToTowers( p_3DClusters->at(it) );  
/*
cout<<"  After Cluster splitting: tower size "<<p_3DClusters->at(it)->getTowers().size()<<endl;
for(int itw=0; itw<p_3DClusters->at(it)->getTowers().size(); itw++){
const Calo3DCluster* p_tower = p_3DClusters->at(it)->getTowers()[itw];
cout<<"  In tower #"<<itw<<": Energy "<<p_tower->getEnergy()<<", towerID size "<<p_tower->getTowerID().size()<<endl;
printf("    TowerID: [%d, %d, %d], \n", p_tower->getTowerID()[0][0], p_tower->getTowerID()[0][1], p_tower->getTowerID()[0][2]);

printf("    2DCluster size: %d, print 2DClusters: \n", p_tower->getCluster().size());
for(int ib=0; ib<p_tower->getCluster().size(); ib++){
  printf("      Layer %d, Energy %.3f, Shower size (%d, %d) \n", p_tower->getCluster()[ib]->getDlayer(), p_tower->getCluster()[ib]->getEnergy(), p_tower->getCluster()[ib]->getShowerUCol().size(), p_tower->getCluster()[ib]->getShowerVCol().size() );
  cout<<"      ShowerU pos+E, Nbars, Nseed, cousin size and address: "<<endl;
  for(int is=0; is<p_tower->getCluster()[ib]->getShowerUCol().size(); is++){
    const Calo1DCluster* p_shower = p_tower->getCluster()[ib]->getShowerUCol()[is];
    printf("        #%d: (%.3f, %.3f, %.3f, %.3f), %d, %d, %d, %p \n",is, p_shower->getPos().x(), p_shower->getPos().y(), p_shower->getPos().z(), p_shower->getEnergy(), p_shower->getBars().size(), p_shower->getNseeds(), p_shower->getCousinClusters().size(), p_shower);
  }
  cout<<"      ShowerV pos+E, Nbars, Nseed, cousin size and address: "<<endl;
  for(int is=0; is<p_tower->getCluster()[ib]->getShowerVCol().size(); is++){
    const Calo1DCluster* p_shower = p_tower->getCluster()[ib]->getShowerVCol()[is];
    printf("        #%d: (%.3f, %.3f, %.3f, %.3f), %d, %d, %d, %p \n",is, p_shower->getPos().x(), p_shower->getPos().y(), p_shower->getPos().z(), p_shower->getEnergy(), p_shower->getBars().size(), p_shower->getNseeds(), p_shower->getCousinClusters().size(), p_shower);
  }
}

printf("    LongiCluster size: (%d, %d) \n", p_tower->getLongiClusterUCol("ESLongiCluster").size(), p_tower->getLongiClusterVCol("ESLongiCluster").size());
for(int il=0; il<p_tower->getLongiClusterUCol("ESLongiCluster").size(); il++){
  const LongiCluster* p_Lclus = p_tower->getLongiClusterUCol("ESLongiCluster")[il];
  printf("      LClusterU #%d: hit size %d, cousin size %d \n", il, p_Lclus->getBarShowers().size(), p_Lclus->getCousinClusters().size() );
  for(int is=0; is<p_Lclus->getBarShowers().size(); is++){
    const Calo1DCluster* p_shower = p_Lclus->getBarShowers()[is];
    printf("        #%d: (%.3f, %.3f, %.3f, %.3f), %d, %d, %d, %p \n", is, p_shower->getPos().x(), p_shower->getPos().y(), p_shower->getPos().z(), p_shower->getEnergy(), p_shower->getBars().size(), p_shower->getNseeds(), p_shower->getCousinClusters().size(), p_shower);
  }

}
for(int il=0; il<p_tower->getLongiClusterVCol("ESLongiCluster").size(); il++){
  const LongiCluster* p_Lclus = p_tower->getLongiClusterVCol("ESLongiCluster")[il];
  printf("      LClusterV #%d: hit size %d, cousin size %d \n", il, p_Lclus->getBarShowers().size(), p_Lclus->getCousinClusters().size() );
  for(int is=0; is<p_Lclus->getBarShowers().size(); is++){
    const Calo1DCluster* p_shower = p_Lclus->getBarShowers()[is];
    printf("        #%d: (%.3f, %.3f, %.3f, %.3f), %d, %d, %d, %p \n", is, p_shower->getPos().x(), p_shower->getPos().y(), p_shower->getPos().z(), p_shower->getEnergy(), p_shower->getBars().size(), p_shower->getNseeds(), p_shower->getCousinClusters().size(), p_shower);
  }
}
}
*/


  }

  p_3DClusters = nullptr;
  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::LongiClusterLinking( std::vector<PandoraPlus::Calo2DCluster*>& m_2DClusCol, 
                                                    std::vector<const PandoraPlus::LongiCluster*>& m_oldClusCol, 
                                                    std::vector<const PandoraPlus::LongiCluster*>& m_outClusCol )
{
  if(m_2DClusCol.size()==0 || m_oldClusCol.size()==0) return StatusCode::SUCCESS;


  bool fl_isXclus = (m_oldClusCol[0]->getSlayer()==0);
  std::vector<PandoraPlus::LongiCluster*> m_cluscol; m_cluscol.clear();
  std::vector<const PandoraPlus::Calo1DCluster*> m_showers; m_showers.clear(); 
  std::vector<const PandoraPlus::Calo1DCluster*> m_showerswotrk; m_showerswotrk.clear(); 

  for(int ib=0; ib<m_2DClusCol.size(); ib++){
    std::vector<const PandoraPlus::Calo1DCluster*> tmp_showersinblock; tmp_showersinblock.clear(); 
    if(fl_isXclus) tmp_showersinblock = m_2DClusCol[ib]->getShowerUCol();
    else tmp_showersinblock = m_2DClusCol[ib]->getShowerVCol();

    int dlayer = m_2DClusCol[ib]->getDlayer();
    bool fl_hastrk = false; 
    for(int icl=0; icl<m_oldClusCol.size(); icl++) 
      if(dlayer<=m_oldClusCol[icl]->getEndDlayer() && dlayer>=m_oldClusCol[icl]->getBeginningDlayer()) { fl_hastrk=true; break; }

    if(fl_hastrk) m_showers.insert(m_showers.end(), tmp_showersinblock.begin(), tmp_showersinblock.end());
    else m_showerswotrk.insert(m_showerswotrk.end(), tmp_showersinblock.begin(), tmp_showersinblock.end());
  }


  //Update old Hough clusters
  for(int ic=0; ic<m_oldClusCol.size(); ic++){
    PandoraPlus::LongiCluster* m_newClus = new PandoraPlus::LongiCluster();

    for(int is=0; is<m_oldClusCol[ic]->getBarShowers().size(); is++){
      const PandoraPlus::Calo1DCluster* m_shower = m_oldClusCol[ic]->getBarShowers()[is];

      const PandoraPlus::Calo1DCluster* m_selshower = NULL;
      bool fl_foundshower = false; 
      for(int js=0; js<m_showers.size(); js++){
        bool fl_inTower = false;
        for(int itw=0; itw<m_showers[js]->getTowerID().size() && !fl_inTower; itw++){
        for(int jtw=0; jtw<m_shower->getTowerID().size(); jtw++){
          if(m_showers[js]->getTowerID()[itw]==m_shower->getTowerID()[jtw]) {fl_inTower=true; break;}
        }}


        if( fl_inTower &&
            m_showers[js]->getDlayer() ==m_shower->getDlayer() && 
            m_showers[js]->getSeeds().size()==1 && 
            (m_showers[js]->getSeeds()[0]->getPosition()-m_shower->getPos()).Mag()<10 ) 
        {m_selshower = m_showers[js]; fl_foundshower=true; break; }
      }
      if(fl_foundshower && m_selshower!=NULL) m_newClus->addBarShower( m_selshower, 0);
    }

    m_newClus->setHoughPars(m_oldClusCol[ic]->getHoughAlpha(), m_oldClusCol[ic]->getHoughRho() );
    m_newClus->setIntercept(m_oldClusCol[ic]->getHoughIntercept());
    //m_newClus->FitAxis();
    m_cluscol.push_back( m_newClus );
  }


  //Merge showers into closest Hough cluster
  std::sort( m_showerswotrk.begin(), m_showerswotrk.end(), compLayer );
  for(int is=0; is<m_showerswotrk.size(); is++) 
    MergeToClosestCluster( m_showerswotrk[is], m_cluscol );


/*
cout<<"    After merging shower: printout LongiClusters "<<endl;
for(int aa=0; aa<m_cluscol.size(); aa++){
printf("      LongiCluster pos: (%.2f, %.2f, %.2f) \n", m_cluscol[aa]->getPos().x(), m_cluscol[aa]->getPos().y(), m_cluscol[aa]->getPos().z() );
for(int ii=0; ii<m_cluscol[aa]->getBarShowers().size(); ii++)
printf("    #%d shower: pos/E=(%.2f, %.2f, %.2f, %.3f), cellID=(%d, %d, %d, %d, %d) \n", ii,
    m_cluscol[aa]->getBarShowers()[ii]->getPos().x(), m_cluscol[aa]->getBarShowers()[ii]->getPos().y(), m_cluscol[aa]->getBarShowers()[ii]->getPos().z(),
    m_cluscol[aa]->getBarShowers()[ii]->getEnergy(),
    m_cluscol[aa]->getBarShowers()[ii]->getModule(), m_cluscol[aa]->getBarShowers()[ii]->getStave(), m_cluscol[aa]->getBarShowers()[ii]->getPart(),
    m_cluscol[aa]->getBarShowers()[ii]->getDlayer(), m_cluscol[aa]->getBarShowers()[ii]->getSlayer() );

}
*/

  //convert to const
  for(int icl=0; icl<m_cluscol.size(); icl++) m_outClusCol.push_back(m_cluscol[icl]);

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::Split3DClusterToTowers( PandoraPlus::Calo3DCluster* m_3dcluster ){
  
  std::vector<PandoraPlus::Calo3DCluster*> m_towerCol; m_towerCol.clear(); 

  std::vector< std::vector<int> > m_towerID = m_3dcluster->getTowerID();
  std::vector<const PandoraPlus::LongiCluster*>  m_LongiClusterUCol = m_3dcluster->getLongiClusterUCol("tmpCluster");
  std::vector<const PandoraPlus::LongiCluster*>  m_LongiClusterVCol = m_3dcluster->getLongiClusterVCol("tmpCluster");
/*
cout<<"Split3DClusterToTowers: Input Info. "<<endl;
cout<<"Split3DClusterToTowers: 2DCluster size: "<<m_3dcluster->getCluster().size()<<endl;
for(int ic=0; ic<m_3dcluster->getCluster().size(); ic++){
  printf("    Check 2DCluster %d: cover tower %d: ", ic, m_3dcluster->getCluster()[ic]->getTowerID().size());
  for(int it=0; it< m_3dcluster->getCluster()[ic]->getTowerID().size(); it++)
    printf("[%d, %d, %d], ", m_3dcluster->getCluster()[ic]->getTowerID()[it][0], m_3dcluster->getCluster()[ic]->getTowerID()[it][1], m_3dcluster->getCluster()[ic]->getTowerID()[it][2]);
  cout<<endl;
  printf("    BarShower size: (%d, %d) \n", m_3dcluster->getCluster()[ic]->getShowerUCol().size(), m_3dcluster->getCluster()[ic]->getShowerVCol().size());
  cout<<"    BarShowerU: ";
  for(int is=0; is<m_3dcluster->getCluster()[ic]->getShowerUCol().size(); is++)
    printf(" [%d, %d, %p], ", m_3dcluster->getCluster()[ic]->getShowerUCol()[is]->getDlayer(), 
                              m_3dcluster->getCluster()[ic]->getShowerUCol()[is]->getTowerID().size(), 
                              m_3dcluster->getCluster()[ic]->getShowerUCol()[is] );
  cout<<endl;
  cout<<"    BarShowerV: ";
  for(int is=0; is<m_3dcluster->getCluster()[ic]->getShowerVCol().size(); is++)
    printf(" [%d, %d, %p], ", m_3dcluster->getCluster()[ic]->getShowerVCol()[is]->getDlayer(),
                              m_3dcluster->getCluster()[ic]->getShowerVCol()[is]->getTowerID().size(),
                              m_3dcluster->getCluster()[ic]->getShowerVCol()[is] );
  cout<<endl;
}
*/

cout<<"Split3DClusterToTowers: Get LongiCluster U: "<<m_LongiClusterUCol.size()<<", V: "<<m_LongiClusterVCol.size()<<endl;
for(int il=0; il<m_LongiClusterUCol.size(); il++){
  printf("  Check LongiClusterU #%d: shower size %d \n", il, m_LongiClusterUCol[il]->getBarShowers().size());
  cout<<"    BarShower: "<<endl;
  for(int is=0; is<m_LongiClusterUCol[il]->getBarShowers().size(); is++)
    printf(" [%d, %d (%d, %d, %d, %d), %p], \n", m_LongiClusterUCol[il]->getBarShowers()[is]->getDlayer(),
                              m_LongiClusterUCol[il]->getBarShowers()[is]->getTowerID().size(),
                              m_LongiClusterUCol[il]->getBarShowers()[is]->getTowerID()[0][0],
                              m_LongiClusterUCol[il]->getBarShowers()[is]->getTowerID()[0][1],
                              m_LongiClusterUCol[il]->getBarShowers()[is]->getTowerID()[0][2],
                              m_LongiClusterUCol[il]->getBarShowers()[is]->getDlayer(),
                              m_LongiClusterUCol[il]->getBarShowers()[is] );
  cout<<endl;
}

for(int il=0; il<m_LongiClusterVCol.size(); il++){
  printf("  Check LongiClusterV #%d: shower size %d \n ", il, m_LongiClusterVCol[il]->getBarShowers().size());
  cout<<"    BarShower: "<<endl;
  for(int is=0; is<m_LongiClusterVCol[il]->getBarShowers().size(); is++)
    printf(" [%d, %d (%d, %d, %d, %d), %p], \n", m_LongiClusterVCol[il]->getBarShowers()[is]->getDlayer(),
                              m_LongiClusterVCol[il]->getBarShowers()[is]->getTowerID().size(),
                              m_LongiClusterVCol[il]->getBarShowers()[is]->getTowerID()[0][0],
                              m_LongiClusterVCol[il]->getBarShowers()[is]->getTowerID()[0][1],
                              m_LongiClusterVCol[il]->getBarShowers()[is]->getTowerID()[0][2],
                              m_LongiClusterVCol[il]->getBarShowers()[is]->getDlayer(),
                              m_LongiClusterVCol[il]->getBarShowers()[is] );
  cout<<endl;
}



  std::map<std::vector<int>, std::vector<const PandoraPlus::Calo2DCluster*> > map_2DCluster;
  std::map<std::vector<int>, std::vector<const PandoraPlus::LongiCluster*> > map_LongiClusterU; 
  std::map<std::vector<int>, std::vector<const PandoraPlus::LongiCluster*> > map_LongiClusterV; 

//cout<<"  Start Split 2DClusters"<<endl;

  //Split 2DCluster
  for(int ic=0; ic<m_3dcluster->getCluster().size(); ic++){
    std::vector< std::vector<int> > cl_towerID = m_3dcluster->getCluster()[ic]->getTowerID();
//cout<<"    2DCluster #"<<ic<<": cover "<<cl_towerID.size()<<" towers: "<<endl;
//for(int a=0; a<cl_towerID.size(); a++) printf("    [%d, %d, %d] \t", cl_towerID[a][0], cl_towerID[a][1], cl_towerID[a][2]);
//cout<<endl;

    //2DCluster does not cover tower
    if(cl_towerID.size()==1){
      map_2DCluster[cl_towerID[0]].push_back(m_3dcluster->getCluster()[ic]);
//cout<<"    2DCluster #"<<ic<<" is in one tower. "<<endl;
      continue; 
    }

    //2DCluster covers towers:
    std::map<std::vector<int>, std::vector<const CaloUnit*> > barUMap;  //slayer == 0.
    std::map<std::vector<int>, std::vector<const CaloUnit*> > barVMap;  //slayer == 1.
    std::map<std::vector<int>, std::vector<const Calo1DCluster*> > barShowerUMap;
    std::map<std::vector<int>, std::vector<const Calo1DCluster*> > barShowerVMap;
    std::map<std::vector<int>, std::vector<const Calo1DCluster*> > barClusterUMap;
    std::map<std::vector<int>, std::vector<const Calo1DCluster*> > barClusterVMap;

//printf("    2DCluster #%d bar size: (%d, %d). BarShower size: (%d, %d) \n", ic, m_3dcluster->getCluster()[ic]->getBarUCol().size(),
//                                                                                m_3dcluster->getCluster()[ic]->getBarVCol().size(),
//                                                                                m_3dcluster->getCluster()[ic]->getShowerUCol().size(),
//                                                                                m_3dcluster->getCluster()[ic]->getShowerVCol().size()  );

    //Bars: 
    for(int ib=0; ib<m_3dcluster->getCluster()[ic]->getBarUCol().size(); ib++){
      const CaloUnit* p_bar = m_3dcluster->getCluster()[ic]->getBarUCol()[ib];
      if(!p_bar) continue; 

      std::vector<int> towerid(3); 
      towerid[0] = p_bar->getModule(); 
      towerid[1] = p_bar->getPart(); 
      towerid[2] = p_bar->getStave(); 

      barUMap[towerid].push_back(p_bar);
      p_bar = nullptr; 
    }
//cout<<"    Check BarU map: size = "<<barUMap.size()<<", Bar size in each tower: "<<endl;
//for(auto iter:barUMap) printf("    In tower [%d, %d, %d]: %d \n",iter.first[0],iter.first[1],iter.first[2],iter.second.size() );

    for(int ib=0; ib<m_3dcluster->getCluster()[ic]->getBarVCol().size(); ib++){
      const CaloUnit* p_bar = m_3dcluster->getCluster()[ic]->getBarVCol()[ib];
      if(!p_bar) continue;

      std::vector<int> towerid(3);
      towerid[0] = p_bar->getModule();
      towerid[1] = p_bar->getPart();
      towerid[2] = p_bar->getStave();

      barVMap[towerid].push_back(p_bar);
      p_bar = nullptr;
    }
//cout<<"    Check BarV map: size = "<<barVMap.size()<<", Bar size in each tower: "<<endl;
//for(auto iter:barVMap) printf("    In tower [%d, %d, %d]: %d \n",iter.first[0],iter.first[1],iter.first[2],iter.second.size() );
//cout<<endl;

    //BarShowers: 
    for(int is=0; is<m_3dcluster->getCluster()[ic]->getShowerUCol().size(); is++){
      Calo1DCluster* p_shower = const_cast<PandoraPlus::Calo1DCluster*>(m_3dcluster->getCluster()[ic]->getShowerUCol()[is]);
      std::vector< std::vector<int> > towerids = p_shower->getTowerID();

//cout<<"    Split BarShower #"<<is<<": Layer "<<p_shower->getDlayer()<<", it covers "<<towerids.size()<<" towers. address "<<p_shower<<endl;

      if(towerids.size()==1){
        barShowerUMap[towerids[0]].push_back(p_shower);
        continue; 
      }

      //Split BarShowers into several:
      std::map< std::vector<int>, Calo1DCluster* > tmp_showerMap; tmp_showerMap.clear();
      for(int ib=0; ib<p_shower->getBars().size(); ib++){
        std::vector<int> m_id(3); 
        m_id[0] = p_shower->getBars()[ib]->getModule();
        m_id[1] = p_shower->getBars()[ib]->getPart();
        m_id[2] = p_shower->getBars()[ib]->getStave();

        if(tmp_showerMap.find(m_id)==tmp_showerMap.end()){
          Calo1DCluster* tmp_shower = new Calo1DCluster();
          tmp_shower->addUnit(p_shower->getBars()[ib]);
          tmp_showerMap[m_id] = tmp_shower;
        }
        else tmp_showerMap[m_id]->addUnit(p_shower->getBars()[ib]);
      }

      for(auto &iter:tmp_showerMap){
        iter.second->setSeed();
        for(auto &iter1:tmp_showerMap){
          if(iter==iter1) continue;
          iter.second->addCousinCluster( iter1.second );
        }
        p_shower->addChildCluster( iter.second );
      }
      
      for(auto &iter:tmp_showerMap) barShowerUMap[iter.first].push_back( iter.second );
      p_shower = nullptr;
    }
//cout<<"    Check ShowerU map: size = "<<barShowerUMap.size()<<", BarShower size in each tower: "<<endl;
//for(auto iter:barShowerUMap){
//  printf("    In tower [%d, %d, %d]: %d :",iter.first[0],iter.first[1],iter.first[2],iter.second.size() );
//  for(int as=0; as<iter.second.size(); as++){
//    printf(" [%d, %d, %d, %p], ", iter.second[as]->getDlayer(), iter.second[as]->getBars().size(), iter.second[as]->getNseeds(), iter.second[as] );
//  }
//  cout<<endl;
//}

    for(int is=0; is<m_3dcluster->getCluster()[ic]->getShowerVCol().size(); is++){
      Calo1DCluster* p_shower = const_cast<PandoraPlus::Calo1DCluster*>(m_3dcluster->getCluster()[ic]->getShowerVCol()[is]);
      std::vector< std::vector<int> > towerids = p_shower->getTowerID();

      if(towerids.size()==1){
        barShowerVMap[towerids[0]].push_back(p_shower);
        continue;
      }

      //Split BarShowers into several:
      std::map< std::vector<int>, Calo1DCluster* > tmp_showerMap; tmp_showerMap.clear();
      for(int ib=0; ib<p_shower->getBars().size(); ib++){
        std::vector<int> m_id(3);
        m_id[0] = p_shower->getBars()[ib]->getModule();
        m_id[1] = p_shower->getBars()[ib]->getPart();
        m_id[2] = p_shower->getBars()[ib]->getStave();

        if(tmp_showerMap.find(m_id)==tmp_showerMap.end()){
          Calo1DCluster* tmp_shower = new Calo1DCluster();
          tmp_shower->addUnit(p_shower->getBars()[ib]);
          tmp_showerMap[m_id] = tmp_shower;
        }
        else tmp_showerMap[m_id]->addUnit(p_shower->getBars()[ib]);
      }

      for(auto &iter: tmp_showerMap){
        for(auto &iter1: tmp_showerMap){
          if(iter==iter1) continue;
          iter.second->addCousinCluster( iter1.second );
        }
        iter.second->setSeed();
        p_shower->addChildCluster( iter.second );
      }

      for(auto &iter: tmp_showerMap) barShowerVMap[iter.first].push_back( iter.second );
      p_shower=nullptr; 
    }
/*
cout<<"    Check ShowerV map: size = "<<barShowerVMap.size()<<", BarShower size in each tower: "<<endl;
for(auto iter:barShowerVMap) {
  printf("    In tower [%d, %d, %d]: %d : ",iter.first[0],iter.first[1],iter.first[2],iter.second.size() );
  for(int as=0; as<iter.second.size(); as++){
    printf(" [%d, %d, %d, %p], ", iter.second[as]->getDlayer(), iter.second[as]->getBars().size(), iter.second[as]->getNseeds(), iter.second[as] );
  }
  cout<<endl;
}
*/
    for(int it=0; it<cl_towerID.size(); it++){
      PandoraPlus::Calo2DCluster* m_2dclus = new PandoraPlus::Calo2DCluster();

      m_2dclus->addTowerID( cl_towerID[it] );
      m_2dclus->setBarUCol( barUMap[cl_towerID[it]] );
      m_2dclus->setBarVCol( barVMap[cl_towerID[it]] );
      m_2dclus->setShowerUCol( barShowerUMap[cl_towerID[it]] );
      m_2dclus->setShowerVCol( barShowerVMap[cl_towerID[it]] );

      map_2DCluster[cl_towerID[it]].push_back( m_2dclus );
    }

  }

cout<<"  Start Split LongiClusterU"<<endl;
  //Split LongiClusterU
  for(int il=0; il<m_LongiClusterUCol.size(); il++){
    if(m_LongiClusterUCol[il]->getBarShowers().size()==0) {std::cout<<"WARNING: Have an empty LongiCluster! Skip it! "<<std::endl; continue;}

//printf("    LCluster #%d: shower size %d. ", il, m_LongiClusterUCol[il]->getBarShowers().size());

    bool fl_coverTower = false; 
    std::vector< std::vector<int> > tmp_towerIDCol; tmp_towerIDCol.clear();
    for(int ib=0; ib<m_LongiClusterUCol[il]->getBarShowers().size(); ib++){
      std::vector< std::vector<int> > tmp_id = m_LongiClusterUCol[il]->getBarShowers()[ib]->getTowerID();
      int m_ntower = tmp_id.size();
      if( m_ntower>1 )  fl_coverTower=true; 
      tmp_towerIDCol.insert( tmp_towerIDCol.end(), tmp_id.begin(), tmp_id.end() );
    }
    std::sort( tmp_towerIDCol.begin(), tmp_towerIDCol.end() );
    auto iter_id = std::unique( tmp_towerIDCol.begin(), tmp_towerIDCol.end() );
    tmp_towerIDCol.erase( iter_id, tmp_towerIDCol.end() );
    if(tmp_towerIDCol.size()>1) fl_coverTower=true;

cout<<" Cover tower: "<<fl_coverTower<<", tower size "<<tmp_towerIDCol.size()<<endl;

    //LongiCluster does not cover tower: 
    if(!fl_coverTower){
      std::vector<int> cl_towerID =  m_LongiClusterUCol[0]->getBarShowers()[0]->getTowerID()[0]; 
      map_LongiClusterU[cl_towerID].push_back(m_LongiClusterUCol[il]);
      continue; 
    }

    //LongiCluster covers towers:
    std::map<std::vector<int>, PandoraPlus::LongiCluster* > tmp_LongiClusMaps; tmp_LongiClusMaps.clear();
    for(int ib=0; ib<m_LongiClusterUCol[il]->getBarShowers().size(); ib++){
      const PandoraPlus::Calo1DCluster* p_shower = m_LongiClusterUCol[il]->getBarShowers()[ib];

//printf("      Shower %d: address %p, Children size %d, tower size %d \n ", ib, p_shower, p_shower->getChildClusters().size(), p_shower->getTowerID().size() );

      if(p_shower->getChildClusters().size()==0 && p_shower->getTowerID().size()==1){
        std::vector<int> b_towerID = p_shower->getTowerID()[0];
        if( tmp_LongiClusMaps.find(b_towerID)==tmp_LongiClusMaps.end() ){
          PandoraPlus::LongiCluster* tmp_clus = new PandoraPlus::LongiCluster();
          tmp_clus->addBarShower( p_shower, 0);
          tmp_LongiClusMaps[b_towerID] = tmp_clus;
        }
        else tmp_LongiClusMaps[b_towerID]->addBarShower( p_shower, 0);
      }

      else{
        for(int ic=0; ic<p_shower->getChildClusters().size(); ic++){    //Loop for ChildClusters
          if( p_shower->getChildClusters()[ic]->getTowerID().size()!=1 ) { std::cout<<"WARNING: ChildCluster still cover towers! "<<std::endl; continue; }
          std::vector<int> b_towerID = p_shower->getChildClusters()[ic]->getTowerID()[0];

          if( tmp_LongiClusMaps.find(b_towerID)==tmp_LongiClusMaps.end() ){
            PandoraPlus::LongiCluster* tmp_clus = new PandoraPlus::LongiCluster();
            tmp_clus->addBarShower( p_shower->getChildClusters()[ic], 0);
            tmp_LongiClusMaps[b_towerID] = tmp_clus;
          }
          else tmp_LongiClusMaps[b_towerID]->addBarShower( p_shower->getChildClusters()[ic], 0);
        }
      }

      p_shower = nullptr;
    }

    //Connect cousins
    if(tmp_LongiClusMaps.size()>1){
      std::map<std::vector<int>, PandoraPlus::LongiCluster* >::iterator iter = tmp_LongiClusMaps.begin();
      for(iter; iter!=tmp_LongiClusMaps.end(); iter++){
        std::map<std::vector<int>, PandoraPlus::LongiCluster* >::iterator iter1 = tmp_LongiClusMaps.begin();
        for(iter1; iter1!=tmp_LongiClusMaps.end(); iter1++)
          if(iter!=iter1) iter->second->addCousinCluster(iter1->second);
      }
    }


    for(auto &iter : tmp_LongiClusMaps) map_LongiClusterU[iter.first].push_back(iter.second);
  }


cout<<"  Start Split LongiClusterV"<<endl;
  //Split LongiClusterV
  for(int il=0; il<m_LongiClusterVCol.size(); il++){
    if(m_LongiClusterVCol[il]->getBarShowers().size()==0) {std::cout<<"WARNING: Have an empty LongiCluster! Skip it! "<<std::endl; continue;}

    //Check if the LongiCluster covers towers. 
    bool fl_coverTower = false;
    std::vector< std::vector<int> > tmp_towerIDCol; tmp_towerIDCol.clear();
    for(int ib=0; ib<m_LongiClusterVCol[il]->getBarShowers().size(); ib++){
      std::vector< std::vector<int> > tmp_id = m_LongiClusterVCol[il]->getBarShowers()[ib]->getTowerID();
      int m_ntower = tmp_id.size();
      if( m_ntower>1 )  fl_coverTower=true;
      tmp_towerIDCol.insert( tmp_towerIDCol.end(), tmp_id.begin(), tmp_id.end() );
    }
    std::sort( tmp_towerIDCol.begin(), tmp_towerIDCol.end() );
    auto iter_id = std::unique( tmp_towerIDCol.begin(), tmp_towerIDCol.end() );
    tmp_towerIDCol.erase( iter_id, tmp_towerIDCol.end() );
    if(tmp_towerIDCol.size()>1) fl_coverTower=true;

cout<<" Cover tower: "<<fl_coverTower<<", tower size "<<tmp_towerIDCol.size()<<endl;

    //LongiCluster does not cover tower:
    if(!fl_coverTower){
      std::vector<int> cl_towerID =  m_LongiClusterVCol[0]->getBarShowers()[0]->getTowerID()[0];
      map_LongiClusterV[cl_towerID].push_back(m_LongiClusterVCol[il]);
      continue;
    }

    //LongiCluster covers towers:
    std::map<std::vector<int>, PandoraPlus::LongiCluster* > tmp_LongiClusMaps; tmp_LongiClusMaps.clear();
    for(int ib=0; ib<m_LongiClusterVCol[il]->getBarShowers().size(); ib++){
      const PandoraPlus::Calo1DCluster* p_shower = m_LongiClusterVCol[il]->getBarShowers()[ib];

      if(p_shower->getChildClusters().size()==0 && p_shower->getTowerID().size()==1){
        std::vector<int> b_towerID = p_shower->getTowerID()[0];
        if( tmp_LongiClusMaps.find(b_towerID)==tmp_LongiClusMaps.end() ){
          PandoraPlus::LongiCluster* tmp_clus = new PandoraPlus::LongiCluster();
          tmp_clus->addBarShower( p_shower, 0);
          tmp_LongiClusMaps[b_towerID] = tmp_clus;
        }
        else tmp_LongiClusMaps[b_towerID]->addBarShower( p_shower, 0);
      }

      else{
        for(int ic=0; ic<p_shower->getChildClusters().size(); ic++){    //Loop for ChildClusters
          if( p_shower->getChildClusters()[ic]->getTowerID().size()!=1 ) { std::cout<<"WARNING: ChildCluster still cover towers! "<<std::endl; continue; }
          std::vector<int> b_towerID = p_shower->getChildClusters()[ic]->getTowerID()[0];

          if( tmp_LongiClusMaps.find(b_towerID)==tmp_LongiClusMaps.end() ){
            PandoraPlus::LongiCluster* tmp_clus = new PandoraPlus::LongiCluster();
            tmp_clus->addBarShower( p_shower->getChildClusters()[ic], 0);
            tmp_LongiClusMaps[b_towerID] = tmp_clus;
          }
          else tmp_LongiClusMaps[b_towerID]->addBarShower( p_shower->getChildClusters()[ic], 0);
        }
      }

      p_shower = nullptr;
    }

    if(tmp_LongiClusMaps.size()>1){
      std::map<std::vector<int>, PandoraPlus::LongiCluster* >::iterator iter = tmp_LongiClusMaps.begin();
      for(iter; iter!=tmp_LongiClusMaps.end(); iter++){
        std::map<std::vector<int>, PandoraPlus::LongiCluster* >::iterator iter1 = tmp_LongiClusMaps.begin();
        for(iter1; iter1!=tmp_LongiClusMaps.end(); iter1++) 
          if(iter!=iter1) iter->second->addCousinCluster(iter1->second);
      }
    }

    for(auto &iter : tmp_LongiClusMaps) map_LongiClusterV[iter.first].push_back(iter.second);
  }


//cout<<"  Form a new tower "<<endl;

  //Form a new tower:
  for(int it=0; it<m_towerID.size(); it++){
    PandoraPlus::Calo3DCluster* m_tower = new PandoraPlus::Calo3DCluster();
    m_tower->addTowerID( m_towerID[it] );
    for(int i2d=0; i2d<map_2DCluster[m_towerID[it]].size(); i2d++) m_tower->addUnit(map_2DCluster[m_towerID[it]][i2d]);
    m_tower->setLongiClusters( settings.map_stringPars["OutputLongiClusName"], map_LongiClusterU[m_towerID[it]], 
                               settings.map_stringPars["OutputLongiClusName"], map_LongiClusterV[m_towerID[it]] );
    m_towerCol.push_back(m_tower);
  }

//cout<<"  After splitting: tower size "<<m_towerCol.size()<<endl;

  for(int it=0; it<m_towerCol.size(); it++) m_3dcluster->addTower( m_towerCol[it] );
  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::Clustering( std::vector<const PandoraPlus::CaloUnit*>& barCol, 
                                           std::vector<PandoraPlus::Calo1DCluster*>& outClus,  
                                           std::vector<const PandoraPlus::LongiCluster*>& m_longiClusCol )
{
  if(barCol.size()==0) return StatusCode::SUCCESS;

//double totE = 0;
//for(int i=0; i<barCol.size(); i++) totE += barCol[i]->getEnergy();
//cout<<"Clustering: Raw bar energy "<<totE<<endl;

  int slayer = barCol[0]->getSlayer();
  //Neighbor clustering
  std::sort(barCol.begin(), barCol.end());
  std::vector<PandoraPlus::Calo1DCluster*> m_clusCol; m_clusCol.clear();
  for(int i=0;i<barCol.size();i++){
    if( (m_clusCol.size()!=0) && (m_clusCol[m_clusCol.size()-1]->isNeighbor(barCol[i])) &&
        ( (slayer==0 && m_clusCol[m_clusCol.size()-1]->getTowerID()[0][1]==barCol[i]->getPart()) ||
          (slayer==1 && m_clusCol[m_clusCol.size()-1]->getTowerID()[0][2]==barCol[i]->getStave()) ) ){

      m_clusCol[m_clusCol.size()-1]->addUnit( barCol[i] );
      continue;
    }

    PandoraPlus::Calo1DCluster* clus = new PandoraPlus::Calo1DCluster();
    clus->addUnit(barCol[i]);
    m_clusCol.push_back(clus);
  }

//totE = 0;
//cout<<"Clustering: After Neighbor clustering. Cluster size: "<<m_clusCol.size()<<endl;
//cout<<"  Bar size and energy in cluster: ";
//for(int j=0; j<m_clusCol.size(); j++){ printf("%d (%.3f), \t",m_clusCol[j]->getBars().size(), m_clusCol[j]->getEnergy() ); totE += m_clusCol[j]->getEnergy(); }
//cout<<endl;
//cout<<"total E: "<<totE<<endl;
//cout<<endl;

  //Save out CaloBars in longiClusCol as seed. 
  std::vector<const PandoraPlus::CaloUnit*> m_seedbars; m_seedbars.clear(); 
  for(int il=0; il<m_longiClusCol.size(); il++){
  for(int is=0; is<m_longiClusCol[il]->getBarShowers().size(); is++){
    std::vector<const PandoraPlus::CaloUnit*> m_bars = m_longiClusCol[il]->getBarShowers()[is]->getBars();
    m_seedbars.insert(m_seedbars.end(), m_bars.begin(), m_bars.end());
  }}

//cout<<"Clustering: seed size from LongiCluster: "<<m_seedbars.size()<<endl;

  //Find seed with LongiCluster
  for(int ic=0; ic<m_clusCol.size(); ic++){
  for(int ib=0; ib<m_clusCol[ic]->getBars().size(); ib++){
      std::vector<const PandoraPlus::CaloUnit*>::iterator iter = find(m_seedbars.begin(), m_seedbars.end(), m_clusCol[ic]->getBars()[ib]);
      if(iter==m_seedbars.end()) continue; 
      m_clusCol[ic]->addSeed( *iter );
    
  }}

/*
  //Cross-check 0-seed clusters: if E > Eth then find seed with localMax. 
  for(int i=0;i<m_clusCol.size();i++){
    if(m_clusCol[i]->getSeeds().size()>0) continue;
    if(m_clusCol[i]->getEnergy()<settings.map_floatPars["Eth_Shower"]) continue; 
    std::vector<const PandoraPlus::CaloUnit*> m_seedVec; m_seedVec.clear();
    findSeeds(m_clusCol[i], m_seedVec);

    m_clusCol[i]->setSeeds(m_seedVec);
    m_clusCol[i]->getScndMoment();
  }
*/

/*
cout<<"Clustering: After seed finding. Cluster size: "<<m_clusCol.size()<<endl;
cout<<"  Bar size in cluster: ";
for(int j=0; j<m_clusCol.size(); j++) cout<<m_clusCol[j]->getBars().size()<<'\t';
cout<<endl;
cout<<"  Seed size in cluster: ";
for(int j=0; j<m_clusCol.size(); j++) cout<<m_clusCol[j]->getSeeds().size()<<'\t';
cout<<endl;
cout<<"    Seed position: "<<endl;
for(int j=0; j<m_clusCol.size(); j++){
  for(int a=0; a<m_clusCol[j]->getSeeds().size(); a++) printf("    (%.3f, %.3f, %.3f ), ", m_clusCol[j]->getSeeds()[a]->getPosition().x(),  m_clusCol[j]->getSeeds()[a]->getPosition().y(),  m_clusCol[j]->getSeeds()[a]->getPosition().z() );
  cout<<endl;
}
*/

  //Merge clusters without seed
  for(int ic=0; ic<m_clusCol.size(); ic++){
    if(m_clusCol[ic]->getNseeds()==0){
      MergeToClosestCluster( m_clusCol[ic], m_clusCol );
      delete m_clusCol[ic]; m_clusCol[ic]=NULL;
      m_clusCol.erase(m_clusCol.begin()+ic);
      ic--;
    }
  }

//totE = 0;
//cout<<"Clustering: Final cluster energy"<<endl;
//for(int j=0; j<m_clusCol.size(); j++){ printf("%d (%.3f), \t",m_clusCol[j]->getBars().size(), m_clusCol[j]->getEnergy() ); totE += m_clusCol[j]->getEnergy(); }
//cout<<endl;

  outClus = m_clusCol; 

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::Clustering( std::vector<const PandoraPlus::CaloUnit*>& barCol,
                                           std::vector<PandoraPlus::Calo1DCluster*>& outClus){
  if(barCol.size()==0) return StatusCode::SUCCESS;

  int slayer = barCol[0]->getSlayer();
  //Neighbor clustering
  std::sort(barCol.begin(), barCol.end());
  std::vector<PandoraPlus::Calo1DCluster*> m_clusCol; m_clusCol.clear();
  for(int i=0;i<barCol.size();i++){
    if( (m_clusCol.size()!=0) && (m_clusCol[m_clusCol.size()-1]->isNeighbor(barCol[i])) && 
        ( (slayer==0 && m_clusCol[m_clusCol.size()-1]->getTowerID()[0][1]==barCol[i]->getPart()) || 
          (slayer==1 && m_clusCol[m_clusCol.size()-1]->getTowerID()[0][2]==barCol[i]->getStave()) ) ){

      m_clusCol[m_clusCol.size()-1]->addUnit( barCol[i] );
      continue;
    }

    PandoraPlus::Calo1DCluster* clus = new PandoraPlus::Calo1DCluster();
    clus->addUnit(barCol[i]);
    m_clusCol.push_back(clus);
  }

  //Find seed in cluster
  for(int i=0;i<m_clusCol.size();i++){
    std::vector<const PandoraPlus::CaloUnit*> m_seedVec; m_seedVec.clear();
    findSeeds(m_clusCol[i], m_seedVec);

    m_clusCol[i]->setSeeds(m_seedVec);
    m_clusCol[i]->getScndMoment();
  }

  //Merge clusters without seed
  for(int ic=0; ic<m_clusCol.size(); ic++){
    if(m_clusCol[ic]->getNseeds()==0){
      MergeToClosestCluster( m_clusCol[ic], m_clusCol );
      delete m_clusCol[ic]; m_clusCol[ic]=NULL;
      m_clusCol.erase(m_clusCol.begin()+ic);
      ic--;
    }
  }
  outClus = m_clusCol;


  return StatusCode::SUCCESS;
}



StatusCode EnergySplittingAlg::ClusterSplitting( PandoraPlus::Calo1DCluster* m_cluster, std::vector<const PandoraPlus::Calo1DCluster*>& outshCol ){
//cout<<"ClusterSplitting: input cluster seed size = "<<m_cluster->getSeeds().size()<<endl;
//cout<<"Seed position: ";
//for(int a=0; a<m_cluster->getNseeds(); a++) printf(" (%.2f, %.2f, %.2f) \t", m_cluster->getSeeds()[a]->getPosition().x(), 
//                                                                             m_cluster->getSeeds()[a]->getPosition().y(),
//                                                                             m_cluster->getSeeds()[a]->getPosition().z() );
//cout<<endl;

  //No seed in cluster: return empty vector.
  if(m_cluster->getNseeds()==0) { std::cout<<"WARNING: Still have no-seed cluster!!"<<std::endl; return StatusCode::SUCCESS; }

  //1 seed or second moment less than threshold: Not split. Turn cluster to shower and return
  else if(m_cluster->getNseeds()<2 || m_cluster->getScndMoment()<settings.map_floatPars["th_split"]){
    PandoraPlus::Calo1DCluster* shower = new PandoraPlus::Calo1DCluster();
    shower->setBars(m_cluster->getBars());
    if(m_cluster->getNseeds()!=0) shower->addSeed(m_cluster->getSeeds()[0]);
    else shower->addSeed(nullptr);
    shower->setIDInfo();
    outshCol.push_back(shower);
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
                                                      std::vector<PandoraPlus::LongiCluster*>& m_clusters ){
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
      m_distance = (sh_pos-m_clusters[ic]->getBarShowersInLayer(m_clusters[ic]->getBeginningDlayer())[0]->getPos()).Mag();     
    else if(dlayer>=maxLayer)
      m_distance = (sh_pos-m_clusters[ic]->getBarShowersInLayer(m_clusters[ic]->getEndDlayer())[0]->getPos()).Mag();    
    else
      m_distance = (sh_pos-m_clusters[ic]->getPos()).Mag();

//cout<<"    Cluster #"<<ic<<": distance with shower = "<<m_distance<<endl;
    
    if( m_distance<minR )  { minR=m_distance; index_cluster=ic; }
  }

//printf("  minR = %.3f, in #cl %d \n", minR, index_cluster);

  if(index_cluster>=0) m_clusters[index_cluster]->addBarShower( m_shower, 1 );
  //delete m_shower; m_shower = NULL;

  return StatusCode::SUCCESS;
}

StatusCode EnergySplittingAlg::findSeeds( PandoraPlus::Calo1DCluster* m_cluster, std::vector<const PandoraPlus::CaloUnit*>& seedCol){

  for(int i=0;i<m_cluster->getBars().size();i++){
    const PandoraPlus::CaloUnit* ibar = m_cluster->getBars()[i];
    std::vector<const PandoraPlus::CaloUnit*> m_neighbor = getNeighbors(m_cluster, ibar);
    if(m_neighbor.size()==0 && ibar->getEnergy()>settings.map_floatPars["Eth_Seed"])
      seedCol.push_back(ibar); 

    else {
      bool isLocalMax=true; 
      for(int j=0;j<m_neighbor.size();j++){
        if(m_neighbor[j]->getEnergy()>ibar->getEnergy()) isLocalMax=false;
      }
      if(ibar->getEnergy()>settings.map_floatPars["Eth_Seed"] && isLocalMax) seedCol.push_back(ibar);
    }
  }

  return StatusCode::SUCCESS;
}

std::vector<const PandoraPlus::CaloUnit*>  EnergySplittingAlg::getNeighbors(PandoraPlus::Calo1DCluster* m_cluster, const PandoraPlus::CaloUnit* seed){

  std::vector<const PandoraPlus::CaloUnit*> m_neighbor; m_neighbor.clear();
  std::vector<const PandoraPlus::CaloUnit*> barCol = m_cluster->getBars();
  for(int i=0;i<barCol.size();i++){
    bool fl_neighbor = false;
    if( seed->getModule()==barCol[i]->getModule() &&
        seed->getPart()==barCol[i]->getPart() &&
        seed->getStave()==barCol[i]->getStave() &&
        seed->getDlayer()==barCol[i]->getDlayer() &&
        seed->getSlayer()==barCol[i]->getSlayer() &&
        abs( seed->getBar()-barCol[i]->getBar() )==1 ) fl_neighbor=true;
    else if( seed->getModule()==barCol[i]->getModule() &&
             seed->getStave()==barCol[i]->getStave() &&
             ( ( seed->getPart()-barCol[i]->getPart()==1 && seed->isAtLowerEdgePhi() && barCol[i]->isAtUpperEdgePhi() ) ||
               ( barCol[i]->getPart()-seed->getPart()==1 && seed->isAtUpperEdgePhi() && barCol[i]->isAtLowerEdgePhi() ) ) ) fl_neighbor=true;
    else if( seed->getModule()==barCol[i]->getModule() &&
             seed->getPart()==barCol[i]->getPart() &&
             ( ( seed->getStave()-barCol[i]->getStave()==1 && seed->isAtLowerEdgeZ() && barCol[i]->isAtUpperEdgeZ() ) ||
               ( barCol[i]->getStave()-seed->getStave()==1 && seed->isAtUpperEdgeZ() && barCol[i]->isAtLowerEdgeZ() ) ) ) fl_neighbor=true;


    if(fl_neighbor) m_neighbor.push_back( barCol[i] );

  }
  if(m_neighbor.size()>2) std::cout<<"WARNING: more than 2 hits in neighborCol!!"<<std::endl;

  return m_neighbor;
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
