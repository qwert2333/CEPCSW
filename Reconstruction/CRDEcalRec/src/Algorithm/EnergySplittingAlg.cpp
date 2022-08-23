#ifndef _ENERGYSPLITTING_ALG_C
#define _ENERGYSPLITTING_ALG_C

#include "Algorithm/EnergySplittingAlg.h"

StatusCode EnergySplittingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  if(settings.map_floatPars.find("th_split")==settings.map_floatPars.end()) settings.map_floatPars["th_split"] = -1;
  if(settings.map_floatPars.find("Eth_SeedAbs")==settings.map_floatPars.end()) settings.map_floatPars["Eth_SeedAbs"] = 0.005;
  if(settings.map_stringPars.find("OutputLongiClusName")==settings.map_stringPars.end()) settings.map_stringPars["OutputLongiClusName"] = "ESLongiCluster";

  return StatusCode::SUCCESS;
};


StatusCode EnergySplittingAlg::Initialize(){

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

cout<<"EnergySplittingAlg: TowerCol size = "<<m_datacol.TowerCol.size()<<endl;

  std::vector<PandoraPlus::Calo3DCluster*>* p_3DClusters = &(m_datacol.Cluster3DCol);
  if( !p_3DClusters || p_3DClusters->size()==0){ std::cout<<"Warning: Empty or invalid input in EnergySplittingAlg! Please check previous algorithm!"<<std::endl; return StatusCode::SUCCESS; }
/*
for(int it=0; it<p_3DClusters->size(); it++){
cout<<"  Check tower #"<<it<<endl;
cout<<"  Block size: "<<p_3DClusters->at(it)->getBlocks().size()<<endl;
for(int ib=0; ib<p_3DClusters->at(it)->getBlocks().size(); ib++ ){

  const PandoraPlus::CaloBlock* p_block = p_3DClusters->at(it)->getBlocks()[ib];
  printf("    Block #%d: Layer = %d, bar size (%d, %d), shower size (%d, %d) \n", ib, p_block->getDlayer(),
              p_block->getBarUCol().size(), p_block->getBarVCol().size(), p_block->getShowerXCol().size(), p_block->getShowerYCol().size());

  std::vector<const CaloUnit *> barXCol = p_block->getBarUCol();
  std::vector<const CaloUnit *> barYCol = p_block->getBarVCol();
  std::vector<const CaloBarShower*> barShowerXCol = p_block->getShowerXCol();
  std::vector<const CaloBarShower*> barShowerYCol = p_block->getShowerYCol();

cout<<"Print Bar X:"<<endl;
for(int i=0; i<barXCol.size(); i++){
cout<<"  Address: "<<barXCol[i];
printf(", pos (%.2f, %.2f, %.2f, %.2f) \n", barXCol[i]->getPosition().x(), barXCol[i]->getPosition().y(), barXCol[i]->getPosition().z(), barXCol[i]->getEnergy() );
}
cout<<endl;
cout<<"Print Bar Y:"<<endl;
for(int i=0; i<barYCol.size(); i++){
cout<<"  Address: "<<barYCol[i];
printf(",  pos (%.2f, %.2f, %.2f, %.2f) \n", barYCol[i]->getPosition().x(), barYCol[i]->getPosition().y(), barYCol[i]->getPosition().z(), barYCol[i]->getEnergy() );
}
cout<<endl;
cout<<"Print shower X: "<<endl;
for(int i=0; i<barShowerXCol.size(); i++){
cout<<"  Address: "<<barShowerXCol[i];
printf(", Nbar = %d, pos (%.2f, %.2f, %.2f, %.2f) \n", barShowerXCol[i]->getBars().size(),
      barShowerXCol[i]->getPos().x(), barShowerXCol[i]->getPos().y(), barShowerXCol[i]->getPos().z(), barShowerXCol[i]->getEnergy());
}
cout<<endl;
cout<<"Print shower Y: "<<endl;
for(int i=0; i<barShowerYCol.size(); i++){
cout<<"  Address: "<<barShowerYCol[i];
printf(", Nbar = %d, pos (%.2f, %.2f, %.2f, %.2f) \n", barShowerYCol[i]->getBars().size(),
      barShowerYCol[i]->getPos().x(), barShowerYCol[i]->getPos().y(), barShowerYCol[i]->getPos().z(), barShowerYCol[i]->getEnergy());
}
cout<<endl;

//  for(int i=0; i<barXCol.size(); i++)
//    if(!barXCol[i]) { cout<<"WARNING: null ptr for barX "<<i<<endl; }
//  for(int i=0; i<barYCol.size(); i++)
//    if(!barYCol[i]) { cout<<"WARNING: null ptr for barY "<<i<<endl; }
//  for(int i=0; i<barShowerXCol.size(); i++)
//    if(!barShowerXCol[i]) { cout<<"WARNING: null ptr for showerX "<<i<<endl; }
//  for(int i=0; i<barShowerYCol.size(); i++)
//    if(!barShowerYCol[i]) { cout<<"WARNING: null ptr for showerY "<<i<<endl; }

}
}
*/

  for(int it=0; it<p_3DClusters->size(); it++){

    //Need to transfer const blocks to mutable. 
    std::vector<PandoraPlus::Calo2DCluster*> m_2Dclus; m_2Dclus.clear();
    for(int ib=0; ib<p_3DClusters->at(it)->getCluster().size(); ib++ )
      m_2Dclus.push_back( const_cast<PandoraPlus::Calo2DCluster *>(p_3DClusters->at(it)->getCluster()[ib]) );


    //for(int ib=0; ib<p_3DClusters->at(it)->getBlocks().size(); ib++){
    //  PandoraPlus::CaloBlock* m_block = new PandoraPlus::CaloBlock();  m_datacol.bk_BlockCol.push_back(m_block);
    //  *m_block = *(p_3DClusters->at(it)->getBlocks()[ib]);
    //  m_block->ClearShower();
    //  m_2Dclus.push_back(m_block);
    //}


    std::vector<const PandoraPlus::LongiCluster*> m_clusUCol; m_clusUCol.clear(); 
    std::vector<const PandoraPlus::LongiCluster*> m_clusVCol; m_clusVCol.clear(); 
    std::map<std::string, std::vector<const PandoraPlus::LongiCluster*> > tmp_mapU = p_3DClusters->at(it)->getLongiClusterUMap();
    std::map<std::string, std::vector<const PandoraPlus::LongiCluster*> > tmp_mapV = p_3DClusters->at(it)->getLongiClusterVMap();
    for(auto &iter : tmp_mapU){
      for(int icl=0; icl<iter.second.size(); icl++) m_clusUCol.push_back(iter.second[icl]);
    }
    for(auto &iter : tmp_mapV){
      for(int icl=0; icl<iter.second.size(); icl++) m_clusVCol.push_back(iter.second[icl]);
    }
    tmp_mapU.clear(); tmp_mapV.clear(); 
    if(m_clusUCol.size()==0 || m_clusVCol.size()==0) continue;


    int minLayerU = 99;
    int minLayerV = 99;
    int maxLayerU = -99;
    int maxLayerV = -99;
    for(int ic=0; ic<m_clusUCol.size(); ic++){
      if(minLayerU>m_clusUCol[ic]->getBeginningDlayer()) minLayerU = m_clusUCol[ic]->getBeginningDlayer();
      if(maxLayerU<m_clusUCol[ic]->getEndDlayer()) maxLayerU = m_clusUCol[ic]->getEndDlayer();
    }
    for(int ic=0; ic<m_clusVCol.size(); ic++){
      if(minLayerV>m_clusVCol[ic]->getBeginningDlayer()) minLayerV = m_clusVCol[ic]->getBeginningDlayer();
      if(maxLayerV<m_clusVCol[ic]->getEndDlayer()) maxLayerV = m_clusVCol[ic]->getEndDlayer();
    }

printf("  In Tower #%d: Block size = %d \n", it, m_2Dclus.size());
printf("  In Tower #%d: HoughCLX range: (%d, %d), HoughCLY range: (%d, %d) \n", it, minLayerU, maxLayerU, minLayerV, maxLayerV);

    //Make clusters and shwoers in blocks
    for(int ib=0; ib<m_2Dclus.size(); ib++ ){

      std::vector<PandoraPlus::Calo1DCluster*> m_barClusUCol; m_barClusUCol.clear();
      std::vector<PandoraPlus::Calo1DCluster*> m_barClusVCol; m_barClusVCol.clear();
      std::vector<const PandoraPlus::CaloUnit*> m_barColU = m_2Dclus[ib]->getBarUCol();
      std::vector<const PandoraPlus::CaloUnit*> m_barColV = m_2Dclus[ib]->getBarVCol();

      //Neighbor clustering
      int dlayer = m_2Dclus[ib]->getDlayer();


      if( dlayer>maxLayerU || dlayer<minLayerU ) Clustering( m_barColU, m_barClusUCol );
      else Clustering( m_barColU, m_barClusUCol, m_clusUCol );

      if( dlayer>maxLayerV || dlayer<minLayerV ) Clustering( m_barColV, m_barClusVCol );
      else Clustering( m_barColV, m_barClusVCol, m_clusVCol );

      m_datacol.bk_Cluster1DCol.insert( m_datacol.bk_Cluster1DCol.end(), m_barClusUCol.begin(), m_barClusUCol.end() );
      m_datacol.bk_Cluster1DCol.insert( m_datacol.bk_Cluster1DCol.end(), m_barClusVCol.begin(), m_barClusVCol.end() );

      //m_2Dclus[ib].setClusterXCol(m_barClusUCol);
      //m_2Dclus[ib].setClusterYCol(m_barClusVCol);

      //Split cluster to showers
      std::vector<const PandoraPlus::CaloBarShower*> m_barShowerUCol; m_barShowerUCol.clear();
      for(int ic=0; ic<m_barClusUCol.size(); ic++){
        std::vector<const PandoraPlus::CaloBarShower*> m_showers; m_showers.clear();
        ClusterSplitting( m_barClusUCol[ic], m_showers );
        if(m_showers.size()==0) continue;
        m_barShowerUCol.insert( m_barShowerUCol.end(), m_showers.begin(), m_showers.end() );
      }
      std::vector<const PandoraPlus::CaloBarShower*> m_barShowerVCol; m_barShowerVCol.clear();
      for(int ic=0; ic<m_barClusVCol.size(); ic++){
        std::vector<const PandoraPlus::CaloBarShower*> m_showers; m_showers.clear();
        ClusterSplitting( m_barClusVCol[ic], m_showers );
        if(m_showers.size()==0) continue;
        m_barShowerVCol.insert( m_barShowerVCol.end(), m_showers.begin(), m_showers.end() );
      }
     
      m_2Dclus[ib]->setShowerUCol( m_barShowerUCol );
      m_2Dclus[ib]->setShowerVCol( m_barShowerVCol );     

      for(int is=0; is<m_barShowerUCol.size(); is++) m_datacol.bk_BarShowerCol.push_back( const_cast<PandoraPlus::CaloBarShower*>(m_barShowerUCol[is]) );
      for(int is=0; is<m_barShowerVCol.size(); is++) m_datacol.bk_BarShowerCol.push_back( const_cast<PandoraPlus::CaloBarShower*>(m_barShowerVCol[is]) );

/*
printf("    Clustering in Block #%d: Layer = %d, bar size (%d, %d), cluster size (%d, %d), shower size (%d, %d) \n", ib, dlayer, m_barColU.size(), m_barColV.size(), m_barClusUCol.size(), m_barClusVCol.size(), m_barShowerUCol.size(), m_barShowerVCol.size());

cout<<"Print Bar X: "<<endl;
for(int i=0; i<m_barColU.size(); i++){
  cout<<"  Address: "<<m_barColU[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.3f) \n", m_barColU[i]->getPosition().x(), m_barColU[i]->getPosition().y(), m_barColU[i]->getPosition().z(), m_barColU[i]->getEnergy() );
}
cout<<"Print Bar X: "<<endl;
for(int i=0; i<m_barColV.size(); i++){
  cout<<"  Address: "<<m_barColV[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.3f) \n", m_barColV[i]->getPosition().x(), m_barColV[i]->getPosition().y(), m_barColV[i]->getPosition().z(), m_barColV[i]->getEnergy() );
}
cout<<"Print shower X: "<<endl;
for(int i=0; i<m_barShowerUCol.size(); i++){
  cout<<"  Address: "<<m_barShowerUCol[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.3f) \n", m_barShowerUCol[i]->getPos().x(), m_barShowerUCol[i]->getPos().y(), m_barShowerUCol[i]->getPos().z(), m_barShowerUCol[i]->getEnergy() );
}
cout<<"Print shower Y: "<<endl;
for(int i=0; i<m_barShowerVCol.size(); i++){
  cout<<"  Address: "<<m_barShowerVCol[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.3f) \n", m_barShowerVCol[i]->getPos().x(), m_barShowerVCol[i]->getPos().y(), m_barShowerVCol[i]->getPos().z(), m_barShowerVCol[i]->getEnergy() );
}
cout<<"Print cluster X: "<<endl;
for(int i=0; i<m_barClusUCol.size(); i++){
  cout<<"  Address: "<<m_barClusUCol[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.3f) \n", m_barClusUCol[i]->getPos().x(), m_barClusUCol[i]->getPos().y(), m_barClusUCol[i]->getPos().z(), m_barClusUCol[i]->getEnergy() );
}
cout<<"Print cluster Y: "<<endl;
for(int i=0; i<m_barClusVCol.size(); i++){
  cout<<"  Address: "<<m_barClusVCol[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.3f) \n", m_barClusVCol[i]->getPos().x(), m_barClusVCol[i]->getPos().y(), m_barClusVCol[i]->getPos().z(), m_barClusVCol[i]->getEnergy() );
}
*/
    }

    //Convert block to const
    //std::vector<const PandoraPlus::CaloBlock*> m_blocks; m_blocks.clear();
    //for(int ib=0; ib<m_2Dclus.size(); ib++)  m_blocks.push_back(m_2Dclus[ib]);

cout<<"  EnergySplittingAlg: Update clusters energy"<<endl;
    std::vector<const LongiCluster*> m_longiClusUCol; m_longiClusUCol.clear();
    std::vector<const LongiCluster*> m_longiClusVCol; m_longiClusVCol.clear();
    LongiClusterLinking(m_2Dclus, m_clusUCol, m_longiClusUCol);
    LongiClusterLinking(m_2Dclus, m_clusVCol, m_longiClusVCol);

/*
cout<<"  Loop check blocks AFTERLONGICLUSTERING:"<<endl;
for(int ib=0; ib<m_blocks.size(); ib++){
int dlayer = m_blocks[ib]->getDlayer();
std::vector<const CaloUnit *> m_barColU = m_blocks[ib]->getBarUCol();
std::vector<const CaloUnit *> m_barColV = m_blocks[ib]->getBarVCol();
std::vector<const CaloBarShower *> m_barShowerUCol = m_blocks[ib]->getShowerXCol();
std::vector<const CaloBarShower *> m_barShowerVCol = m_blocks[ib]->getShowerYCol();

printf("    Clustering in Block #%d: Layer = %d, bar size (%d, %d), shower size (%d, %d) \n", ib, dlayer, m_barColU.size(), m_barColV.size(), m_barShowerUCol.size(), m_barShowerVCol.size());

cout<<"Print Bar X: "<<endl;
for(int i=0; i<m_barColU.size(); i++){
  cout<<"  Address: "<<m_barColU[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.2f) \n", m_barColU[i]->getPosition().x(), m_barColU[i]->getPosition().y(), m_barColU[i]->getPosition().z(), m_barColU[i]->getEnergy() );
}
cout<<"Print Bar X: "<<endl;
for(int i=0; i<m_barColV.size(); i++){
  cout<<"  Address: "<<m_barColV[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.2f) \n", m_barColV[i]->getPosition().x(), m_barColV[i]->getPosition().y(), m_barColV[i]->getPosition().z(), m_barColV[i]->getEnergy() );
}
cout<<"Print shower X: "<<endl;
for(int i=0; i<m_barShowerUCol.size(); i++){
  cout<<"  Address: "<<m_barShowerUCol[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.2f) \n", m_barShowerUCol[i]->getPos().x(), m_barShowerUCol[i]->getPos().y(), m_barShowerUCol[i]->getPos().z(), m_barShowerUCol[i]->getEnergy() );
}
cout<<"Print shower Y: "<<endl;
for(int i=0; i<m_barShowerVCol.size(); i++){
  cout<<"  Address: "<<m_barShowerVCol[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.2f) \n", m_barShowerVCol[i]->getPos().x(), m_barShowerVCol[i]->getPos().y(), m_barShowerVCol[i]->getPos().z(), m_barShowerVCol[i]->getEnergy() );
}

}
*/

    //p_3DClusters->at(it)->CleanBlock();
    //p_3DClusters->at(it)->CleanLongiClusters();
    p_3DClusters->at(it)->setLongiClusters( settings.map_stringPars["OutputLongiClusName"], m_longiClusUCol, 
                                            settings.map_stringPars["OutputLongiClusName"], m_longiClusVCol );
    //p_3DClusters->at(it)->setBlocks( m_blocks );
  }


cout<<"End EnergySplittingAlg"<<endl;
  p_3DClusters = nullptr;
  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::LongiClusterLinking( std::vector<PandoraPlus::Calo2DCluster*>& m_blocks, 
                                                    std::vector<const PandoraPlus::LongiCluster*>& m_oldClusCol, 
                                                    std::vector<const PandoraPlus::LongiCluster*>& m_outClusCol )
{
  if(m_blocks.size()==0 || m_oldClusCol.size()==0) return StatusCode::SUCCESS;

  std::vector<PandoraPlus::LongiCluster*> m_cluscol; m_cluscol.clear();

  bool fl_isXclus = (m_oldClusCol[0]->getSlayer()==0);

  std::vector<const PandoraPlus::CaloBarShower*> m_showers; m_showers.clear(); 
  std::vector<const PandoraPlus::CaloBarShower*> m_showerswotrk; m_showerswotrk.clear(); 
//cout<<"  LongiClusterLinking: input block Dlayer: ";

  for(int ib=0; ib<m_blocks.size(); ib++){
//cout<<m_blocks[ib]->getDlayer()<<"(";

    std::vector<const PandoraPlus::CaloBarShower*> tmp_showersinblock; tmp_showersinblock.clear(); 
    if(fl_isXclus) tmp_showersinblock = m_blocks[ib]->getShowerUCol();
    else tmp_showersinblock = m_blocks[ib]->getShowerVCol();

    int dlayer = m_blocks[ib]->getDlayer();
    bool fl_hastrk = false; 
    for(int icl=0; icl<m_oldClusCol.size(); icl++) 
      if(dlayer<=m_oldClusCol[icl]->getEndDlayer() && dlayer>=m_oldClusCol[icl]->getBeginningDlayer()) { fl_hastrk=true; break; }
//cout<<fl_hastrk<<", shower size="<<tmp_showersinblock.size()<<")"<<'\t';

    if(fl_hastrk) m_showers.insert(m_showers.end(), tmp_showersinblock.begin(), tmp_showersinblock.end());
    else m_showerswotrk.insert(m_showerswotrk.end(), tmp_showersinblock.begin(), tmp_showersinblock.end());
  }

//cout<<endl;


//cout<<"  LongiClusterLinking: #normal shower = "<<m_showers.size()<<", #shower wo trk = "<<m_showerswotrk.size()<<endl;

//cout<<"  Normal shower: "<<endl;
//for(int is=0; is<m_showers.size(); is++){ 
//  cout<<"    Address: "<<m_showers[is];
//  printf(" pos (%.2f, %.2f, %.2f) \n", m_showers[is]->getPos().x(), m_showers[is]->getPos().y(), m_showers[is]->getPos().z());
//}

  //Update old Hough clusters
  for(int ic=0; ic<m_oldClusCol.size(); ic++){
    PandoraPlus::LongiCluster* m_newClus = new PandoraPlus::LongiCluster();
//cout<<"  Showers in cluster "<<ic<<endl;
    for(int is=0; is<m_oldClusCol[ic]->getBarShowers().size(); is++){
      const PandoraPlus::CaloBarShower* m_shower = m_oldClusCol[ic]->getBarShowers()[is];
//printf("    (%.2f, %.2f, %.2f) \n", m_shower->getPos().x(), m_shower->getPos().y(), m_shower->getPos().z());

      const PandoraPlus::CaloBarShower* m_selshower = NULL;
      bool fl_foundshower = false; 
      for(int js=0; js<m_showers.size(); js++){
        if( m_showers[js]->getModule() ==m_shower->getModule() && 
            m_showers[js]->getStave()  ==m_shower->getStave() && 
            m_showers[js]->getPart()   ==m_shower->getPart() && 
            m_showers[js]->getDlayer() ==m_shower->getDlayer() && 
            (m_showers[js]->getSeed()->getPosition()-m_shower->getPos()).Mag()<10 ) {m_selshower = m_showers[js]; fl_foundshower=true; break; }
      }
      if(fl_foundshower && m_selshower!=NULL) m_newClus->addBarShower( m_selshower, 1);
//else{ printf(" In Cluster #%d: shower #%d does not find the new shower. Layer = %d. \n ", ic, is, m_shower->getDlayer()); }
    }

    m_newClus->setHoughPars(m_oldClusCol[ic]->getHoughAlpha(), m_oldClusCol[ic]->getHoughRho() );
    m_newClus->setIntercept(m_oldClusCol[ic]->getHoughIntercept());
    //m_newClus->FitAxis();
    m_cluscol.push_back( m_newClus );
  }

/*
cout<<"    Before merging shower: printout LongiClusters "<<endl;
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

  //Merge showers into closest Hough cluster
  std::sort( m_showerswotrk.begin(), m_showerswotrk.end(), compLayer );
  for(int is=0; is<m_showerswotrk.size(); is++) 
    MergeToClosestCluster( m_showerswotrk[is], m_cluscol );

  //for(int ib=0; ib<m_blocks.size(); ib++) m_blocks[ib]->Check();

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


StatusCode EnergySplittingAlg::Clustering( std::vector<const PandoraPlus::CaloUnit*>& barCol, 
                                           std::vector<PandoraPlus::Calo1DCluster*>& outClus,  
                                           std::vector<const PandoraPlus::LongiCluster*>& m_longiClusCol )
{
  if(barCol.size()==0) return StatusCode::SUCCESS;

  int slayer = barCol[0]->getSlayer();
  //Neighbor clustering
  std::sort(barCol.begin(), barCol.end());
  std::vector<PandoraPlus::Calo1DCluster*> m_clusCol; m_clusCol.clear();
  for(int i=0;i<barCol.size();i++){
    if( (m_clusCol.size()!=0) && (m_clusCol[m_clusCol.size()-1]->isNeighbor(barCol[i])) &&
        ( (slayer==0 && m_clusCol[m_clusCol.size()-1]->getParts()[0]==barCol[i]->getPart()) ||
          (slayer==1 && m_clusCol[m_clusCol.size()-1]->getStaves()[0]==barCol[i]->getStave()) ) ){

      m_clusCol[m_clusCol.size()-1]->addCluster( barCol[i] );
      continue;
    }

    PandoraPlus::Calo1DCluster* clus = new PandoraPlus::Calo1DCluster();
    clus->addCluster(barCol[i]);
    m_clusCol.push_back(clus);
  }

//cout<<"Clustering: After Neighbor clustering. Cluster size: "<<m_clusCol.size()<<endl;
//cout<<"  Bar size in cluster: ";
//for(int j=0; j<m_clusCol.size(); j++) cout<<m_clusCol[j]->getBars().size()<<'\t';
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
        ( (slayer==0 && m_clusCol[m_clusCol.size()-1]->getParts()[0]==barCol[i]->getPart()) || 
          (slayer==1 && m_clusCol[m_clusCol.size()-1]->getStaves()[0]==barCol[i]->getStave()) ) ){

      m_clusCol[m_clusCol.size()-1]->addCluster( barCol[i] );
      continue;
    }

    PandoraPlus::Calo1DCluster* clus = new PandoraPlus::Calo1DCluster();
    clus->addCluster(barCol[i]);
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



StatusCode EnergySplittingAlg::ClusterSplitting( PandoraPlus::Calo1DCluster* m_cluster, std::vector<const PandoraPlus::CaloBarShower*>& outshCol ){
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
    PandoraPlus::CaloBarShower* shower = new PandoraPlus::CaloBarShower();
    shower->setBars(m_cluster->getBars());
    if(m_cluster->getNseeds()!=0) shower->setSeed(m_cluster->getSeeds()[0]);
    else shower->setSeed(nullptr);
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
    PandoraPlus::CaloBarShower* shower = new PandoraPlus::CaloBarShower();
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
    shower->setSeed(Bars[iseed]);
    shower->setIDInfo(); 

/*
cout<<"  Print shower #"<<is<<endl;
printf("    Shower Nbars = %d, pos/E (%.2f, %.2f, %.2f, %.3f) \n", shower->getBars().size(), 
                                                                   shower->getPos().x(), shower->getPos().y(), shower->getPos().z(), shower->getEnergy() );
printf("    Seed pos/E (%.2f, %.2f, %.2f, %.3f) \n", shower->getSeed()->getPosition().x(), shower->getSeed()->getPosition().y(), shower->getSeed()->getPosition().z(), shower->getSeed()->getEnergy() );
cout<<"    Bars: "<<endl;
for(int a=0; a<shower->getBars().size(); a++) 
  printf("      Bar pos/E (%.2f, %.2f, %.2f, %.3f) \n",  shower->getBars()[a]->getPosition().x(), 
                                                         shower->getBars()[a]->getPosition().y(),
                                                         shower->getBars()[a]->getPosition().z(),
                                                         shower->getBars()[a]->getEnergy() );
*/

    outshCol.push_back(shower);
  }
  
  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::MergeToClosestCluster( PandoraPlus::Calo1DCluster* iclus, std::vector<PandoraPlus::Calo1DCluster*>& clusvec ){

  int cLedge = iclus->getLeftEdge();
  int cRedge = iclus->getRightEdge();

  //Find the closest cluster with iclus. 
  int minD = 99;
  int index = -1;
  for(int icl=0; icl<clusvec.size(); icl++){
    if(clusvec[icl]->getNseeds()==0) continue;
    int iLedge = clusvec[icl]->getLeftEdge();
    int iRedge = clusvec[icl]->getRightEdge();

    int dis = (cLedge-iRedge>0 ? cLedge-iRedge : iLedge-cRedge );
    if(dis>5) continue; //Don't merge to a too far cluster. 
    if(dis<minD && clusvec[icl]->getEnergy()>2.*iclus->getEnergy() ){ minD = dis; index=icl; } //Don't merge to a too small cluster. 
  }
  if(index<0) return StatusCode::FAILURE;

  //Merge to the selected cluster
  for(int icl=0; icl<iclus->getBars().size(); icl++)  clusvec[index]->addCluster(iclus->getBars()[icl]);

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::MergeToClosestCluster( const PandoraPlus::CaloBarShower* m_shower, 
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
    if(m_neighbor.size()==0 && ibar->getEnergy()>settings.map_floatPars["Eth_SeedAbs"])
      seedCol.push_back(ibar); 

    else {
      bool isLocalMax=true; 
      for(int j=0;j<m_neighbor.size();j++){
        if(m_neighbor[j]->getEnergy()>ibar->getEnergy()) isLocalMax=false;
      }
      if(ibar->getEnergy()>settings.map_floatPars["Eth_SeedAbs"] && isLocalMax) seedCol.push_back(ibar);
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
