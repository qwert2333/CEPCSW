#ifndef _ENERGYSPLITTING_ALG_C
#define _ENERGYSPLITTING_ALG_C

#include "Algorithm/EnergySplittingAlg.h"

StatusCode EnergySplittingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  if(settings.map_floatPars.find("th_split")==settings.map_floatPars.end()) settings.map_floatPars["th_split"] = -1;
  if(settings.map_floatPars.find("Eth_SeedAbs")==settings.map_floatPars.end()) settings.map_floatPars["Eth_SeedAbs"] = 0.005;

  return StatusCode::SUCCESS;
};


StatusCode EnergySplittingAlg::Initialize(){

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

cout<<"EnergySplittingAlg: TowerCol size = "<<m_datacol.TowerCol.size()<<endl;

  std::vector<PandoraPlus::CaloTower*>* p_towerCol = &(m_datacol.TowerCol);
/*
for(int it=0; it<p_towerCol->size(); it++){
cout<<"  Check tower #"<<it<<endl;
cout<<"  Block size: "<<p_towerCol->at(it)->getBlocks().size()<<endl;
for(int ib=0; ib<p_towerCol->at(it)->getBlocks().size(); ib++ ){

  const PandoraPlus::CaloBlock* p_block = p_towerCol->at(it)->getBlocks()[ib];
  printf("    Block #%d: Layer = %d, bar size (%d, %d), shower size (%d, %d) \n", ib, p_block->getDlayer(),
              p_block->getBarXCol().size(), p_block->getBarYCol().size(), p_block->getShowerXCol().size(), p_block->getShowerYCol().size());

  std::vector<const CaloUnit *> barXCol = p_block->getBarXCol();
  std::vector<const CaloUnit *> barYCol = p_block->getBarYCol();
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
      barShowerXCol[i]->getPos().x(), barShowerXCol[i]->getPos().y(), barShowerXCol[i]->getPos().z(), barShowerXCol[i]->getE());
}
cout<<endl;
cout<<"Print shower Y: "<<endl;
for(int i=0; i<barShowerYCol.size(); i++){
cout<<"  Address: "<<barShowerYCol[i];
printf(", Nbar = %d, pos (%.2f, %.2f, %.2f, %.2f) \n", barShowerYCol[i]->getBars().size(),
      barShowerYCol[i]->getPos().x(), barShowerYCol[i]->getPos().y(), barShowerYCol[i]->getPos().z(), barShowerYCol[i]->getE());
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

  for(int it=0; it<p_towerCol->size(); it++){

    //Need to transfer const blocks to mutable. 
    std::vector<PandoraPlus::CaloBlock*> m_blocksInTower; m_blocksInTower.clear();
    for(int ib=0; ib<p_towerCol->at(it)->getBlocks().size(); ib++){
      PandoraPlus::CaloBlock* m_block = new PandoraPlus::CaloBlock();  m_datacol.bk_BlockCol.push_back(m_block);
      *m_block = *(p_towerCol->at(it)->getBlocks()[ib]);
      m_block->ClearShower();
      m_blocksInTower.push_back(m_block);
    }

    std::vector<const PandoraPlus::LongiCluster*> m_clusXCol = p_towerCol->at(it)->getLongiClusterXCol();
    std::vector<const PandoraPlus::LongiCluster*> m_clusYCol = p_towerCol->at(it)->getLongiClusterYCol();
    if(m_clusXCol.size()==0 || m_clusYCol.size()==0) continue;

    int minLayerX = 99;
    int minLayerY = 99;
    int maxLayerX = -99;
    int maxLayerY = -99;
    for(int ic=0; ic<m_clusXCol.size(); ic++){
      if(minLayerX>m_clusXCol[ic]->getBeginningDlayer()) minLayerX = m_clusXCol[ic]->getBeginningDlayer();
      if(maxLayerX<m_clusXCol[ic]->getEndDlayer()) maxLayerX = m_clusXCol[ic]->getEndDlayer();
    }
    for(int ic=0; ic<m_clusYCol.size(); ic++){
      if(minLayerY>m_clusYCol[ic]->getBeginningDlayer()) minLayerY = m_clusYCol[ic]->getBeginningDlayer();
      if(maxLayerY<m_clusYCol[ic]->getEndDlayer()) maxLayerY = m_clusYCol[ic]->getEndDlayer();
    }

printf("  In Tower #%d: Block size = %d \n", it, m_blocksInTower.size());
printf("  In Tower #%d: HoughCLX range: (%d, %d), HoughCLY range: (%d, %d) \n", it, minLayerX, maxLayerX, minLayerY, maxLayerY);

    //Make clusters and shwoers in blocks
    for(int ib=0; ib<m_blocksInTower.size(); ib++ ){

      std::vector<PandoraPlus::CaloBarCluster*> m_barClusXCol; m_barClusXCol.clear();
      std::vector<PandoraPlus::CaloBarCluster*> m_barClusYCol; m_barClusYCol.clear();
      std::vector<const PandoraPlus::CaloUnit*> m_barColX = m_blocksInTower[ib]->getBarXCol();
      std::vector<const PandoraPlus::CaloUnit*> m_barColY = m_blocksInTower[ib]->getBarYCol();

      //Neighbor clustering
      int dlayer = m_blocksInTower[ib]->getDlayer();


      if( dlayer>maxLayerX || dlayer<minLayerX ) Clustering( m_barColX, m_barClusXCol );
      else Clustering( m_barColX, m_barClusXCol, m_clusXCol );

      if( dlayer>maxLayerY || dlayer<minLayerY ) Clustering( m_barColY, m_barClusYCol );
      else Clustering( m_barColY, m_barClusYCol, m_clusYCol );

      m_datacol.bk_BarClusCol.insert( m_datacol.bk_BarClusCol.end(), m_barClusXCol.begin(), m_barClusXCol.end() );
      m_datacol.bk_BarClusCol.insert( m_datacol.bk_BarClusCol.end(), m_barClusYCol.begin(), m_barClusYCol.end() );

      //m_blocksInTower[ib].setClusterXCol(m_barClusXCol);
      //m_blocksInTower[ib].setClusterYCol(m_barClusYCol);
/*
if(dlayer==6){
printf("    Clustering in Block #%d: Layer = %d, bar size (%d, %d), cluster size (%d, %d) \n", ib, dlayer, m_barColX.size(), m_barColY.size(), m_barClusXCol.size(), m_barClusYCol.size() );
//cout<<"Print Bar X: "<<endl;
//for(int i=0; i<m_barColX.size(); i++){
//  cout<<"  Address: "<<m_barColX[i];
//  printf(",  pos (%.2f, %.2f, %.2f, %.2f) \n", m_barColX[i]->getPosition().x(), m_barColX[i]->getPosition().y(), m_barColX[i]->getPosition().z(), m_barColX[i]->getEnergy() );
//}
cout<<"Print Bar Y: "<<endl;
for(int i=0; i<m_barColY.size(); i++){
  cout<<"  Address: "<<m_barColY[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.2f) \n", m_barColY[i]->getPosition().x(), m_barColY[i]->getPosition().y(), m_barColY[i]->getPosition().z(), m_barColY[i]->getEnergy() );
}
cout<<"Print LongiClusterY: "<<endl;
for(int i=0; i<m_clusYCol.size(); i++){
  cout<<"  LongiClusterY #"<<i<<" covered layer: ";
  for(int al=0; al<m_clusYCol[i]->getBarShowers().size(); al++ ) cout<<m_clusYCol[i]->getBarShowers()[al]->getDlayer()<<'\t';
  cout<<endl;
}
//cout<<"Print cluster X: "<<endl;
//for(int i=0; i<m_barClusXCol.size(); i++){
//  cout<<"  Address: "<<m_barClusXCol[i];
//  printf(",  pos (%.2f, %.2f, %.2f, %.2f) \n", m_barClusXCol[i]->getPos().x(), m_barClusXCol[i]->getPos().y(), m_barClusXCol[i]->getPos().z(), m_barClusXCol[i]->getE() );
//}
cout<<"Print cluster Y: "<<endl;
for(int i=0; i<m_barClusYCol.size(); i++){
  cout<<"  Address: "<<m_barClusYCol[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.2f) \n", m_barClusYCol[i]->getPos().x(), m_barClusYCol[i]->getPos().y(), m_barClusYCol[i]->getPos().z(), m_barClusYCol[i]->getE() );
}

}
*/
      //Split cluster to showers
      std::vector<const PandoraPlus::CaloBarShower*> m_barShowerXCol; m_barShowerXCol.clear();
      for(int ic=0; ic<m_barClusXCol.size(); ic++){
        std::vector<const PandoraPlus::CaloBarShower*> m_showers; m_showers.clear();
        ClusterSplitting( m_barClusXCol[ic], m_showers );
        if(m_showers.size()==0) continue;
        m_barShowerXCol.insert( m_barShowerXCol.end(), m_showers.begin(), m_showers.end() );
      }
      std::vector<const PandoraPlus::CaloBarShower*> m_barShowerYCol; m_barShowerYCol.clear();
      for(int ic=0; ic<m_barClusYCol.size(); ic++){
        std::vector<const PandoraPlus::CaloBarShower*> m_showers; m_showers.clear();
        ClusterSplitting( m_barClusYCol[ic], m_showers );
        if(m_showers.size()==0) continue;
        m_barShowerYCol.insert( m_barShowerYCol.end(), m_showers.begin(), m_showers.end() );
      }
     
      m_blocksInTower[ib]->setShowerXCol( m_barShowerXCol );
      m_blocksInTower[ib]->setShowerYCol( m_barShowerYCol );     

      for(int is=0; is<m_barShowerXCol.size(); is++) m_datacol.bk_BarShowerCol.push_back( const_cast<PandoraPlus::CaloBarShower*>(m_barShowerXCol[is]) );
      for(int is=0; is<m_barShowerYCol.size(); is++) m_datacol.bk_BarShowerCol.push_back( const_cast<PandoraPlus::CaloBarShower*>(m_barShowerYCol[is]) );

/*
printf("    Clustering in Block #%d: Layer = %d, bar size (%d, %d), cluster size (%d, %d), shower size (%d, %d) \n", ib, dlayer, m_barColX.size(), m_barColY.size(), m_barClusXCol.size(), m_barClusYCol.size(), m_barShowerXCol.size(), m_barShowerYCol.size());

cout<<"Print Bar X: "<<endl;
for(int i=0; i<m_barColX.size(); i++){
  cout<<"  Address: "<<m_barColX[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.3f) \n", m_barColX[i]->getPosition().x(), m_barColX[i]->getPosition().y(), m_barColX[i]->getPosition().z(), m_barColX[i]->getEnergy() );
}
cout<<"Print Bar X: "<<endl;
for(int i=0; i<m_barColY.size(); i++){
  cout<<"  Address: "<<m_barColY[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.3f) \n", m_barColY[i]->getPosition().x(), m_barColY[i]->getPosition().y(), m_barColY[i]->getPosition().z(), m_barColY[i]->getEnergy() );
}
cout<<"Print shower X: "<<endl;
for(int i=0; i<m_barShowerXCol.size(); i++){
  cout<<"  Address: "<<m_barShowerXCol[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.3f) \n", m_barShowerXCol[i]->getPos().x(), m_barShowerXCol[i]->getPos().y(), m_barShowerXCol[i]->getPos().z(), m_barShowerXCol[i]->getE() );
}
cout<<"Print shower Y: "<<endl;
for(int i=0; i<m_barShowerYCol.size(); i++){
  cout<<"  Address: "<<m_barShowerYCol[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.3f) \n", m_barShowerYCol[i]->getPos().x(), m_barShowerYCol[i]->getPos().y(), m_barShowerYCol[i]->getPos().z(), m_barShowerYCol[i]->getE() );
}
cout<<"Print cluster X: "<<endl;
for(int i=0; i<m_barClusXCol.size(); i++){
  cout<<"  Address: "<<m_barClusXCol[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.3f) \n", m_barClusXCol[i]->getPos().x(), m_barClusXCol[i]->getPos().y(), m_barClusXCol[i]->getPos().z(), m_barClusXCol[i]->getE() );
}
cout<<"Print cluster Y: "<<endl;
for(int i=0; i<m_barClusYCol.size(); i++){
  cout<<"  Address: "<<m_barClusYCol[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.3f) \n", m_barClusYCol[i]->getPos().x(), m_barClusYCol[i]->getPos().y(), m_barClusYCol[i]->getPos().z(), m_barClusYCol[i]->getE() );
}
*/
    }

    //Convert block to const
    std::vector<const PandoraPlus::CaloBlock*> m_blocks; m_blocks.clear();
    for(int ib=0; ib<m_blocksInTower.size(); ib++)  m_blocks.push_back(m_blocksInTower[ib]);

cout<<"  EnergySplittingAlg: Update clusters energy"<<endl;
    std::vector<const LongiCluster*> m_longiClusXCol; m_longiClusXCol.clear();
    std::vector<const LongiCluster*> m_longiClusYCol; m_longiClusYCol.clear();
    LongiClusterLinking(m_blocks, m_clusXCol, m_longiClusXCol);
    LongiClusterLinking(m_blocks, m_clusYCol, m_longiClusYCol);

/*
cout<<"  Loop check blocks AFTERLONGICLUSTERING:"<<endl;
for(int ib=0; ib<m_blocks.size(); ib++){
int dlayer = m_blocks[ib]->getDlayer();
std::vector<const CaloUnit *> m_barColX = m_blocks[ib]->getBarXCol();
std::vector<const CaloUnit *> m_barColY = m_blocks[ib]->getBarYCol();
std::vector<const CaloBarShower *> m_barShowerXCol = m_blocks[ib]->getShowerXCol();
std::vector<const CaloBarShower *> m_barShowerYCol = m_blocks[ib]->getShowerYCol();

printf("    Clustering in Block #%d: Layer = %d, bar size (%d, %d), shower size (%d, %d) \n", ib, dlayer, m_barColX.size(), m_barColY.size(), m_barShowerXCol.size(), m_barShowerYCol.size());

cout<<"Print Bar X: "<<endl;
for(int i=0; i<m_barColX.size(); i++){
  cout<<"  Address: "<<m_barColX[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.2f) \n", m_barColX[i]->getPosition().x(), m_barColX[i]->getPosition().y(), m_barColX[i]->getPosition().z(), m_barColX[i]->getEnergy() );
}
cout<<"Print Bar X: "<<endl;
for(int i=0; i<m_barColY.size(); i++){
  cout<<"  Address: "<<m_barColY[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.2f) \n", m_barColY[i]->getPosition().x(), m_barColY[i]->getPosition().y(), m_barColY[i]->getPosition().z(), m_barColY[i]->getEnergy() );
}
cout<<"Print shower X: "<<endl;
for(int i=0; i<m_barShowerXCol.size(); i++){
  cout<<"  Address: "<<m_barShowerXCol[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.2f) \n", m_barShowerXCol[i]->getPos().x(), m_barShowerXCol[i]->getPos().y(), m_barShowerXCol[i]->getPos().z(), m_barShowerXCol[i]->getE() );
}
cout<<"Print shower Y: "<<endl;
for(int i=0; i<m_barShowerYCol.size(); i++){
  cout<<"  Address: "<<m_barShowerYCol[i];
  printf(",  pos (%.2f, %.2f, %.2f, %.2f) \n", m_barShowerYCol[i]->getPos().x(), m_barShowerYCol[i]->getPos().y(), m_barShowerYCol[i]->getPos().z(), m_barShowerYCol[i]->getE() );
}

}
*/

    p_towerCol->at(it)->CleanBlock();
    p_towerCol->at(it)->CleanLongiClusters();
    p_towerCol->at(it)->setLongiClusters( m_longiClusXCol, m_longiClusYCol );
    p_towerCol->at(it)->setBlocks( m_blocks );
  }


cout<<"End EnergySplittingAlg"<<endl;
  p_towerCol = nullptr;
  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::LongiClusterLinking( std::vector<const PandoraPlus::CaloBlock*>& m_blocks, 
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
    if(fl_isXclus) tmp_showersinblock = m_blocks[ib]->getShowerXCol();
    else tmp_showersinblock = m_blocks[ib]->getShowerYCol();

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
      if(fl_foundshower && m_selshower!=NULL) m_newClus->addBarShower( m_selshower );
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
    m_cluscol[aa]->getBarShowers()[ii]->getE(),
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
    m_cluscol[aa]->getBarShowers()[ii]->getE(),
    m_cluscol[aa]->getBarShowers()[ii]->getModule(), m_cluscol[aa]->getBarShowers()[ii]->getStave(), m_cluscol[aa]->getBarShowers()[ii]->getPart(),
    m_cluscol[aa]->getBarShowers()[ii]->getDlayer(), m_cluscol[aa]->getBarShowers()[ii]->getSlayer() );

}
*/

  //convert to const
  for(int icl=0; icl<m_cluscol.size(); icl++) m_outClusCol.push_back(m_cluscol[icl]);

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::Clustering( std::vector<const PandoraPlus::CaloUnit*>& barCol, 
                                           std::vector<PandoraPlus::CaloBarCluster*>& outClus,  
                                           std::vector<const PandoraPlus::LongiCluster*>& m_longiClusCol )
{
  if(barCol.size()==0) return StatusCode::SUCCESS;

  //Neighbor clustering
  std::sort(barCol.begin(), barCol.end(), compBar);
  std::vector<PandoraPlus::CaloBarCluster*> m_clusCol; m_clusCol.clear();
  for(int i=0;i<barCol.size();i++){
    if((m_clusCol.size()!=0) && (m_clusCol[m_clusCol.size()-1]->isNeighbor(barCol[i])) ){
      m_clusCol[m_clusCol.size()-1]->addBar( barCol[i] );
      continue;
    }

    PandoraPlus::CaloBarCluster* clus = new PandoraPlus::CaloBarCluster();
    clus->addBar(barCol[i]);
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
                                           std::vector<PandoraPlus::CaloBarCluster*>& outClus){
  if(barCol.size()==0) return StatusCode::SUCCESS;

  //Neighbor clustering
  std::sort(barCol.begin(), barCol.end(), compBar);
  std::vector<PandoraPlus::CaloBarCluster*> m_clusCol; m_clusCol.clear();
  for(int i=0;i<barCol.size();i++){
    if((m_clusCol.size()!=0) && (m_clusCol[m_clusCol.size()-1]->isNeighbor(barCol[i])) ){
      m_clusCol[m_clusCol.size()-1]->addBar( barCol[i] );
      continue;
    }

    PandoraPlus::CaloBarCluster* clus = new PandoraPlus::CaloBarCluster();
    clus->addBar(barCol[i]);
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



StatusCode EnergySplittingAlg::ClusterSplitting( PandoraPlus::CaloBarCluster* m_cluster, std::vector<const PandoraPlus::CaloBarShower*>& outshCol ){
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
                                                                   shower->getPos().x(), shower->getPos().y(), shower->getPos().z(), shower->getE() );
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


StatusCode EnergySplittingAlg::MergeToClosestCluster( PandoraPlus::CaloBarCluster* iclus, std::vector<PandoraPlus::CaloBarCluster*>& clusvec ){

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
    if(dis<minD && clusvec[icl]->getE()>2.*iclus->getE() ){ minD = dis; index=icl; } //Don't merge to a too small cluster. 
  }
  if(index<0) return StatusCode::FAILURE;

  //Merge to the selected cluster
  for(int icl=0; icl<iclus->getBars().size(); icl++)  clusvec[index]->addBar(iclus->getBars()[icl]);

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
//printf(" (%.3f, %.3f, %.3f, %.3f), layer #%d \n", sh_pos.X(), sh_pos.Y(), sh_pos.Z(), m_shower->getE(), dlayer);

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

  if(index_cluster>=0) m_clusters[index_cluster]->addBarShower( m_shower );
  //delete m_shower; m_shower = NULL;

  return StatusCode::SUCCESS;
}

StatusCode EnergySplittingAlg::findSeeds( PandoraPlus::CaloBarCluster* m_cluster, std::vector<const PandoraPlus::CaloUnit*>& seedCol){

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

std::vector<const PandoraPlus::CaloUnit*>  EnergySplittingAlg::getNeighbors(PandoraPlus::CaloBarCluster* m_cluster, const PandoraPlus::CaloUnit* seed){

  std::vector<const PandoraPlus::CaloUnit*> m_neighbor;
  for(int i=0;i<m_cluster->getBars().size();i++)
    if( seed->isNeighbor(m_cluster->getBars()[i]) ) m_neighbor.push_back(m_cluster->getBars()[i]);
  
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
