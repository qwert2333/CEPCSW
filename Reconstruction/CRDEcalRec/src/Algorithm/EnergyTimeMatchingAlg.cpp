#ifndef ETMATCHING_ALG_C
#define ETMATCHING_ALG_C

#include "Algorithm/EnergyTimeMatchingAlg.h"
StatusCode EnergyTimeMatchingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Set initial value
  if(settings.map_floatPars.find("chi2Wi_E")==settings.map_floatPars.end())         settings.map_floatPars["chi2Wi_E"] = 1.;
  if(settings.map_floatPars.find("chi2Wi_T")==settings.map_floatPars.end())         settings.map_floatPars["chi2Wi_T"] = 10.;
  if(settings.map_floatPars.find("sigmaE")==settings.map_floatPars.end())           settings.map_floatPars["sigmaE"] = 0.10;
  if(settings.map_floatPars.find("sigmaPos")==settings.map_floatPars.end())         settings.map_floatPars["sigmaPos"] = 34.89;
  if(settings.map_floatPars.find("nMat")==settings.map_floatPars.end())             settings.map_floatPars["nMat"] = 2.15;
  if(settings.map_floatPars.find("fl_UseChi2")==settings.map_floatPars.end())       settings.map_floatPars["fl_UseChi2"] = 1;
  if(settings.map_boolPars.find("fl_WriteCluster")==settings.map_boolPars.end())    settings.map_boolPars["fl_WriteCluster"] = 1;
  if(settings.map_floatPars.find("Debug")==settings.map_floatPars.end())            settings.map_floatPars["Debug"] = 0;

  return StatusCode::SUCCESS;
};

StatusCode EnergyTimeMatchingAlg::Initialize(){

	return StatusCode::SUCCESS;
};

StatusCode EnergyTimeMatchingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  std::vector<PandoraPlus::TransShower*> m_transhowerCol; m_transhowerCol.clear(); 
  std::vector<PandoraPlus::CaloCluster*> m_clusterCol; m_clusterCol.clear();

  //Basic unit: tower. 
  std::vector<PandoraPlus::CaloTower*>* p_towerCol = &(m_datacol.TowerCol);
  if(p_towerCol->size()==0){ std::cout<<"Warning: Empty input in EnergyTimeMatchingAlg! Please check previous algorithm!"<<std::endl; return StatusCode::SUCCESS; }

cout<<"Start EnergyTimeMatchingAlg. Tower size: "<<p_towerCol->size()<<endl;

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
printf(", pos (%.2f, %.2f, %.2f, %.2f), ID (%d, %d, %d) \n", barXCol[i]->getPosition().x(), barXCol[i]->getPosition().y(), barXCol[i]->getPosition().z(), barXCol[i]->getEnergy(), 
      barXCol[i]->getModule(), barXCol[i]->getStave(), barXCol[i]->getPart() );
}
cout<<endl;
cout<<"Print Bar Y:"<<endl;
for(int i=0; i<barYCol.size(); i++){
cout<<"  Address: "<<barYCol[i];
printf(",  pos (%.2f, %.2f, %.2f, %.2f), ID (%d, %d, %d) \n", barYCol[i]->getPosition().x(), barYCol[i]->getPosition().y(), barYCol[i]->getPosition().z(), barYCol[i]->getEnergy(),
      barYCol[i]->getModule(), barYCol[i]->getStave(), barYCol[i]->getPart() );
}
cout<<endl;

cout<<"Print shower X: "<<endl;
for(int i=0; i<barShowerXCol.size(); i++){
cout<<"  Address: "<<barShowerXCol[i];
printf(", Nbar = %d, pos (%.2f, %.2f, %.2f, %.2f), ID (%d, %d, %d) \n", barShowerXCol[i]->getBars().size(), 
      barShowerXCol[i]->getPos().x(), barShowerXCol[i]->getPos().y(), barShowerXCol[i]->getPos().z(), barShowerXCol[i]->getE(), 
      barShowerXCol[i]->getModule(), barShowerXCol[i]->getStave(), barShowerXCol[i]->getPart() );
}
cout<<endl;
cout<<"Print shower Y: "<<endl;
for(int i=0; i<barShowerYCol.size(); i++){
cout<<"  Address: "<<barShowerYCol[i];
printf(", Nbar = %d, pos (%.2f, %.2f, %.2f, %.2f), ID (%d, %d, %d) \n", barShowerYCol[i]->getBars().size(),
      barShowerYCol[i]->getPos().x(), barShowerYCol[i]->getPos().y(), barShowerYCol[i]->getPos().z(), barShowerYCol[i]->getE(), 
      barShowerYCol[i]->getModule(), barShowerYCol[i]->getStave(), barShowerYCol[i]->getPart() );
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

    //Case1: Don't use shadow cluster or chi2 matching, Save all combinations for 3D clustering.
    if(settings.map_floatPars["fl_UseChi2"]==0){
cout<<"  EnergyTimeMatchingAlg: Tower #"<<it<<", Not use chi2 matching"<<endl;
      std::vector<const PandoraPlus::CaloBlock*> m_blocks = p_towerCol->at(it)->getBlocks(); 
      if(m_blocks.size()==0) continue; 

      for(int ib=0;ib<m_blocks.size();ib++){
        std::vector<const PandoraPlus::CaloBarShower*> showerXCol = m_blocks[ib]->getShowerXCol();
        std::vector<const PandoraPlus::CaloBarShower*> showerYCol = m_blocks[ib]->getShowerYCol();
   
        std::vector<PandoraPlus::TransShower*> m_showerinlayer; m_showerinlayer.clear();

        GetFullMatchedShowers( showerXCol, showerYCol, m_showerinlayer );
        m_datacol.bk_TransShowerCol.insert( m_datacol.bk_TransShowerCol.end(), m_showerinlayer.begin(), m_showerinlayer.end() );
        
        m_transhowerCol.insert( m_transhowerCol.end(), m_showerinlayer.begin(), m_showerinlayer.end() );   
      }
    }
    

    //Case2: chi2 matching (combined chi2 in all layers)
    else{
cout<<"  EnergyTimeMatchingAlg: Tower #"<<it<<", Use chi2 matching"<<endl;
      std::vector<const PandoraPlus::LongiCluster*> m_longiClXCol = p_towerCol->at(it)->getLongiClusterXCol();
      std::vector<const PandoraPlus::LongiCluster*> m_longiClYCol = p_towerCol->at(it)->getLongiClusterYCol();

cout<<"  EnergyTimeMatchingAlg: Longitudinal cluster size: X "<<m_longiClXCol.size()<<", Y: "<<m_longiClYCol.size()<<endl;

      const int NclusX = m_longiClXCol.size(); 
      const int NclusY = m_longiClYCol.size(); 
      if(NclusX==0 || NclusY==0) continue; 

      std::vector<PandoraPlus::CaloCluster*> tmp_clusters; tmp_clusters.clear(); 
      
      //Case 2.1: 1*1
      if(NclusX==1 && NclusY==1){
        PandoraPlus::CaloCluster* tmp_clus = new PandoraPlus::CaloCluster();
        XYClusterMatchingL0(m_longiClXCol[0], m_longiClYCol[0], tmp_clus); 
        tmp_clusters.push_back(tmp_clus);
        //continue;
      }
      //Case 2.2: 1*N
      else if(NclusX==1){ XYClusterMatchingL1(m_longiClXCol[0], m_longiClYCol, tmp_clusters); }
      else if(NclusY==1){ XYClusterMatchingL1(m_longiClYCol[0], m_longiClXCol, tmp_clusters); }

      //Case 2.3: N*N
      else if( NclusX==NclusY ){ 
        XYClusterMatchingL2(m_longiClXCol, m_longiClYCol, tmp_clusters);
      }
      //Case 2.4
      else{ 
        XYClusterMatchingL3(m_longiClXCol, m_longiClYCol, tmp_clusters);
      }

cout<<"  EnergyTimeMatchingAlg: After cluster matching. 3D cluster size: "<<tmp_clusters.size()<<endl;

      for(int ic=0; ic<tmp_clusters.size(); ic++){ 
        m_clusterCol.push_back( tmp_clusters[ic] );
        std::vector<const PandoraPlus::TransShower*> m_showersinclus = tmp_clusters[ic]->getShowers(); 
cout<<"  EnergyTimeMatchingAlg: In Cluster #"<<ic<<": shower size = "<<m_showersinclus.size()<<endl;
        //m_transhowerCol.insert(m_transhowerCol.end(), m_showersinclus.begin(), m_showersinclus.end());
        for(int is=0; is<m_showersinclus.size(); is++){ 
          m_transhowerCol.push_back( const_cast<PandoraPlus::TransShower *>(m_showersinclus[is]) );
          m_datacol.bk_TransShowerCol.push_back( const_cast<PandoraPlus::TransShower *>(m_showersinclus[is]) );
        }

        m_datacol.bk_ClusterCol.push_back( tmp_clusters[ic] );
      }
    }
cout<<"  EnergyTimeMatchingAlg: present 2D shower size: "<<m_transhowerCol.size()<<", 3D cluster size: "<<m_clusterCol.size()<<endl;

  }

  p_towerCol = nullptr; 
  m_datacol.TransShowerCol = m_transhowerCol;
  if(settings.map_boolPars["fl_WriteCluster"]) m_datacol.map_CaloCluster["EcalCluster"] = m_clusterCol;

cout<<"End EnergyTimeMatchingAlg"<<endl;
  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}

//Longitudinal cluster: 1*1
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL0( const PandoraPlus::LongiCluster* m_longiClX, 
                                                       const PandoraPlus::LongiCluster* m_longiClY, 
                                                       PandoraPlus::CaloCluster* m_clus )
{
cout<<"  Cluster matching for case: 1 * 1"<<endl;

  std::vector<int> layerindex; layerindex.clear(); 
  std::map<int, std::vector<const PandoraPlus::CaloBarShower*> > map_showersXinlayer; map_showersXinlayer.clear();
  std::map<int, std::vector<const PandoraPlus::CaloBarShower*> > map_showersYinlayer; map_showersYinlayer.clear();

  for(int is=0; is<m_longiClX->getBarShowers().size(); is++){
    int m_layer = m_longiClX->getBarShowers()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);  
    map_showersXinlayer[m_layer].push_back(m_longiClX->getBarShowers()[is]);
  }
  for(int is=0; is<m_longiClY->getBarShowers().size(); is++){
    int m_layer = m_longiClY->getBarShowers()[is]->getDlayer(); 
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
    map_showersYinlayer[m_layer].push_back(m_longiClY->getBarShowers()[is]);
  }

  for(int il=0; il<layerindex.size(); il++){
    std::vector<const PandoraPlus::CaloBarShower*> m_showerXcol = map_showersXinlayer[layerindex[il]];
    std::vector<const PandoraPlus::CaloBarShower*> m_showerYcol = map_showersYinlayer[layerindex[il]];

printf("  In Layer %d: shower size (%d, %d) \n", layerindex[il], m_showerXcol.size(), m_showerYcol.size());

    std::vector<PandoraPlus::TransShower*> m_showerinlayer; m_showerinlayer.clear();

    if(m_showerXcol.size()==0 || m_showerYcol.size()==0) continue;
    else if(m_showerXcol.size()==1 && m_showerYcol.size()==1){ 
      PandoraPlus::TransShower* tmp_shower = new PandoraPlus::TransShower();
      GetMatchedShowersL0(m_showerXcol[0], m_showerYcol[0], tmp_shower); 
      m_showerinlayer.push_back(tmp_shower); 
    }
    else if(m_showerXcol.size()==1) GetMatchedShowersL1(m_showerXcol[0], m_showerYcol, m_showerinlayer );
    else if(m_showerYcol.size()==1) GetMatchedShowersL1(m_showerYcol[0], m_showerXcol, m_showerinlayer );
    else if(m_showerXcol.size()== m_showerYcol.size()) GetMatchedShowersL2(m_showerXcol, m_showerYcol, m_showerinlayer );
    else GetFullMatchedShowers(m_showerXcol, m_showerYcol, m_showerinlayer);
    //else GetMatchedShowersL3(m_showerXcol, m_showerYcol, m_showerinlayer);
cout<<"  After matching: shower size = "<<m_showerinlayer.size()<<endl;

    for(int is=0; is<m_showerinlayer.size(); is++) m_clus->addShower(m_showerinlayer[is]);
  }

  return StatusCode::SUCCESS;
}


//Longitudinal cluster: 1*N
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL1( const PandoraPlus::LongiCluster* m_longiCl1, 
                                                       std::vector<const PandoraPlus::LongiCluster*>& m_longiClN, 
                                                       std::vector<PandoraPlus::CaloCluster*>& m_clusters )
{

cout<<"  Input LongiCluster: "<<endl;
for(int ic=0; ic<m_longiClN.size(); ic++){
for(int is=0; is<m_longiClN[ic]->getBarShowers().size(); is++){
printf("    (%.3f, %.3f, %.3f, %.3f) \t", m_longiClN[ic]->getBarShowers()[is]->getPos().X(), 
                                          m_longiClN[ic]->getBarShowers()[is]->getPos().Y(),
                                          m_longiClN[ic]->getBarShowers()[is]->getPos().Z(), 
                                          m_longiClN[ic]->getBarShowers()[is]->getE() );
}
cout<<endl;
}

  m_clusters.clear(); 
  m_clusters.resize(m_longiClN.size());

  int slayer = m_longiCl1->getSlayer(); 
  std::vector<int> layerindex; layerindex.clear();
  std::map<int, std::vector<const PandoraPlus::CaloBarShower*> > map_showersXinlayer; map_showersXinlayer.clear();
  std::map<int, std::vector<const PandoraPlus::CaloBarShower*> > map_showersYinlayer; map_showersYinlayer.clear();

  for(int is=0; is<m_longiCl1->getBarShowers().size(); is++){
    int m_layer = m_longiCl1->getBarShowers()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);

    if(slayer==0) map_showersXinlayer[m_layer].push_back(m_longiCl1->getBarShowers()[is]);
    else map_showersYinlayer[m_layer].push_back(m_longiCl1->getBarShowers()[is]);
  }
  for(int ic=0; ic<m_longiClN.size(); ic++){
  for(int is=0; is<m_longiClN[ic]->getBarShowers().size(); is++){
    int m_layer = m_longiClN[ic]->getBarShowers()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);

    if(slayer==0) map_showersYinlayer[m_layer].push_back(m_longiClN[ic]->getBarShowers()[is]);
    else map_showersXinlayer[m_layer].push_back(m_longiClN[ic]->getBarShowers()[is]);
  }}
  
  std::map<int, std::vector<const PandoraPlus::TransShower*> > map_2Dshowersinlayer; map_2Dshowersinlayer.clear(); 

  sort(layerindex.begin(), layerindex.end());
  for(int il=0; il<layerindex.size(); il++){

    std::vector<const PandoraPlus::CaloBarShower*> m_showerXcol = map_showersXinlayer[layerindex[il]];
    std::vector<const PandoraPlus::CaloBarShower*> m_showerYcol = map_showersYinlayer[layerindex[il]];
cout<<"    XYClusterMatchingL1: In Layer #"<<layerindex[il]<<": (X, Y) = ("<<m_showerXcol.size()<<", "<<m_showerYcol.size()<<") "<<endl;

    std::vector<PandoraPlus::TransShower*> m_showerinlayer; m_showerinlayer.clear();

    if(m_showerXcol.size()==0 || m_showerYcol.size()==0) continue;
    else if(m_showerXcol.size()==1 && m_showerYcol.size()==1){ 
      PandoraPlus::TransShower* tmp_shower = new PandoraPlus::TransShower();
      GetMatchedShowersL0(m_showerXcol[0], m_showerYcol[0], tmp_shower); 
      m_showerinlayer.push_back(tmp_shower); 
    }
    else if(m_showerXcol.size()==1) GetMatchedShowersL1(m_showerXcol[0], m_showerYcol, m_showerinlayer );
    else if(m_showerYcol.size()==1) GetMatchedShowersL1(m_showerYcol[0], m_showerXcol, m_showerinlayer );
    else if(m_showerXcol.size() == m_showerYcol.size()) GetMatchedShowersL2(m_showerXcol, m_showerYcol, m_showerinlayer );
    //else{
    //  GetMatchedShowersL1( , , m_showerinlayer );
    //}
    else GetFullMatchedShowers(m_showerXcol, m_showerYcol, m_showerinlayer);
    //else GetMatchedShowersL3(m_showerXcol, m_showerYcol, m_showerinlayer);

    //map_2Dshowersinlayer[layerindex[il]] = m_showerinlayer; 
    for(int is=0; is<m_showerinlayer.size(); is++) map_2Dshowersinlayer[layerindex[il]].push_back( m_showerinlayer[is] );
  }


  for(int il=0; il<layerindex.size(); il++){
cout<<"    In Layer #"<<layerindex[il]<<":  TransShower size = "<<map_2Dshowersinlayer[layerindex[il]].size()<<endl;

  for(int is=0; is<map_2Dshowersinlayer[layerindex[il]].size(); is++){

printf("      TransShower #%d: pos+E=(%.3f, %.3f, %.3f) \n", is, map_2Dshowersinlayer[layerindex[il]][is]->getPos().z(), map_2Dshowersinlayer[layerindex[il]][is]->getPos().y(), map_2Dshowersinlayer[layerindex[il]][is]->getShowerE() );

    const PandoraPlus::CaloBarShower* m_barshower; 
    if(slayer==0) m_barshower=map_2Dshowersinlayer[layerindex[il]][is]->getShowerY();
    else          m_barshower=map_2Dshowersinlayer[layerindex[il]][is]->getShowerX();

//printf("      Selected barshower: DLayer=%d, sLayer=%d, pos+E=(%.3f, %.3f, %.3f) \n", m_barshower.getDlayer(), m_barshower.getSlayer(), m_barshower.getPos().Y(), m_barshower.getPos().Z(), m_barshower.getE() );

    int index_longiclus=-1; 
    bool fl_foundsh = false; 
    for(int ic=0; ic<m_longiClN.size() && !fl_foundsh; ic++){
    for(int jc=0; jc<m_longiClN[ic]->getBarShowers().size() && !fl_foundsh; jc++){
      if(m_barshower==m_longiClN[ic]->getBarShowers()[jc]) {index_longiclus=ic; fl_foundsh=true; break; }
    }}
//cout<<"      Found LongiShower: "<<index_longiclus<<endl;

    if(index_longiclus<0){ std::cout<<"WARNING: did not find properate longitudinal cluster! "<<std::endl; continue; }

    if( !m_clusters[index_longiclus] ){
      PandoraPlus::CaloCluster* p_newclus = new PandoraPlus::CaloCluster();
      p_newclus->addShower(map_2Dshowersinlayer[layerindex[il]][is]);
      m_clusters[index_longiclus] = p_newclus;
    }
    else m_clusters[index_longiclus]->addShower(map_2Dshowersinlayer[layerindex[il]][is]);
  }}

  return StatusCode::SUCCESS;
}


//Longitudinal cluster: N*N
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL2( std::vector<const PandoraPlus::LongiCluster*>& m_longiClXCol, 
                                                       std::vector<const PandoraPlus::LongiCluster*>& m_longiClYCol, 
                                                       std::vector<PandoraPlus::CaloCluster*>& m_clusters  )
{
  if(m_longiClXCol.size()==0 || m_longiClYCol.size()==0 || m_longiClXCol.size()!=m_longiClYCol.size()) return StatusCode::SUCCESS;

  std::vector<int> layerindex; layerindex.clear();
  std::map<int, std::vector<std::vector<const PandoraPlus::CaloBarShower*>> > map_showersXinlayer; map_showersXinlayer.clear();
  std::map<int, std::vector<std::vector<const PandoraPlus::CaloBarShower*>> > map_showersYinlayer; map_showersYinlayer.clear();

  //Find layers need to match.
//cout<<"  XYClusterMatchingL2: Find layers need to match "<<endl;
  for(int ic=0; ic<m_longiClXCol.size(); ic++){
  for(int is=0; is<m_longiClXCol[ic]->getBarShowers().size(); is++){
    int m_layer = m_longiClXCol[ic]->getBarShowers()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
  }}
  for(int ic=0; ic<m_longiClYCol.size(); ic++){
  for(int is=0; is<m_longiClYCol[ic]->getBarShowers().size(); is++){
    int m_layer = m_longiClYCol[ic]->getBarShowers()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
  }}
  sort(layerindex.begin(), layerindex.end());

  //Fill shower maps
//cout<<"  XYClusterMatchingL2: Fill shower maps"<<endl;  
  for(int il=0; il<layerindex.size(); il++){
    for(int ic=0; ic<m_longiClXCol.size(); ic++)
      map_showersXinlayer[layerindex[il]].push_back( m_longiClXCol[ic]->getBarShowersInLayer(layerindex[il]) );
    for(int ic=0; ic<m_longiClYCol.size(); ic++)
      map_showersYinlayer[layerindex[il]].push_back( m_longiClYCol[ic]->getBarShowersInLayer(layerindex[il]) );
  }
//cout<<"  XYClusterMatchingL2: map size X "<<map_showersXinlayer.size()<<", Y "<<map_showersYinlayer.size()<<endl;

  //Get the chi2 map for N*N
//cout<<"  XYClusterMatchingL2: Get chi2 maps"<<endl;
  const int Nclus = m_longiClXCol.size();
  double  **map_chi2[14] = {NULL};  //TODO: hardcoded layer number here! 
  double sumchi2[Nclus][Nclus] = {0};

  for(int il=0; il<layerindex.size(); il++){
//cout<<"  XYClusterMatchingL2: at layer #"<<layerindex[il]<<endl;
    std::vector<std::vector<const PandoraPlus::CaloBarShower*>> m_showerXcol = map_showersXinlayer[layerindex[il]];
    std::vector<std::vector<const PandoraPlus::CaloBarShower*>> m_showerYcol = map_showersYinlayer[layerindex[il]];

    //double **chi2inlayer = GetClusterChi2Map(m_showerXcol, m_showerYcol);
    //if(chi2inlayer!=nullptr) map_chi2[layerindex[il]-1] = chi2inlayer;
   map_chi2[layerindex[il]-1] = GetClusterChi2Map(m_showerXcol, m_showerYcol);
  }

//cout<<"  XYClusterMatchingL2: Print Sumchi2 matrix"<<endl;
  for(int ic=0; ic<Nclus; ic++){
  for(int jc=0; jc<Nclus; jc++){
    for(int il=0; il<14; il++){ 
      if(map_chi2[il]==nullptr) continue; 
      sumchi2[ic][jc] += map_chi2[il][ic][jc];
    }
//cout<<sumchi2[ic][jc]<<'\t';
  }
//cout<<endl;
  }

  //Get the chi2 of N! combinations
  int Ncomb=1;
  for(int i=Nclus; i>0; i--) Ncomb = Ncomb*i;

  map<double, vector<pair<int, int>> > matchingMap;
  int num[Nclus];
  int num_init[Nclus];
  for(int i=0;i<Nclus;i++){ num[i]=i; num_init[i]=i;}

  for(int icont=0;icont<Ncomb;icont++){
    vector<pair<int, int>> Index;
    for(int i=0;i<Nclus;i++){
       pair<int, int> p1(num_init[i], num[i]);
       Index.push_back(p1);
    }
    double chi2_tot=0;
    for(int i=0;i<Index.size();i++) chi2_tot += sumchi2[Index[i].first][Index[i].second];
    matchingMap[chi2_tot] = Index;

    Index.clear();
    if(!next_permutation(num, num+Nclus)) break;
  }
//cout<<"  XYClusterMatchingL2: Got chi2 of combination"<<endl;

  //map is ordered with [double] value, first element has the smallest chi2 value. 
  map<double, vector<pair<int, int>> >::iterator iter = matchingMap.begin();
  vector<pair<int, int>> Index = iter->second;

  for(int ii=0; ii<Index.size(); ii++){
//cout<<"  XYClusterMatchingL2: Selected combination: "<<Index[ii].first<<", "<<Index[ii].second<<endl;
    const PandoraPlus::LongiCluster* m_clusX = m_longiClXCol[Index[ii].first];
    const PandoraPlus::LongiCluster* m_clusY = m_longiClYCol[Index[ii].second];
    PandoraPlus::CaloCluster* tmp_clus = new PandoraPlus::CaloCluster();

    XYClusterMatchingL0(m_clusX, m_clusY, tmp_clus);
    m_clusters.push_back(tmp_clus);
  }

  return StatusCode::SUCCESS;
}


//Longitudinal cluster: M*N
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL3( std::vector<const PandoraPlus::LongiCluster*>& m_longiClXCol, 
                                                       std::vector<const PandoraPlus::LongiCluster*>& m_longiClYCol, 
                                                       std::vector<PandoraPlus::CaloCluster*>& m_clusters )
{



  return StatusCode::SUCCESS;
}



StatusCode EnergyTimeMatchingAlg::GetFullMatchedShowers(  std::vector<const PandoraPlus::CaloBarShower*>& barShowerXCol,
                                                          std::vector<const PandoraPlus::CaloBarShower*>& barShowerYCol,
                                                          std::vector<PandoraPlus::TransShower*>& outshCol )
{
  outshCol.clear();

  for(int is=0;is<barShowerXCol.size();is++){
  for(int js=0;js<barShowerYCol.size();js++){
    if(barShowerXCol[is]->getBars().size()==0 || barShowerYCol[js]->getBars().size()==0) continue;
    PandoraPlus::TransShower* tmp_shower = new PandoraPlus::TransShower();
    GetMatchedShowersL0(barShowerXCol[is], barShowerYCol[js], tmp_shower);
    outshCol.push_back(tmp_shower);
  }}

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL0( const PandoraPlus::CaloBarShower* barShowerX,   
                                                       const PandoraPlus::CaloBarShower* barShowerY,
                                                       PandoraPlus::TransShower* outsh )
{

  std::vector<const PandoraPlus::CaloHit*> m_digiCol; m_digiCol.clear();
  int NbarsX = barShowerX->getBars().size();
  int NbarsY = barShowerY->getBars().size();
  if(NbarsX==0 || NbarsY==0){ std::cout<<"WARNING: empty DigiHitsCol returned!"<<std::endl; return StatusCode::SUCCESS; }

  int _module = barShowerX->getBars()[0]->getModule();
  int _dlayer = barShowerX->getBars()[0]->getDlayer();
  int _part   = barShowerX->getBars()[0]->getPart();
  int _stave  = barShowerX->getBars()[0]->getStave();


  float rotAngle = -_module*PI/4.;

  for(int ibar=0;ibar<NbarsX;ibar++){
    double En_x = barShowerX->getBars()[ibar]->getEnergy();
    TVector3 m_vecx = barShowerX->getBars()[ibar]->getPosition();
    m_vecx.RotateZ(rotAngle);

    for(int jbar=0;jbar<NbarsY;jbar++){
      double En_y = barShowerY->getBars()[jbar]->getEnergy();
      TVector3 m_vecy = barShowerY->getBars()[jbar]->getPosition();
      m_vecy.RotateZ(rotAngle);

      TVector3 p_hit(m_vecy.x(), (m_vecx.y()+m_vecy.y())/2., m_vecx.z() );
      p_hit.RotateZ(-rotAngle);
      double m_Ehit = En_x*En_y/barShowerY->getE() + En_x*En_y/barShowerX->getE();
      //Create new CaloHit
      PandoraPlus::CaloHit* hit = new PandoraPlus::CaloHit();
      //hit.setCellID(0);
      hit->setPosition(p_hit);
      hit->setEnergy(m_Ehit);
      m_digiCol.push_back(hit);
    }
  }

  outsh->setBarShowers( barShowerX, barShowerY );
  outsh->setCaloHits( m_digiCol );
  outsh->setIDInfo( _module, _stave, _dlayer, _part);
cout<<"    End output shower"<<endl;

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL1( const PandoraPlus::CaloBarShower* shower1,
                                                       std::vector<const PandoraPlus::CaloBarShower*>& showerNCol, 
                                                       std::vector<PandoraPlus::TransShower*>& outshCol )
{

  outshCol.clear();

  int _slayer = shower1->getBars()[0]->getSlayer();

  const int NshY = showerNCol.size();
  double totE_shY = 0;
  double EshY[NshY] = {0};
  for(int is=0;is<NshY;is++){ EshY[is] = showerNCol[is]->getE(); totE_shY += EshY[is]; }
  for(int is=0;is<NshY;is++){
    double wi_E = EshY[is]/totE_shY;
    PandoraPlus::CaloBarShower* m_splitshower1 = new PandoraPlus::CaloBarShower();

    PandoraPlus::CaloUnit* m_wiseed = shower1->getSeed()->Clone();
    m_wiseed->setQ( wi_E*m_wiseed->getQ1(), wi_E*m_wiseed->getQ2() );

    std::vector<const PandoraPlus::CaloUnit*> m_wibars;
    for(int ib=0;ib<shower1->getBars().size();ib++){
      PandoraPlus::CaloUnit* m_wibar = shower1->getBars()[ib]->Clone();
      m_wibar->setQ(wi_E*m_wibar->getQ1(), wi_E*m_wibar->getQ2());
      m_wibars.push_back(m_wibar);
    }

    m_splitshower1->setBars( m_wibars );
    m_splitshower1->setSeed( m_wiseed );

    PandoraPlus::TransShower* m_shower = new PandoraPlus::TransShower();
    if(_slayer==0 ) GetMatchedShowersL0( m_splitshower1, showerNCol[is], m_shower);
    else            GetMatchedShowersL0( showerNCol[is], m_splitshower1, m_shower);
    outshCol.push_back( m_shower );
  }

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL2( std::vector<const PandoraPlus::CaloBarShower*>& barShowerXCol,
                                                       std::vector<const PandoraPlus::CaloBarShower*>& barShowerYCol,
                                                       std::vector<PandoraPlus::TransShower*>& outshCol )
{
  outshCol.clear();
  if(barShowerXCol.size() != barShowerYCol.size() ) return StatusCode::FAILURE;

  const int Nshower = barShowerXCol.size();
  double chi2[Nshower][Nshower];
  double chi2_E[Nshower][Nshower];
  double chi2_tx[Nshower][Nshower];
  double chi2_ty[Nshower][Nshower];

  double wi_E = settings.map_floatPars["chi2Wi_E"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);
  double wi_T = settings.map_floatPars["chi2Wi_T"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);

  TVector3 m_vec(0,0,0);
  double rotAngle = -(barShowerXCol[0]->getBars())[0]->getModule()*PI/4.;
  TVector3 Cblock((barShowerXCol[0]->getBars())[0]->getPosition().x(), (barShowerXCol[0]->getBars())[0]->getPosition().y(), (barShowerYCol[0]->getBars())[0]->getPosition().z());
  Cblock.RotateZ(rotAngle);

  for(int ix=0;ix<Nshower;ix++){
  for(int iy=0;iy<Nshower;iy++){
    const PandoraPlus::CaloBarShower* showerX = barShowerXCol[ix];
    const PandoraPlus::CaloBarShower* showerY = barShowerYCol[iy];

    double Ex = showerX->getE();
    double Ey = showerY->getE();
    chi2_E[ix][iy] = pow(fabs(Ex-Ey)/settings.map_floatPars["sigmaE"], 2);
    double PosTx = C*(showerY->getT1()-showerY->getT2())/(2*settings.map_floatPars["nMat"]) + showerY->getPos().z();
    chi2_tx[ix][iy] = pow( fabs(PosTx-showerX->getPos().z())/settings.map_floatPars["sigmaPos"], 2 );

    double PosTy = C*(showerX->getT1()-showerX->getT2())/(2*settings.map_floatPars["nMat"]);
    m_vec.SetXYZ(showerY->getPos().x(), showerY->getPos().y(), showerY->getPos().z());
    m_vec.RotateZ(rotAngle);
    chi2_ty[ix][iy] = pow( fabs(PosTy - (m_vec-Cblock).x() )/settings.map_floatPars["sigmaPos"], 2);

    chi2[ix][iy] = chi2_E[ix][iy]*wi_E + (chi2_tx[ix][iy]+chi2_ty[ix][iy])*wi_T ;

  }}

  int Ncomb=1;
  for(int i=Nshower; i>0; i--) Ncomb = Ncomb*i;
  
  map<double, vector<pair<int, int>> > matchingMap;
  int num[Nshower];
  int num_init[Nshower];
  for(int i=0;i<Nshower;i++){ num[i]=i; num_init[i]=i;}

  for(int icont=0;icont<Ncomb;icont++){
    vector<pair<int, int>> Index;
    for(int i=0;i<Nshower;i++){
       pair<int, int> p1(num_init[i], num[i]);
       Index.push_back(p1);
    }
    double chi2_tot=0;
    for(int i=0;i<Index.size();i++) chi2_tot += chi2[Index[i].first][Index[i].second];
    matchingMap[chi2_tot] = Index;

    Index.clear();
    if(!next_permutation(num, num+Nshower)) break;
  }

  map<double, vector<pair<int, int>> >::iterator iter = matchingMap.begin();
  vector<pair<int, int>> Index = iter->second;

  for(int i=0;i<Index.size();i++){
    const PandoraPlus::CaloBarShower* showerX = barShowerXCol[Index[i].first];
    const PandoraPlus::CaloBarShower* showerY = barShowerYCol[Index[i].second];

    PandoraPlus::TransShower* tmp_shower = new PandoraPlus::TransShower();
    GetMatchedShowersL0(showerX, showerY, tmp_shower);
    outshCol.push_back(tmp_shower);
  }

  return StatusCode::SUCCESS;
}

StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL3(  std::vector<const PandoraPlus::CaloBarShower*>& barShowerXCol, 
                                                        std::vector<const PandoraPlus::CaloBarShower*>& barShowerYCol, 
                                                        std::vector<PandoraPlus::TransShower*>& outshCol )
{


  return StatusCode::SUCCESS;
}


double** EnergyTimeMatchingAlg::GetClusterChi2Map( std::vector<std::vector<const PandoraPlus::CaloBarShower*>>& barShowerXCol, 
                                                   std::vector<std::vector<const PandoraPlus::CaloBarShower*>>& barShowerYCol )
{
  const int NclusX = barShowerXCol.size(); 
  const int NclusY = barShowerYCol.size(); 

  if(NclusX==0 || NclusY==0) return nullptr; 
  //If one longidutinal cluster is empty in this layer: skip this layer. 
  for(int icx=0; icx<NclusX; icx++) if(barShowerXCol[icx].size()==0) return nullptr;
  for(int icy=0; icy<NclusY; icy++) if(barShowerYCol[icy].size()==0) return nullptr;

  double  **chi2map = new double*[NclusX];

  double chi2map_E[NclusX][NclusY]={0};
  double chi2map_tx[NclusX][NclusY]={0};
  double chi2map_ty[NclusX][NclusY]={0};

  double wi_E = settings.map_floatPars["chi2Wi_E"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);
  double wi_T = settings.map_floatPars["chi2Wi_T"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);


  TVector3 m_vec(0,0,0);
  double rotAngle = -(barShowerXCol[0][0]->getBars())[0]->getModule()*PI/4.;
  TVector3 Ctower((barShowerXCol[0][0]->getBars())[0]->getPosition().x(), (barShowerXCol[0][0]->getBars())[0]->getPosition().y(), (barShowerYCol[0][0]->getBars())[0]->getPosition().z());
  Ctower.RotateZ(rotAngle);

  for(int ix=0;ix<NclusX;ix++){
  chi2map[ix] = new double[NclusY];
  for(int iy=0;iy<NclusY;iy++){
    std::vector<const PandoraPlus::CaloBarShower*> clusterX = barShowerXCol[ix];
    std::vector<const PandoraPlus::CaloBarShower*> clusterY = barShowerYCol[iy];

    double min_chi2E = 999;
    double min_chi2tx = 999;
    double min_chi2ty = 999;

    for(int icx=0; icx<clusterX.size(); icx++){
    for(int icy=0; icy<clusterY.size(); icy++){
      const PandoraPlus::CaloBarShower* showerX = clusterX[icx];      
      const PandoraPlus::CaloBarShower* showerY = clusterY[icy];      


      double Ex = showerX->getE();
      double Ey = showerY->getE();
      double chi2_E = pow(fabs(Ex-Ey)/settings.map_floatPars["sigmaE"], 2);
      double PosTx = C*(showerY->getT1()-showerY->getT2())/(2*settings.map_floatPars["nMat"]) + showerY->getPos().z();
      double chi2_tx = pow( fabs(PosTx-showerX->getPos().z())/settings.map_floatPars["sigmaPos"], 2 );
   
      double PosTy = C*(showerX->getT1()-showerX->getT2())/(2*settings.map_floatPars["nMat"]);
      m_vec = showerY->getPos();
      m_vec.RotateZ(rotAngle);
      double chi2_ty = pow( fabs(PosTy - (m_vec-Ctower).x() )/settings.map_floatPars["sigmaPos"], 2);

      if(chi2_E<min_chi2E) min_chi2E=chi2_E; 
      if(chi2_tx<min_chi2tx) min_chi2tx=chi2_tx; 
      if(chi2_ty<min_chi2ty) min_chi2ty=chi2_ty;
    }}

    chi2map_E[ix][iy] = min_chi2E; 
    chi2map_tx[ix][iy] = min_chi2tx;
    chi2map_ty[ix][iy] = min_chi2ty;
    chi2map[ix][iy] = chi2map_E[ix][iy]*wi_E + (chi2map_tx[ix][iy]+chi2map_ty[ix][iy])*wi_T ;

  }}

  return chi2map; 

};


#endif
