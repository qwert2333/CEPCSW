#ifndef ETMATCHING_ALG_C
#define ETMATCHING_ALG_C

#include "Algorithm/EnergyTimeMatchingAlg.h"
StatusCode EnergyTimeMatchingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Set initial value
  if(settings.map_floatPars.find("chi2Wi_E")==settings.map_floatPars.end())          settings.map_floatPars["chi2Wi_E"] = 1.;
  if(settings.map_floatPars.find("chi2Wi_T")==settings.map_floatPars.end())          settings.map_floatPars["chi2Wi_T"] = 10.;
  if(settings.map_floatPars.find("sigmaE")==settings.map_floatPars.end())            settings.map_floatPars["sigmaE"] = 0.10;
  if(settings.map_floatPars.find("sigmaPos")==settings.map_floatPars.end())          settings.map_floatPars["sigmaPos"] = 34.89;
  if(settings.map_floatPars.find("nMat")==settings.map_floatPars.end())              settings.map_floatPars["nMat"] = 2.15;
  if(settings.map_floatPars.find("fl_UseChi2")==settings.map_floatPars.end())        settings.map_floatPars["fl_UseChi2"] = 1;
  if(settings.map_boolPars.find("fl_WriteCluster")==settings.map_boolPars.end())     settings.map_boolPars["fl_WriteCluster"] = 1;
  if(settings.map_floatPars.find("Debug")==settings.map_floatPars.end())             settings.map_floatPars["Debug"] = 0;
  if(settings.map_stringPars.find("ReadinHFClusterName")==settings.map_stringPars.end()) settings.map_stringPars["ReadinHFClusterName"] = "ESHalfCluster";
  if(settings.map_stringPars.find("ReadinTowerName")==settings.map_stringPars.end()) settings.map_stringPars["ReadinTowerName"] = "ESTower";
  

  return StatusCode::SUCCESS;
};


StatusCode EnergyTimeMatchingAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  m_HFClusUCol.clear();
  m_HFClusVCol.clear();
  m_transhowerCol.clear();
  m_clusterCol.clear();
  m_towerCol.clear();

  m_towerCol = m_datacol.map_CaloCluster[settings.map_stringPars["ReadinTowerName"]];
cout<<"  EnergyTimeMatchingAlg: Readin tower size: "<<m_towerCol.size()<<endl;
	return StatusCode::SUCCESS;
};


StatusCode EnergyTimeMatchingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  m_clusterCol.clear();
  //Loop for towers:  
  for(int it=0; it<m_towerCol.size(); it++){
    m_HFClusUCol.clear(); m_HFClusVCol.clear();
 
    m_HFClusUCol = m_towerCol.at(it)->getHalfClusterUCol(settings.map_stringPars["ReadinHFClusterName"]+"U");
    m_HFClusVCol = m_towerCol.at(it)->getHalfClusterVCol(settings.map_stringPars["ReadinHFClusterName"]+"V");
//printf("  In Tower #%d: HalfCluster size (%d, %d), input total energy %.3f \n", it, m_HFClusUCol.size(), m_HFClusVCol.size(), m_towerCol[it]->getEnergy());
//cout<<"  Check track association: "<<endl;
//for(int icl=0; icl<m_HFClusUCol.size(); icl++) printf("    In HFClusU #%d: track size = %d \n", icl, m_HFClusUCol[icl]->getAssociatedTracks().size() );
//for(int icl=0; icl<m_HFClusVCol.size(); icl++) printf("    In HFClusV #%d: track size = %d \n", icl, m_HFClusVCol[icl]->getAssociatedTracks().size() );

    const int NclusX = m_HFClusUCol.size(); 
    const int NclusY = m_HFClusVCol.size(); 

    if(NclusX==0 || NclusY==0) continue;    
    std::vector<PandoraPlus::Calo3DCluster*> tmp_clusters; tmp_clusters.clear(); 
    
    //Case 2.1: 1*1
    if(NclusX==1 && NclusY==1){
      PandoraPlus::Calo3DCluster* tmp_clus = new PandoraPlus::Calo3DCluster();
      XYClusterMatchingL0(m_HFClusUCol[0], m_HFClusVCol[0], tmp_clus); 
      tmp_clusters.push_back(tmp_clus);
      //continue;
    }
    //Case 2.2: 1*N
    else if(NclusX==1){ XYClusterMatchingL1(m_HFClusUCol[0], m_HFClusVCol, tmp_clusters); }
    else if(NclusY==1){ XYClusterMatchingL1(m_HFClusVCol[0], m_HFClusUCol, tmp_clusters); }
    //Case 2.3: N*N
    else if( NclusX==NclusY ){ 
      XYClusterMatchingL2(m_HFClusUCol, m_HFClusVCol, tmp_clusters);
    }
    //Case 2.4 M*N
    else{ 
      XYClusterMatchingL3(m_HFClusUCol, m_HFClusVCol, tmp_clusters);
    }

    //Clean empty Calo3DClusters
    for(int ic=0; ic<tmp_clusters.size(); ic++){
      if(!tmp_clusters[ic]){
        delete tmp_clusters[ic]; tmp_clusters[ic]=NULL;
        tmp_clusters.erase(tmp_clusters.begin()+ic);
        ic--;
    }}

    //Save Calo3DClusters and Calo2DClusters into dataCol backupCol.
    for(int ic=0; ic<tmp_clusters.size(); ic++){
      m_clusterCol.push_back( tmp_clusters[ic] );
      std::vector<const PandoraPlus::Calo2DCluster*> m_showersinclus = tmp_clusters[ic]->getCluster();
      for(int is=0; is<m_showersinclus.size(); is++){
        m_transhowerCol.push_back( const_cast<PandoraPlus::Calo2DCluster *>(m_showersinclus[is]) );
        m_datacol.bk_Cluster2DCol.push_back( const_cast<PandoraPlus::Calo2DCluster *>(m_showersinclus[is]) );
      }
      m_datacol.bk_Cluster3DCol.push_back( tmp_clusters[ic] );
    }   

  }//End loop towers

cout<<"  Print reconstructed cluster energy and trk: "<<endl;
for(int i=0; i<m_clusterCol.size(); i++) cout<<"En = "<<m_clusterCol[i]->getEnergy()<<", trk size "<<m_clusterCol[i]->getAssociatedTracks().size()<<endl;


//cout<<"Cluster Merge Type1: from 1DShower aspect"<<endl;
  //Merge reconstructed cluster
  //  Type1: from 1DShower aspect
  for(int ic=0; ic<m_clusterCol.size() && m_clusterCol.size()>1; ic++){
    for(int jc=ic+1; jc<m_clusterCol.size(); jc++){
 
      int NlinkedU = 0; 
      int NlinkedV = 0;

      //Save out BarShowers(Calo1DCluster) in Cluster[jc]
      std::vector<const Calo1DCluster*> m_1DshowerUCol; m_1DshowerUCol.clear();
      std::vector<const Calo1DCluster*> m_1DshowerVCol; m_1DshowerVCol.clear();
      for(int is=0; is<m_clusterCol[jc]->getCluster().size(); is++){
        for(int js=0; js<m_clusterCol[jc]->getCluster()[is]->getShowerUCol().size(); js++)
          m_1DshowerUCol.push_back( m_clusterCol[jc]->getCluster()[is]->getShowerUCol()[js] );
        for(int js=0; js<m_clusterCol[jc]->getCluster()[is]->getShowerVCol().size(); js++)
          m_1DshowerVCol.push_back( m_clusterCol[jc]->getCluster()[is]->getShowerVCol()[js] );
      }

//printf("  1DShower size in Cl #%d: [%d, %d] \n", jc, m_1DshowerUCol.size(), m_1DshowerVCol.size());
//cout<<"  Check for cousins: "<<endl;
      //Loop in Cluster[ic], count the linked showers
      for(int is=0; is<m_clusterCol[ic]->getCluster().size(); is++){
        const Calo2DCluster* p_shower = m_clusterCol[ic]->getCluster()[is];

//printf("    In 3D Cluster %d 2DShower %d: CousinU size %d, CousinV size %d \n", ic, is, p_shower->getShowerUCol()[0]->getCousinClusters().size(), p_shower->getShowerVCol()[0]->getCousinClusters().size());

        for(int icc=0; icc<p_shower->getShowerUCol()[0]->getCousinClusters().size(); icc++)
          if( find(m_1DshowerUCol.begin(), m_1DshowerUCol.end(), p_shower->getShowerUCol()[0]->getCousinClusters()[icc])!=m_1DshowerUCol.end() ){ NlinkedU++; break; }

        for(int icc=0; icc<p_shower->getShowerVCol()[0]->getCousinClusters().size(); icc++)
          if( find(m_1DshowerVCol.begin(), m_1DshowerVCol.end(), p_shower->getShowerVCol()[0]->getCousinClusters()[icc])!=m_1DshowerVCol.end() ){ NlinkedV++; break; }

        p_shower = nullptr;
      }

//cout<<"  Found linked shower: "<<NlinkedU<<", "<<NlinkedV<<endl;
      if(NlinkedU/(float)m_clusterCol[ic]->getCluster().size()>0.7 || NlinkedV/(float)m_clusterCol[ic]->getCluster().size()>0.7){
        m_clusterCol[ic]->mergeCluster( m_clusterCol[jc] );
        m_clusterCol.erase(m_clusterCol.begin()+jc);
        jc--;
        if(jc<ic) jc=ic;
      }

  }}


//cout<<"Cluster Merge Type2: from HalfCluster aspect"<<endl;
  //  Type2: from longiCluster aspect
  for(int ic=0; ic<m_clusterCol.size() && m_clusterCol.size()>1; ic++){
    if( m_clusterCol[ic]->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster").size()==0 && 
        m_clusterCol[ic]->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster").size()==0 ) continue;
//cout<<"  Cluster ic has cousin"<<endl;
    for(int jc=ic+1; jc<m_clusterCol.size(); jc++){
      if( m_clusterCol[jc]->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster").size()==0 && 
          m_clusterCol[jc]->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster").size()==0 ) continue;
//cout<<"  Cluster jc also has cousin, start check"<<endl;

      std::vector<const PandoraPlus::CaloHalfCluster*> m_cousinU = m_clusterCol[ic]->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster"); 
      std::vector<const PandoraPlus::CaloHalfCluster*> m_cousinV = m_clusterCol[ic]->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster"); 
      if( find(m_cousinU.begin(), m_cousinU.end(), m_clusterCol[jc]->getHalfClusterUCol("LinkedLongiCluster")[0] )!=m_cousinU.end() ||
          find(m_cousinV.begin(), m_cousinV.end(), m_clusterCol[jc]->getHalfClusterVCol("LinkedLongiCluster")[0] )!=m_cousinV.end() ){
        m_clusterCol[ic]->mergeCluster( m_clusterCol[jc] );
        m_clusterCol.erase(m_clusterCol.begin()+jc);
        jc--;
        if(jc<ic) jc=ic;
      }        
  }}


  //  Type3: merge clusters linked to the same track. 
  for(int ic=0; ic<m_clusterCol.size() && m_clusterCol.size()>1; ic++){
    if(m_clusterCol[ic]->getAssociatedTracks().size()==0) continue;
    std::vector<const PandoraPlus::Track*> m_trkCol = m_clusterCol[ic]->getAssociatedTracks();

    for(int jc=ic+1; jc<m_clusterCol.size(); jc++){
      if(m_clusterCol[jc]->getAssociatedTracks().size()==0) continue;

      for(int itrk=0; itrk<m_clusterCol[jc]->getAssociatedTracks().size(); itrk++){
        if( find(m_trkCol.begin(), m_trkCol.end(), m_clusterCol[jc]->getAssociatedTracks()[itrk])!= m_trkCol.end() ){
          m_clusterCol[ic]->mergeCluster( m_clusterCol[jc] );
          m_clusterCol.erase(m_clusterCol.begin()+jc);
          jc--;
          if(jc<ic) jc=ic;
        }
      }
    }
  }

//cout<<"    After merge: cluster size "<<m_clusterCol.size()<<endl;



  m_datacol.map_2DCluster["EcalShowerInLayer"] = m_transhowerCol;
  if(settings.map_boolPars["fl_WriteCluster"]) m_datacol.map_CaloCluster["EcalCluster"] = m_clusterCol;

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}

//Longitudinal cluster: 1*1
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL0( const PandoraPlus::CaloHalfCluster* m_longiClU, 
                                                       const PandoraPlus::CaloHalfCluster* m_longiClV, 
                                                       PandoraPlus::Calo3DCluster* m_clus )
{
//cout<<"  Cluster matching for case: 1 * 1. Input HalfCluster En: "<<m_longiClU->getEnergy()<<", "<<m_longiClV->getEnergy()<<endl;
//cout<<"  Print 1DShower En in HalfClusterU: "<<endl;
//for(int i=0; i<m_longiClU->getCluster().size(); i++) 
//  cout<<m_longiClU->getCluster()[i]->getDlayer()<<'\t'<<m_longiClU->getCluster()[i]->getEnergy()<<endl;
//cout<<"  Print 1DShower En in HalfClusterV: "<<endl;
//for(int i=0; i<m_longiClV->getCluster().size(); i++) 
//  cout<<m_longiClV->getCluster()[i]->getDlayer()<<'\t'<<m_longiClV->getCluster()[i]->getEnergy()<<endl;


  std::vector<int> layerindex; layerindex.clear(); 
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > map_showersXinlayer; map_showersXinlayer.clear();
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > map_showersYinlayer; map_showersYinlayer.clear();

  for(int is=0; is<m_longiClU->getCluster().size(); is++){
    int m_layer = m_longiClU->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);  
    map_showersXinlayer[m_layer].push_back(m_longiClU->getCluster()[is]);
  }
  for(int is=0; is<m_longiClV->getCluster().size(); is++){
    int m_layer = m_longiClV->getCluster()[is]->getDlayer(); 
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
    map_showersYinlayer[m_layer].push_back(m_longiClV->getCluster()[is]);
  }

  for(int il=0; il<layerindex.size(); il++){
    std::vector<const PandoraPlus::Calo1DCluster*> m_showerXcol = map_showersXinlayer[layerindex[il]];
    std::vector<const PandoraPlus::Calo1DCluster*> m_showerYcol = map_showersYinlayer[layerindex[il]];

//printf("  In Layer %d: shower size (%d, %d) \n", layerindex[il], m_showerXcol.size(), m_showerYcol.size());

    std::vector<PandoraPlus::Calo2DCluster*> m_showerinlayer; m_showerinlayer.clear();

    if(m_showerXcol.size()==0 || m_showerYcol.size()==0) continue;
    else if(m_showerXcol.size()==1 && m_showerYcol.size()==1){ 
      PandoraPlus::Calo2DCluster* tmp_shower = new PandoraPlus::Calo2DCluster();
      GetMatchedShowersL0(m_showerXcol[0], m_showerYcol[0], tmp_shower); 
      m_showerinlayer.push_back(tmp_shower); 
    }
    else if(m_showerXcol.size()==1) GetMatchedShowersL1(m_showerXcol[0], m_showerYcol, m_showerinlayer );
    else if(m_showerYcol.size()==1) GetMatchedShowersL1(m_showerYcol[0], m_showerXcol, m_showerinlayer );
    else if(m_showerXcol.size()== m_showerYcol.size()) GetMatchedShowersL2(m_showerXcol, m_showerYcol, m_showerinlayer );
    //else GetFullMatchedShowers(m_showerXcol, m_showerYcol, m_showerinlayer);
    else GetMatchedShowersL3(m_showerXcol, m_showerYcol, m_showerinlayer);
//cout<<"    After matching: shower size = "<<m_showerinlayer.size()<<endl;
//for(int is=0; is<m_showerinlayer.size(); is++) cout<<m_showerinlayer[is]->getEnergy()<<'\t';
//cout<<endl;

    for(int is=0; is<m_showerinlayer.size(); is++) m_clus->addUnit(m_showerinlayer[is]);
  }
  m_clus->addHalfClusterU( "LinkedLongiCluster", m_longiClU );
  m_clus->addHalfClusterV( "LinkedLongiCluster", m_longiClV );
  for(auto itrk : m_longiClU->getAssociatedTracks()){
    if( find(m_clus->getAssociatedTracks().begin(), m_clus->getAssociatedTracks().end(), itrk)==m_clus->getAssociatedTracks().end() )
      m_clus->addAssociatedTrack(itrk);
  }
  for(auto itrk : m_longiClU->getAssociatedTracks()){
    if( find(m_clus->getAssociatedTracks().begin(), m_clus->getAssociatedTracks().end(), itrk)==m_clus->getAssociatedTracks().end() )
      m_clus->addAssociatedTrack(itrk);
  }

  return StatusCode::SUCCESS;
}


//Longitudinal cluster: 1*N
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL1( const PandoraPlus::CaloHalfCluster* m_longiCl1, 
                                                       std::vector<const PandoraPlus::CaloHalfCluster*>& m_longiClN, 
                                                       std::vector<PandoraPlus::Calo3DCluster*>& m_clusters )
{
  m_clusters.clear(); 
  m_clusters.resize(m_longiClN.size());

  int slayer = m_longiCl1->getSlayer(); 
  std::vector<int> layerindex; layerindex.clear();
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > map_showersXinlayer; map_showersXinlayer.clear();
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > map_showersYinlayer; map_showersYinlayer.clear();

  for(int is=0; is<m_longiCl1->getCluster().size(); is++){
    int m_layer = m_longiCl1->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);

    if(slayer==0) map_showersXinlayer[m_layer].push_back(m_longiCl1->getCluster()[is]);
    else map_showersYinlayer[m_layer].push_back(m_longiCl1->getCluster()[is]);
  }
  for(int ic=0; ic<m_longiClN.size(); ic++){
  for(int is=0; is<m_longiClN[ic]->getCluster().size(); is++){
    int m_layer = m_longiClN[ic]->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);

    if(slayer==0) map_showersYinlayer[m_layer].push_back(m_longiClN[ic]->getCluster()[is]);
    else map_showersXinlayer[m_layer].push_back(m_longiClN[ic]->getCluster()[is]);
  }}
  
  std::map<int, std::vector<const PandoraPlus::Calo2DCluster*> > map_2Dshowersinlayer; map_2Dshowersinlayer.clear(); 

  sort(layerindex.begin(), layerindex.end());
  for(int il=0; il<layerindex.size(); il++){

    std::vector<const PandoraPlus::Calo1DCluster*> m_showerXcol = map_showersXinlayer[layerindex[il]];
    std::vector<const PandoraPlus::Calo1DCluster*> m_showerYcol = map_showersYinlayer[layerindex[il]];

//cout<<"    XYClusterMatchingL1: In Layer #"<<layerindex[il]<<": (X, Y) = ("<<m_showerXcol.size()<<", "<<m_showerYcol.size()<<") "<<endl;

    std::vector<PandoraPlus::Calo2DCluster*> m_showerinlayer; m_showerinlayer.clear();

    if(m_showerXcol.size()==0 || m_showerYcol.size()==0) continue;
    else if(m_showerXcol.size()==1 && m_showerYcol.size()==1){ 
      PandoraPlus::Calo2DCluster* tmp_shower = new PandoraPlus::Calo2DCluster();
      GetMatchedShowersL0(m_showerXcol[0], m_showerYcol[0], tmp_shower); 
      m_showerinlayer.push_back(tmp_shower); 
    }
    else if(m_showerXcol.size()==1) GetMatchedShowersL1(m_showerXcol[0], m_showerYcol, m_showerinlayer );
    else if(m_showerYcol.size()==1) GetMatchedShowersL1(m_showerYcol[0], m_showerXcol, m_showerinlayer );
    else if(m_showerXcol.size() == m_showerYcol.size()) GetMatchedShowersL2(m_showerXcol, m_showerYcol, m_showerinlayer );
    //else{
    //  GetMatchedShowersL1( , , m_showerinlayer );
    //}
    //else GetFullMatchedShowers(m_showerXcol, m_showerYcol, m_showerinlayer);
    else GetMatchedShowersL3(m_showerXcol, m_showerYcol, m_showerinlayer);

    //map_2Dshowersinlayer[layerindex[il]] = m_showerinlayer; 
    for(int is=0; is<m_showerinlayer.size(); is++) map_2Dshowersinlayer[layerindex[il]].push_back( m_showerinlayer[is] );
  }

  //Longitudinal linking for 2D clusters
  for(int il=0; il<layerindex.size(); il++){
//cout<<"    In Layer #"<<layerindex[il]<<":  Calo2DCluster size = "<<map_2Dshowersinlayer[layerindex[il]].size()<<endl;

  for(int is=0; is<map_2Dshowersinlayer[layerindex[il]].size(); is++){

//printf("      Calo2DCluster #%d: pos+E=(%.3f, %.3f, %.3f) \n", is, map_2Dshowersinlayer[layerindex[il]][is]->getPos().z(), map_2Dshowersinlayer[layerindex[il]][is]->getPos().y(), map_2Dshowersinlayer[layerindex[il]][is]->getEnergy() );

    const PandoraPlus::Calo1DCluster* m_barshower; 
    if(slayer==1) m_barshower=map_2Dshowersinlayer[layerindex[il]][is]->getShowerUCol()[0];
    else          m_barshower=map_2Dshowersinlayer[layerindex[il]][is]->getShowerVCol()[0];

//printf("      Selected barshower: DLayer=%d, sLayer=%d, pos+E=(%.3f, %.3f, %.3f) \n", m_barshower->getDlayer(), m_barshower->getSlayer(), m_barshower->getPos().Y(), m_barshower->getPos().Z(), m_barshower->getEnergy() );
//printf("  Selected Shower (slayer = %d): Dlayer %d, En %.3f,  Address %p \n", slayer, m_barshower->getDlayer(), m_barshower->getEnergy(), m_barshower);

    int index_longiclus=-1; 
    bool fl_foundsh = false; 
    for(int ic=0; ic<m_longiClN.size() && !fl_foundsh; ic++){
    for(int jc=0; jc<m_longiClN[ic]->getCluster().size() && !fl_foundsh; jc++){
//const PandoraPlus::Calo1DCluster* tmp_sh = m_longiClN[ic]->getCluster()[jc];
//printf("    Check Shower in LongiCl: slayer %d, Dlater %d, En %.3f, Address %p \n", tmp_sh->getSlayer(), tmp_sh->getDlayer(), tmp_sh->getEnergy(), tmp_sh);
//tmp_sh = nullptr;
      if(m_barshower==m_longiClN[ic]->getCluster()[jc]) {index_longiclus=ic; fl_foundsh=true; break; }
    }}

    if(index_longiclus<0){ std::cout<<"WARNING: did not find properate longitudinal cluster! "<<std::endl; continue; }

    if( !m_clusters[index_longiclus] ){
      PandoraPlus::Calo3DCluster* p_newclus = new PandoraPlus::Calo3DCluster();
      p_newclus->addUnit(map_2Dshowersinlayer[layerindex[il]][is]);
      m_clusters[index_longiclus] = p_newclus;
    }
    else m_clusters[index_longiclus]->addUnit(map_2Dshowersinlayer[layerindex[il]][is]);
  }}

  for(int ic=0; ic<m_clusters.size(); ic++){
    if(!m_clusters[ic]) continue;
    if(slayer==0){ 
      m_clusters[ic]->addHalfClusterU( "LinkedLongiCluster", m_longiCl1 );
      m_clusters[ic]->addHalfClusterV( "LinkedLongiCluster", m_longiClN[ic] );

      for(auto itrk : m_longiCl1->getAssociatedTracks()){
        if( find(m_clusters[ic]->getAssociatedTracks().begin(), m_clusters[ic]->getAssociatedTracks().end(), itrk)==m_clusters[ic]->getAssociatedTracks().end() )
          m_clusters[ic]->addAssociatedTrack(itrk);
      }
      for(auto itrk : m_longiClN[ic]->getAssociatedTracks()){
        if( find(m_clusters[ic]->getAssociatedTracks().begin(), m_clusters[ic]->getAssociatedTracks().end(), itrk)==m_clusters[ic]->getAssociatedTracks().end() )
          m_clusters[ic]->addAssociatedTrack(itrk);
      }
    }
    else{          
      m_clusters[ic]->addHalfClusterU( "LinkedLongiCluster", m_longiClN[ic] );
      m_clusters[ic]->addHalfClusterV( "LinkedLongiCluster", m_longiCl1 );

      for(auto itrk : m_longiCl1->getAssociatedTracks()){
        if( find(m_clusters[ic]->getAssociatedTracks().begin(), m_clusters[ic]->getAssociatedTracks().end(), itrk)==m_clusters[ic]->getAssociatedTracks().end() )
          m_clusters[ic]->addAssociatedTrack(itrk);
      }
      for(auto itrk : m_longiClN[ic]->getAssociatedTracks()){
        if( find(m_clusters[ic]->getAssociatedTracks().begin(), m_clusters[ic]->getAssociatedTracks().end(), itrk)==m_clusters[ic]->getAssociatedTracks().end() )
          m_clusters[ic]->addAssociatedTrack(itrk);
      }
    }
  }

  return StatusCode::SUCCESS;
}


//Longitudinal cluster: N*N
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL2( std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClUCol, 
                                                       std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClVCol, 
                                                       std::vector<PandoraPlus::Calo3DCluster*>& m_clusters  )
{
  if(m_ClUCol.size()==0 || m_ClVCol.size()==0 || m_ClUCol.size()!=m_ClVCol.size()) return StatusCode::SUCCESS;

  std::vector<int> layerindex; layerindex.clear();
  std::map<int, std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> > map_showersXinlayer; map_showersXinlayer.clear();
  std::map<int, std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> > map_showersYinlayer; map_showersYinlayer.clear();

  //Find layers need to match.
//cout<<"  XYClusterMatchingL2: Find layers need to match "<<endl;
  for(int ic=0; ic<m_ClUCol.size(); ic++){
  for(int is=0; is<m_ClUCol[ic]->getCluster().size(); is++){
    int m_layer = m_ClUCol[ic]->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
  }}
  for(int ic=0; ic<m_ClVCol.size(); ic++){
  for(int is=0; is<m_ClVCol[ic]->getCluster().size(); is++){
    int m_layer = m_ClVCol[ic]->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
  }}
  sort(layerindex.begin(), layerindex.end());

  //Fill shower maps
//cout<<"  XYClusterMatchingL2: Fill shower maps"<<endl;  
  for(int il=0; il<layerindex.size(); il++){
    for(int ic=0; ic<m_ClUCol.size(); ic++)
      map_showersXinlayer[layerindex[il]].push_back( m_ClUCol[ic]->getClusterInLayer(layerindex[il]) );
    for(int ic=0; ic<m_ClVCol.size(); ic++)
      map_showersYinlayer[layerindex[il]].push_back( m_ClVCol[ic]->getClusterInLayer(layerindex[il]) );
  }
//cout<<"  XYClusterMatchingL2: map size X "<<map_showersXinlayer.size()<<", Y "<<map_showersYinlayer.size()<<endl;

  //Get the chi2 map for N*N
//cout<<"  XYClusterMatchingL2: Get chi2 maps"<<endl;
  const int Nclus = m_ClUCol.size();
  double  **map_chi2[PandoraPlus::CaloUnit::Nlayer] = {NULL};  
  double sumchi2[Nclus][Nclus] = {0};

  for(int il=0; il<layerindex.size(); il++){
//cout<<"  XYClusterMatchingL2: at layer #"<<layerindex[il]<<endl;
    std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> m_showerXcol = map_showersXinlayer[layerindex[il]];
    std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> m_showerYcol = map_showersYinlayer[layerindex[il]];

    //double **chi2inlayer = GetClusterChi2Map(m_showerXcol, m_showerYcol);
    //if(chi2inlayer!=nullptr) map_chi2[layerindex[il]-1] = chi2inlayer;
   map_chi2[layerindex[il]-1] = GetClusterChi2Map(m_showerXcol, m_showerYcol);
  }

//cout<<"  XYClusterMatchingL2: Print Sumchi2 matrix"<<endl;
  for(int ic=0; ic<Nclus; ic++){
  for(int jc=0; jc<Nclus; jc++){
    for(int il=0; il<14; il++){ 
      if(map_chi2[il]==nullptr) continue; 

      if(m_ClUCol[ic]->getHalfClusterCol("CousinCluster").size()!=0 && m_ClVCol[jc]->getHalfClusterCol("CousinCluster").size()!=0) 
        sumchi2[ic][jc] = 0.;
      if( m_ClUCol[ic]->getAssociatedTracks().size()==1 && m_ClVCol[jc]->getAssociatedTracks().size()==1 &&
          m_ClUCol[ic]->getAssociatedTracks()[0]==m_ClVCol[jc]->getAssociatedTracks()[0] )
        sumchi2[ic][jc] = 0.;
      else sumchi2[ic][jc] += map_chi2[il][ic][jc];
    }
  }
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
    const PandoraPlus::CaloHalfCluster* m_clusX = m_ClUCol[Index[ii].first];
    const PandoraPlus::CaloHalfCluster* m_clusY = m_ClVCol[Index[ii].second];
    PandoraPlus::Calo3DCluster* tmp_clus = new PandoraPlus::Calo3DCluster();

    XYClusterMatchingL0(m_clusX, m_clusY, tmp_clus);
    m_clusters.push_back(tmp_clus);
  }

  return StatusCode::SUCCESS;
}


//Longitudinal cluster: M*N
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL3( std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClUCol, 
                                                       std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClVCol, 
                                                       std::vector<PandoraPlus::Calo3DCluster*>& m_clusters )
{
  if( m_ClUCol.size()==0 || m_ClVCol.size()==0 ) return StatusCode::SUCCESS;
//cout<<"  XYClusterMatchingL3: HalfCluster size "<<m_ClUCol.size()<<", "<<m_ClVCol.size()<<endl;
/*
for(int icl=0; icl<m_ClUCol.size(); icl++){
  cout<<"      In HFClusU #"<<icl<<": shower size = "<<m_ClUCol[icl]->getCluster().size()<<endl;
  for(auto ish : m_ClUCol[icl]->getCluster()){
    printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d \n ", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(),  ish->getNseeds() );
    cout<<endl;
  }
}
cout<<endl;

for(int icl=0; icl<m_ClVCol.size(); icl++){
  cout<<"      In HFClusV #"<<icl<<": shower size = "<<m_ClVCol[icl]->getCluster().size()<<endl;
  for(auto ish : m_ClVCol[icl]->getCluster()){
    printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d \n ", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(),  ish->getNseeds() );
    cout<<endl;
  }
}
*/

  std::vector<int> layerindex; layerindex.clear();
  std::map<int, std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> > map_showersXinlayer; map_showersXinlayer.clear();
  std::map<int, std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> > map_showersYinlayer; map_showersYinlayer.clear();

  //Find layers need to match.
  for(int ic=0; ic<m_ClUCol.size(); ic++){
  for(int is=0; is<m_ClUCol[ic]->getCluster().size(); is++){
    int m_layer = m_ClUCol[ic]->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
  }}
  for(int ic=0; ic<m_ClVCol.size(); ic++){
  for(int is=0; is<m_ClVCol[ic]->getCluster().size(); is++){
    int m_layer = m_ClVCol[ic]->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
  }}
  sort(layerindex.begin(), layerindex.end());  

  //Fill shower maps
  for(int il=0; il<layerindex.size(); il++){
    for(int ic=0; ic<m_ClUCol.size(); ic++)
      map_showersXinlayer[layerindex[il]].push_back( m_ClUCol[ic]->getClusterInLayer(layerindex[il]) );
    for(int ic=0; ic<m_ClVCol.size(); ic++)
      map_showersYinlayer[layerindex[il]].push_back( m_ClVCol[ic]->getClusterInLayer(layerindex[il]) );
  }
//cout<<"  XYClusterMatchingL3: map size X "<<map_showersXinlayer.size()<<", Y "<<map_showersYinlayer.size()<<endl;

  //Get the chi2 map
  const int NclusU = m_ClUCol.size();
  const int NclusV = m_ClVCol.size();
  double  **map_chi2[PandoraPlus::CaloUnit::Nlayer] = {NULL};
  double sumchi2[NclusU][NclusV] = {0};

  for(int il=0; il<layerindex.size(); il++){
    std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> m_showerXcol = map_showersXinlayer[layerindex[il]];
    std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> m_showerYcol = map_showersYinlayer[layerindex[il]];
    map_chi2[layerindex[il]-1] = GetClusterChi2Map(m_showerXcol, m_showerYcol);
  }

//cout<<"  XYClusterMatchingL3: Print Sumchi2 matrix"<<endl;
  map<double, pair<int, int> > m_chi2Map; m_chi2Map.clear();
  for(int ic=0; ic<NclusU; ic++){
  for(int jc=0; jc<NclusV; jc++){
    for(int il=0; il<14; il++){
      if(map_chi2[il]==nullptr) continue;

      if(m_ClUCol[ic]->getHalfClusterCol("CousinCluster").size()!=0 && m_ClVCol[jc]->getHalfClusterCol("CousinCluster").size()!=0)
        sumchi2[ic][jc] = 0.;
      if( m_ClUCol[ic]->getAssociatedTracks().size()==1 && m_ClVCol[jc]->getAssociatedTracks().size()==1 &&
          m_ClUCol[ic]->getAssociatedTracks()[0]==m_ClVCol[jc]->getAssociatedTracks()[0] )
        sumchi2[ic][jc] = 0.;
      else sumchi2[ic][jc] += map_chi2[il][ic][jc];
    }
    pair<int, int> p1(ic, jc);
    m_chi2Map[sumchi2[ic][jc]] = p1;
//cout<<sumchi2[ic][jc]<<'\t';
  }
//cout<<endl;
  }

  //pickout pairs <indexVec> with smallest chi2:
  pair<int, int> lastpair;
  vector<pair<int, int>> indexVec; indexVec.clear();
  auto iter = m_chi2Map.begin();
  for(iter; iter!=m_chi2Map.end(); iter++){
    pair<int, int> indexpair = iter->second;
    bool inLine = false;
    bool isLast = false;
    for(int i=0; i<indexVec.size(); i++){
      if( indexVec.size() == min(NclusU, NclusV)-1 ) {lastpair = indexpair; isLast=true;  break; }
      if( indexpair.first == indexVec[i].first || indexpair.second == indexVec[i].second ) { inLine=true; break; }
    }
    if(isLast) break;
    if(inLine) continue;
    indexVec.push_back(indexpair);
  }
  if( indexVec.size()!=min(NclusU, NclusV)-1 )
    std::cout<<"ERROR in XYClusterMatchingL3: found pair size "<<indexVec.size()<<" does not equal to min shower size -1 "<<min(NclusU, NclusV)-1<<endl;

  //For the left: 1*N
  std::vector<const PandoraPlus::CaloHalfCluster*> leftHFClusCol; leftHFClusCol.clear();
  for(int i=0; i<max(NclusU, NclusV); i++){
    bool fl_exist = false;
    for(int j=0; j<indexVec.size(); j++){
      int m_index = NclusU>NclusV ? indexVec[j].first : indexVec[j].second;
      if(i==m_index){ fl_exist = true; break; }
    }
    if(!fl_exist){
       const PandoraPlus::CaloHalfCluster* m_shower = NclusU>NclusV ? m_ClUCol[i] : m_ClVCol[i] ;
       leftHFClusCol.push_back(m_shower);
    }
  }
  if(leftHFClusCol.size() != fabs( NclusU-NclusV )+1 )
    std::cout<<"ERROR in XYClusterMatchingL3: Last pair number "<<leftHFClusCol.size()<<" does not equal to shower difference "<<fabs( NclusU-NclusV )+1<<std::endl;  

  //Match the pairs in the indexVec: 
//cout<<"  XYClusterMatchingL3: Match the pairs in the indexVec"<<endl;
  for(int ip=0; ip<indexVec.size(); ip++){
    PandoraPlus::Calo3DCluster* tmp_clus = new PandoraPlus::Calo3DCluster();
    XYClusterMatchingL0(m_ClUCol[indexVec[ip].first], m_ClVCol[indexVec[ip].second], tmp_clus);
    m_clusters.push_back(tmp_clus);
  }

//cout<<"  Last pair: "<<lastpair.first<<" "<<lastpair.second<<endl;
  //Match the left 1*N: 
//cout<<"  XYClusterMatchingL3: Match the left 1*N"<<endl;
  std::vector<PandoraPlus::Calo3DCluster*> tmp_3dclus; tmp_3dclus.clear();
  int ilast = NclusU<NclusV ? lastpair.first : lastpair.second;
  const PandoraPlus::CaloHalfCluster* m_clus = NclusU<NclusV ? m_ClUCol[ilast] : m_ClVCol[ilast];
  XYClusterMatchingL1(m_clus, leftHFClusCol, tmp_3dclus );

  m_clusters.insert(m_clusters.end(), tmp_3dclus.begin(), tmp_3dclus.end());
  m_clus = nullptr; 

  return StatusCode::SUCCESS;
}



StatusCode EnergyTimeMatchingAlg::GetFullMatchedShowers(  std::vector<const PandoraPlus::Calo1DCluster*>& barShowerUCol,
                                                          std::vector<const PandoraPlus::Calo1DCluster*>& barShowerVCol,
                                                          std::vector<PandoraPlus::Calo2DCluster*>& outshCol )
{
  outshCol.clear();

  for(int is=0;is<barShowerUCol.size();is++){
  for(int js=0;js<barShowerVCol.size();js++){
    if(barShowerUCol[is]->getBars().size()==0 || barShowerVCol[js]->getBars().size()==0) continue;
    PandoraPlus::Calo2DCluster* tmp_shower = new PandoraPlus::Calo2DCluster();
    GetMatchedShowersL0(barShowerUCol[is], barShowerVCol[js], tmp_shower);
    outshCol.push_back(tmp_shower);
  }}

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL0( const PandoraPlus::Calo1DCluster* barShowerU,   
                                                       const PandoraPlus::Calo1DCluster* barShowerV,
                                                       PandoraPlus::Calo2DCluster* outsh )
{

  std::vector<const PandoraPlus::CaloHit*> m_digiCol; m_digiCol.clear();
  int NbarsX = barShowerU->getBars().size();
  int NbarsY = barShowerV->getBars().size();
  if(NbarsX==0 || NbarsY==0){ std::cout<<"WARNING: empty DigiHitsCol returned!"<<std::endl; return StatusCode::SUCCESS; }
  if(barShowerU->getTowerID().size()==0) { std::cout<<"WARNING:GetMatchedShowersL0  No TowerID in 1DCluster!"<<std::endl; return StatusCode::SUCCESS; }
  //if(barShowerU->getTowerID().size()==0) { barShowerU->setIDInfo(); }

  int _module = barShowerU->getTowerID()[0][0];
  int _part   = barShowerU->getTowerID()[0][1];
  int _stave  = barShowerU->getTowerID()[0][2];
  int _dlayer = barShowerU->getBars()[0]->getDlayer();


  float rotAngle = -_module*PI/4.;

  for(int ibar=0;ibar<NbarsX;ibar++){
    double En_x = barShowerU->getBars()[ibar]->getEnergy();
    TVector3 m_vecx = barShowerU->getBars()[ibar]->getPosition();
    m_vecx.RotateZ(rotAngle);

    for(int jbar=0;jbar<NbarsY;jbar++){
      double En_y = barShowerV->getBars()[jbar]->getEnergy();
      TVector3 m_vecy = barShowerV->getBars()[jbar]->getPosition();
      m_vecy.RotateZ(rotAngle);

      TVector3 p_hit(m_vecy.x(), (m_vecx.y()+m_vecy.y())/2., m_vecx.z() );
      p_hit.RotateZ(-rotAngle);
      double m_Ehit = En_x*En_y/barShowerV->getEnergy() + En_x*En_y/barShowerU->getEnergy();
      //Create new CaloHit
      PandoraPlus::CaloHit* hit = new PandoraPlus::CaloHit();
      //hit.setCellID(0);
      hit->setPosition(p_hit);
      hit->setEnergy(m_Ehit);
      m_digiCol.push_back(hit);
    }
  }

  outsh->addUnit( barShowerU );
  outsh->addUnit( barShowerV );
  outsh->setCaloHits( m_digiCol );
//cout<<"    End output shower"<<endl;

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL1( const PandoraPlus::Calo1DCluster* shower1,
                                                       std::vector<const PandoraPlus::Calo1DCluster*>& showerNCol, 
                                                       std::vector<PandoraPlus::Calo2DCluster*>& outshCol )
{
//cout<<"  GetMatchedShowersL1: input shower size: 1 * "<<showerNCol.size()<<endl;
  outshCol.clear();

  int _slayer = shower1->getBars()[0]->getSlayer();

  const int NshY = showerNCol.size();
  double totE_shY = 0;
  double EshY[NshY] = {0};
  for(int is=0;is<NshY;is++){ EshY[is] = showerNCol[is]->getEnergy(); totE_shY += EshY[is]; }
  for(int is=0;is<NshY;is++){
    double wi_E = EshY[is]/totE_shY;
    PandoraPlus::Calo1DCluster* m_splitshower1 = new PandoraPlus::Calo1DCluster();

    PandoraPlus::CaloUnit* m_wiseed = nullptr; 
    if(shower1->getSeeds().size()>0) m_wiseed = shower1->getSeeds()[0]->Clone();
    else{ cout<<"ERROR: Input shower has no seed! Check! Use the most energitic bar as seed. bar size: "<<shower1->getBars().size()<<endl; 
      double m_maxE = -99;
      int index = -1; 
      for(int ib=0; ib<shower1->getBars().size(); ib++){
        if(shower1->getBars()[ib]->getEnergy()>m_maxE) { m_maxE=shower1->getBars()[ib]->getEnergy(); index=ib; }
      }
      if(index>=0) m_wiseed = shower1->getBars()[index]->Clone(); 
    }
    m_wiseed->setQ( wi_E*m_wiseed->getQ1(), wi_E*m_wiseed->getQ2() );

    std::vector<const PandoraPlus::CaloUnit*> m_wibars; m_wibars.clear(); 
    for(int ib=0;ib<shower1->getBars().size();ib++){
      PandoraPlus::CaloUnit* m_wibar = shower1->getBars()[ib]->Clone();
      m_wibar->setQ(wi_E*m_wibar->getQ1(), wi_E*m_wibar->getQ2());
      m_wibars.push_back(m_wibar);
    }
    m_splitshower1->setBars( m_wibars );
    m_splitshower1->addSeed( m_wiseed );
    m_splitshower1->setIDInfo();
    PandoraPlus::Calo2DCluster* m_shower = new PandoraPlus::Calo2DCluster();
    if(_slayer==0 ) GetMatchedShowersL0( m_splitshower1, showerNCol[is], m_shower);
    else            GetMatchedShowersL0( showerNCol[is], m_splitshower1, m_shower);
    outshCol.push_back( m_shower );
  }
  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL2( std::vector<const PandoraPlus::Calo1DCluster*>& barShowerUCol,
                                                       std::vector<const PandoraPlus::Calo1DCluster*>& barShowerVCol,
                                                       std::vector<PandoraPlus::Calo2DCluster*>& outshCol )
{
  outshCol.clear();
  if(barShowerUCol.size() != barShowerVCol.size() ) return StatusCode::FAILURE;

  const int Nshower = barShowerUCol.size();
  double chi2[Nshower][Nshower];
  double chi2_E[Nshower][Nshower];
  double chi2_tx[Nshower][Nshower];
  double chi2_ty[Nshower][Nshower];

  double wi_E = settings.map_floatPars["chi2Wi_E"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);
  double wi_T = settings.map_floatPars["chi2Wi_T"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);

  TVector3 m_vec(0,0,0);
  double rotAngle = -(barShowerUCol[0]->getBars())[0]->getModule()*PI/4.;
  TVector3 Cblock( (barShowerUCol[0]->getBars())[0]->getPosition().x(), 
                   (barShowerUCol[0]->getBars())[0]->getPosition().y(), 
                   (barShowerVCol[0]->getBars())[0]->getPosition().z() );
  Cblock.RotateZ(rotAngle);

  for(int ix=0;ix<Nshower;ix++){
  for(int iy=0;iy<Nshower;iy++){
    const PandoraPlus::Calo1DCluster* showerX = barShowerUCol[ix];
    const PandoraPlus::Calo1DCluster* showerY = barShowerVCol[iy];

    double Ex = showerX->getEnergy();
    double Ey = showerY->getEnergy();
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
    const PandoraPlus::Calo1DCluster* showerX = barShowerUCol[Index[i].first];
    const PandoraPlus::Calo1DCluster* showerY = barShowerVCol[Index[i].second];

    PandoraPlus::Calo2DCluster* tmp_shower = new PandoraPlus::Calo2DCluster();
    GetMatchedShowersL0(showerX, showerY, tmp_shower);
    outshCol.push_back(tmp_shower);
  }

  return StatusCode::SUCCESS;
}

StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL3(  std::vector<const PandoraPlus::Calo1DCluster*>& barShowerUCol, 
                                                        std::vector<const PandoraPlus::Calo1DCluster*>& barShowerVCol, 
                                                        std::vector<PandoraPlus::Calo2DCluster*>& outshCol )
{
  outshCol.clear();

  const int NshowerU = barShowerUCol.size();
  const int NshowerV = barShowerVCol.size();

  double chi2[NshowerU][NshowerV];
  double chi2_E[NshowerU][NshowerV];
  double chi2_tx[NshowerU][NshowerV];
  double chi2_ty[NshowerU][NshowerV];

  double wi_E = settings.map_floatPars["chi2Wi_E"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);
  double wi_T = settings.map_floatPars["chi2Wi_T"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);

  TVector3 m_vec(0,0,0);
  double rotAngle = -(barShowerUCol[0]->getBars())[0]->getModule()*PI/4.;
  TVector3 Cblock( (barShowerUCol[0]->getBars())[0]->getPosition().x(),
                   (barShowerUCol[0]->getBars())[0]->getPosition().y(),
                   (barShowerVCol[0]->getBars())[0]->getPosition().z() );
  Cblock.RotateZ(rotAngle);
  map<double, pair<int, int> > m_chi2Map; m_chi2Map.clear();

  for(int ix=0;ix<NshowerU;ix++){
  for(int iy=0;iy<NshowerV;iy++){
    const PandoraPlus::Calo1DCluster* showerX = barShowerUCol[ix];
    const PandoraPlus::Calo1DCluster* showerY = barShowerVCol[iy];

    double Ex = showerX->getEnergy();
    double Ey = showerY->getEnergy();
    chi2_E[ix][iy] = pow(fabs(Ex-Ey)/settings.map_floatPars["sigmaE"], 2);
    double PosTx = C*(showerY->getT1()-showerY->getT2())/(2*settings.map_floatPars["nMat"]) + showerY->getPos().z();
    chi2_tx[ix][iy] = pow( fabs(PosTx-showerX->getPos().z())/settings.map_floatPars["sigmaPos"], 2 );

    double PosTy = C*(showerX->getT1()-showerX->getT2())/(2*settings.map_floatPars["nMat"]);
    m_vec.SetXYZ(showerY->getPos().x(), showerY->getPos().y(), showerY->getPos().z());
    m_vec.RotateZ(rotAngle);
    chi2_ty[ix][iy] = pow( fabs(PosTy - (m_vec-Cblock).x() )/settings.map_floatPars["sigmaPos"], 2);

    chi2[ix][iy] = chi2_E[ix][iy]*wi_E + (chi2_tx[ix][iy]+chi2_ty[ix][iy])*wi_T ;

    pair<int, int> p1(ix, iy);
    m_chi2Map[chi2[ix][iy]] = p1;    
  }}

  pair<int, int> lastpair;
  vector<pair<int, int>> indexVec; indexVec.clear();
  map<double, pair<int, int> >::iterator iter = m_chi2Map.begin();

  for(iter; iter!=m_chi2Map.end(); iter++){
    pair<int, int> indexpair = iter->second;
    bool inLine = false;
    bool isLast = false;
    for(int i=0; i<indexVec.size(); i++){
      if( indexVec.size() == min(NshowerU, NshowerV)-1 ) {lastpair = indexpair; isLast=true;  break; }
      if( indexpair.first == indexVec[i].first || indexpair.second == indexVec[i].second ) { inLine=true; break; }
    }
    if(isLast) break;
    if(inLine) continue;
    indexVec.push_back(indexpair);
  }
  if( indexVec.size()!=min(NshowerU, NshowerV)-1 )
    cout<<"ERROR in EnergyTimeMatchingAlg::GetMatchedShowersL3: found pair size "<<indexVec.size()<<" does not equal to min shower size -1 "<<min(NshowerU, NshowerV)-1<<endl;

  vector<const PandoraPlus::Calo1DCluster*> leftShowers; leftShowers.clear();
  for(int i=0; i<max(NshowerU, NshowerV); i++){
    bool fl_exist = false;
    for(int j=0; j<indexVec.size(); j++){
      int m_index = NshowerU>NshowerV ? indexVec[j].first : indexVec[j].second;
      if(i==m_index){ fl_exist = true; break; }
    }
    if(!fl_exist){
       const PandoraPlus::Calo1DCluster* m_shower = NshowerU>NshowerV ? barShowerUCol[i] : barShowerVCol[i] ;
       leftShowers.push_back(m_shower);
    }
  }
  if(leftShowers.size() != fabs( NshowerU-NshowerV )+1 )
    cout<<"ERROR in XYShowerChi2MatchingL1: Last pair number "<<leftShowers.size()<<" does not equal to shower difference "<<fabs( NshowerU-NshowerV )+1<<endl;

  for(int ip=0; ip<indexVec.size(); ip++){
    const PandoraPlus::Calo1DCluster* showerX = barShowerUCol[indexVec[ip].first];
    const PandoraPlus::Calo1DCluster* showerY = barShowerVCol[indexVec[ip].second];

    PandoraPlus::Calo2DCluster* tmp_shower = new PandoraPlus::Calo2DCluster();
    GetMatchedShowersL0(showerX, showerY, tmp_shower);
    outshCol.push_back(tmp_shower);
  }


  std::vector<PandoraPlus::Calo2DCluster*> m_showerinlayer; m_showerinlayer.clear();
  int ilast = NshowerU<NshowerV ? lastpair.first : lastpair.second;
  const PandoraPlus::Calo1DCluster* m_shower = NshowerU<NshowerV ? barShowerUCol[ilast] : barShowerVCol[ilast] ;
  GetMatchedShowersL1( m_shower, leftShowers, m_showerinlayer);

  outshCol.insert(outshCol.end(), m_showerinlayer.begin(), m_showerinlayer.end());


  return StatusCode::SUCCESS;
}


double** EnergyTimeMatchingAlg::GetClusterChi2Map( std::vector<std::vector<const PandoraPlus::Calo1DCluster*>>& barShowerUCol, 
                                                   std::vector<std::vector<const PandoraPlus::Calo1DCluster*>>& barShowerVCol )
{

//  add one function: if ShowerU[i] and ShowerV[j] both have cousin shower in another tower, set chi2[i][j] to 0. 

  const int NclusX = barShowerUCol.size(); 
  const int NclusY = barShowerVCol.size(); 

  if(NclusX==0 || NclusY==0) return nullptr; 
  //If one longidutinal cluster is empty in this layer: skip this layer. 
  for(int icx=0; icx<NclusX; icx++) if(barShowerUCol[icx].size()==0) return nullptr;
  for(int icy=0; icy<NclusY; icy++) if(barShowerVCol[icy].size()==0) return nullptr;

  double  **chi2map = new double*[NclusX];

  double chi2map_E[NclusX][NclusY]={0};
  double chi2map_tx[NclusX][NclusY]={0};
  double chi2map_ty[NclusX][NclusY]={0};

  double wi_E = settings.map_floatPars["chi2Wi_E"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);
  double wi_T = settings.map_floatPars["chi2Wi_T"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);


  TVector3 m_vec(0,0,0);
  double rotAngle = -(barShowerUCol[0][0]->getBars())[0]->getModule()*PI/4.;
  TVector3 Ctower((barShowerUCol[0][0]->getBars())[0]->getPosition().x(), (barShowerUCol[0][0]->getBars())[0]->getPosition().y(), (barShowerVCol[0][0]->getBars())[0]->getPosition().z());
  Ctower.RotateZ(rotAngle);

  for(int ix=0;ix<NclusX;ix++){
  chi2map[ix] = new double[NclusY];
  for(int iy=0;iy<NclusY;iy++){
    std::vector<const PandoraPlus::Calo1DCluster*> clusterU = barShowerUCol[ix];
    std::vector<const PandoraPlus::Calo1DCluster*> clusterV = barShowerVCol[iy];

    double min_chi2E = 999;
    double min_chi2tx = 999;
    double min_chi2ty = 999;

    for(int icx=0; icx<clusterU.size(); icx++){
    for(int icy=0; icy<clusterV.size(); icy++){
      const PandoraPlus::Calo1DCluster* showerX = clusterU[icx];      
      const PandoraPlus::Calo1DCluster* showerY = clusterV[icy];      


      double Ex = showerX->getEnergy();
      double Ey = showerY->getEnergy();
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
