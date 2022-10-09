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

  //Output: Calo2DClusters and Calo3DClusters. 
  std::vector<PandoraPlus::Calo2DCluster*> m_transhowerCol; m_transhowerCol.clear(); 
  std::vector<PandoraPlus::Calo3DCluster*> m_clusterCol; m_clusterCol.clear();

  //Readin: 3DClusters
  std::vector<PandoraPlus::Calo3DCluster*>* p_3DClusters = &(m_datacol.Cluster3DCol);

  //Loop for 3DClusters
  for(int i3d=0; i3d<p_3DClusters->size(); i3d++){
    //Basic unit: tower. 
    std::vector<const Calo3DCluster*> p_towerCol = p_3DClusters->at(i3d)->getTowers();

    //Loop for towers:  
    for(int it=0; it<p_towerCol.size(); it++){
 
      std::vector<const PandoraPlus::LongiCluster*> m_longiClUCol = p_towerCol.at(it)->getLongiClusterUCol("ESLongiCluster");
      std::vector<const PandoraPlus::LongiCluster*> m_longiClVCol = p_towerCol.at(it)->getLongiClusterVCol("ESLongiCluster");

      const int NclusX = m_longiClUCol.size(); 
      const int NclusY = m_longiClVCol.size(); 
      if(NclusX==0 || NclusY==0) continue;    
      std::vector<PandoraPlus::Calo3DCluster*> tmp_clusters; tmp_clusters.clear(); 
      
      //Case 2.1: 1*1
      if(NclusX==1 && NclusY==1){
        PandoraPlus::Calo3DCluster* tmp_clus = new PandoraPlus::Calo3DCluster();
        XYClusterMatchingL0(m_longiClUCol[0], m_longiClVCol[0], tmp_clus); 
        tmp_clusters.push_back(tmp_clus);
        //continue;
      }
      //Case 2.2: 1*N
      else if(NclusX==1){ XYClusterMatchingL1(m_longiClUCol[0], m_longiClVCol, tmp_clusters); }
      else if(NclusY==1){ XYClusterMatchingL1(m_longiClVCol[0], m_longiClUCol, tmp_clusters); }
      //Case 2.3: N*N
      else if( NclusX==NclusY ){ 
        XYClusterMatchingL2(m_longiClUCol, m_longiClVCol, tmp_clusters);
      }
      //Case 2.4
      else{ 
        XYClusterMatchingL3(m_longiClUCol, m_longiClVCol, tmp_clusters);
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


    //Merge reconstructed cluster
    //  Type1: from 1DShower aspect
    for(int ic=0; ic<m_clusterCol.size() && m_clusterCol.size()>1; ic++){
      for(int jc=ic+1; jc<m_clusterCol.size(); jc++){

        int NlinkedU = 0; 
        int NlinkedV = 0;

        //Save out BarShowers in Cluster[jc]
        std::vector<const Calo1DCluster*> m_1DshowerUCol; m_1DshowerUCol.clear();
        std::vector<const Calo1DCluster*> m_1DshowerVCol; m_1DshowerVCol.clear();
        for(int is=0; is<m_clusterCol[jc]->getCluster().size(); is++){
          for(int js=0; js<m_clusterCol[jc]->getCluster()[js]->getClusterU().size(); js++)
            m_1DshowerUCol.push_back( m_clusterCol[jc]->getCluster()[js]->getClusterU()[js] );
          for(int js=0; js<m_clusterCol[jc]->getCluster()[js]->getClusterV().size(); js++)
            m_1DshowerVCol.push_back( m_clusterCol[jc]->getCluster()[js]->getClusterV()[js] );
        }

        //Loop in Cluster[ic], count the linked showers
        for(int is=0; is<m_clusterCol[ic]->getCluster().size(); is++){
          const Calo2DCluster* p_shower = m_clusterCol[ic]->getCluster()[is];

          for(int icc=0; icc<p_shower->getClusterU()[0]->getCousinClusters().size(); icc++)
            if( find(m_1DshowerUCol.begin(), m_1DshowerUCol.end(), p_shower->getClusterU()[0]->getCousinClusters()[icc])!=m_1DshowerUCol.end() ){ NlinkedU++; break; }

          for(int icc=0; icc<p_shower->getClusterV()[0]->getCousinClusters().size(); icc++)
            if( find(m_1DshowerVCol.begin(), m_1DshowerVCol.end(), p_shower->getClusterV()[0]->getCousinClusters()[icc])!=m_1DshowerVCol.end() ){ NlinkedV++; break; }

          p_shower = nullptr;
        }

        if(NlinkedU/(float)m_clusterCol[ic]->getCluster().size()>0.7 || NlinkedV/(float)m_clusterCol[ic]->getCluster().size()>0.7){
          m_clusterCol[ic]->mergeCluster( m_clusterCol[jc] );
          m_clusterCol.erase(m_clusterCol.begin()+jc);
          jc--;
          if(jc<ic) jc=ic;
        }
    }}

    //  Type2: from longiCluster aspect
    for(int ic=0; ic<m_clusterCol.size() && m_clusterCol.size()>1; ic++){
      if( m_clusterCol[ic]->getLongiClusterUCol("LinkedLongiCluster")[0]->getCousinClusters().size()==0 && 
          m_clusterCol[ic]->getLongiClusterVCol("LinkedLongiCluster")[0]->getCousinClusters().size()==0 ) continue;
      for(int jc=ic+1; jc<m_clusterCol.size(); jc++){
        if( m_clusterCol[jc]->getLongiClusterUCol("LinkedLongiCluster")[0]->getCousinClusters().size()==0 && 
            m_clusterCol[jc]->getLongiClusterVCol("LinkedLongiCluster")[0]->getCousinClusters().size()==0 ) continue;

        std::vector<const PandoraPlus::LongiCluster*> m_cousinU = m_clusterCol[ic]->getLongiClusterUCol("LinkedLongiCluster")[0]->getCousinClusters(); 
        std::vector<const PandoraPlus::LongiCluster*> m_cousinV = m_clusterCol[ic]->getLongiClusterVCol("LinkedLongiCluster")[0]->getCousinClusters(); 
        if( find(m_cousinU.begin(), m_cousinU.end(), m_clusterCol[jc]->getLongiClusterUCol("LinkedLongiCluster")[0] )!=m_cousinU.end() ||
            find(m_cousinV.begin(), m_cousinV.end(), m_clusterCol[jc]->getLongiClusterVCol("LinkedLongiCluster")[0] )!=m_cousinV.end() ){
          m_clusterCol[ic]->mergeCluster( m_clusterCol[jc] );
          m_clusterCol.erase(m_clusterCol.begin()+jc);
          jc--;
          if(jc<ic) jc=ic;
        }        
    }}


  }//End loop 3DCluster. 


  m_datacol.map_ShowerInLayer["EcalShowerInLayer"] = m_transhowerCol;
  if(settings.map_boolPars["fl_WriteCluster"]) m_datacol.map_CaloCluster["EcalCluster"] = m_clusterCol;
  p_3DClusters = nullptr;

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}

//Longitudinal cluster: 1*1
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL0( const PandoraPlus::LongiCluster* m_longiClU, 
                                                       const PandoraPlus::LongiCluster* m_longiClV, 
                                                       PandoraPlus::Calo3DCluster* m_clus )
{
//cout<<"  Cluster matching for case: 1 * 1"<<endl;

  std::vector<int> layerindex; layerindex.clear(); 
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > map_showersXinlayer; map_showersXinlayer.clear();
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > map_showersYinlayer; map_showersYinlayer.clear();

  for(int is=0; is<m_longiClU->getBarShowers().size(); is++){
    int m_layer = m_longiClU->getBarShowers()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);  
    map_showersXinlayer[m_layer].push_back(m_longiClU->getBarShowers()[is]);
  }
  for(int is=0; is<m_longiClV->getBarShowers().size(); is++){
    int m_layer = m_longiClV->getBarShowers()[is]->getDlayer(); 
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
    map_showersYinlayer[m_layer].push_back(m_longiClV->getBarShowers()[is]);
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
    else GetFullMatchedShowers(m_showerXcol, m_showerYcol, m_showerinlayer);
    //else GetMatchedShowersL3(m_showerXcol, m_showerYcol, m_showerinlayer);
//cout<<"  After matching: shower size = "<<m_showerinlayer.size()<<endl;

    for(int is=0; is<m_showerinlayer.size(); is++) m_clus->addUnit(m_showerinlayer[is]);
  }
  m_clus->addLongiClusterU( "LinkedLongiCluster", m_longiClU );
  m_clus->addLongiClusterV( "LinkedLongiCluster", m_longiClV );

  return StatusCode::SUCCESS;
}


//Longitudinal cluster: 1*N
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL1( const PandoraPlus::LongiCluster* m_longiCl1, 
                                                       std::vector<const PandoraPlus::LongiCluster*>& m_longiClN, 
                                                       std::vector<PandoraPlus::Calo3DCluster*>& m_clusters )
{

cout<<"  Input LongiCluster: "<<endl;
for(int ic=0; ic<m_longiClN.size(); ic++){
for(int is=0; is<m_longiClN[ic]->getBarShowers().size(); is++){
printf("    (%.3f, %.3f, %.3f, %.3f) \t", m_longiClN[ic]->getBarShowers()[is]->getPos().X(), 
                                          m_longiClN[ic]->getBarShowers()[is]->getPos().Y(),
                                          m_longiClN[ic]->getBarShowers()[is]->getPos().Z(), 
                                          m_longiClN[ic]->getBarShowers()[is]->getEnergy() );
}
cout<<endl;
}

for(int is=0; is<m_longiCl1->getBarShowers().size(); is++){
printf("    (%.3f, %.3f, %.3f, %.3f) \t", m_longiCl1->getBarShowers()[is]->getPos().X(),
                                          m_longiCl1->getBarShowers()[is]->getPos().Y(),
                                          m_longiCl1->getBarShowers()[is]->getPos().Z(),
                                          m_longiCl1->getBarShowers()[is]->getEnergy() );
}
cout<<endl;


  m_clusters.clear(); 
  m_clusters.resize(m_longiClN.size());

  int slayer = m_longiCl1->getSlayer(); 
  std::vector<int> layerindex; layerindex.clear();
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > map_showersXinlayer; map_showersXinlayer.clear();
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > map_showersYinlayer; map_showersYinlayer.clear();

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
  
  std::map<int, std::vector<const PandoraPlus::Calo2DCluster*> > map_2Dshowersinlayer; map_2Dshowersinlayer.clear(); 

  sort(layerindex.begin(), layerindex.end());
  for(int il=0; il<layerindex.size(); il++){

    std::vector<const PandoraPlus::Calo1DCluster*> m_showerXcol = map_showersXinlayer[layerindex[il]];
    std::vector<const PandoraPlus::Calo1DCluster*> m_showerYcol = map_showersYinlayer[layerindex[il]];

cout<<"    XYClusterMatchingL1: In Layer #"<<layerindex[il]<<": (X, Y) = ("<<m_showerXcol.size()<<", "<<m_showerYcol.size()<<") "<<endl;

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
    else GetFullMatchedShowers(m_showerXcol, m_showerYcol, m_showerinlayer);
    //else GetMatchedShowersL3(m_showerXcol, m_showerYcol, m_showerinlayer);

    //map_2Dshowersinlayer[layerindex[il]] = m_showerinlayer; 
    for(int is=0; is<m_showerinlayer.size(); is++) map_2Dshowersinlayer[layerindex[il]].push_back( m_showerinlayer[is] );
  }


  for(int il=0; il<layerindex.size(); il++){
cout<<"    In Layer #"<<layerindex[il]<<":  Calo2DCluster size = "<<map_2Dshowersinlayer[layerindex[il]].size()<<endl;

  for(int is=0; is<map_2Dshowersinlayer[layerindex[il]].size(); is++){

//printf("      Calo2DCluster #%d: pos+E=(%.3f, %.3f, %.3f) \n", is, map_2Dshowersinlayer[layerindex[il]][is]->getPos().z(), map_2Dshowersinlayer[layerindex[il]][is]->getPos().y(), map_2Dshowersinlayer[layerindex[il]][is]->getShowerE() );

    const PandoraPlus::Calo1DCluster* m_barshower; 
    if(slayer==0) m_barshower=map_2Dshowersinlayer[layerindex[il]][is]->getClusterU()[0];
    else          m_barshower=map_2Dshowersinlayer[layerindex[il]][is]->getClusterV()[0];

//printf("      Selected barshower: DLayer=%d, sLayer=%d, pos+E=(%.3f, %.3f, %.3f) \n", m_barshower->getDlayer(), m_barshower->getSlayer(), m_barshower->getPos().Y(), m_barshower->getPos().Z(), m_barshower->getEnergy() );

    int index_longiclus=-1; 
    bool fl_foundsh = false; 
    for(int ic=0; ic<m_longiClN.size() && !fl_foundsh; ic++){
    for(int jc=0; jc<m_longiClN[ic]->getBarShowers().size() && !fl_foundsh; jc++){
      if(m_barshower==m_longiClN[ic]->getBarShowers()[jc]) {index_longiclus=ic; fl_foundsh=true; break; }
    }}
cout<<"      Found LongiShower: "<<index_longiclus<<endl;

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
      m_clusters[ic]->addLongiClusterU( "LinkedLongiCluster", m_longiCl1 );
      m_clusters[ic]->addLongiClusterV( "LinkedLongiCluster", m_longiClN[ic] );
    }
    else{          
      m_clusters[ic]->addLongiClusterU( "LinkedLongiCluster", m_longiClN[ic] );
      m_clusters[ic]->addLongiClusterV( "LinkedLongiCluster", m_longiCl1 );
    }
  }

  return StatusCode::SUCCESS;
}


//Longitudinal cluster: N*N
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL2( std::vector<const PandoraPlus::LongiCluster*>& m_longiClUCol, 
                                                       std::vector<const PandoraPlus::LongiCluster*>& m_longiClVCol, 
                                                       std::vector<PandoraPlus::Calo3DCluster*>& m_clusters  )
{
  if(m_longiClUCol.size()==0 || m_longiClVCol.size()==0 || m_longiClUCol.size()!=m_longiClVCol.size()) return StatusCode::SUCCESS;

  std::vector<int> layerindex; layerindex.clear();
  std::map<int, std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> > map_showersXinlayer; map_showersXinlayer.clear();
  std::map<int, std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> > map_showersYinlayer; map_showersYinlayer.clear();

  //Find layers need to match.
//cout<<"  XYClusterMatchingL2: Find layers need to match "<<endl;
  for(int ic=0; ic<m_longiClUCol.size(); ic++){
  for(int is=0; is<m_longiClUCol[ic]->getBarShowers().size(); is++){
    int m_layer = m_longiClUCol[ic]->getBarShowers()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
  }}
  for(int ic=0; ic<m_longiClVCol.size(); ic++){
  for(int is=0; is<m_longiClVCol[ic]->getBarShowers().size(); is++){
    int m_layer = m_longiClVCol[ic]->getBarShowers()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
  }}
  sort(layerindex.begin(), layerindex.end());

  //Fill shower maps
//cout<<"  XYClusterMatchingL2: Fill shower maps"<<endl;  
  for(int il=0; il<layerindex.size(); il++){
    for(int ic=0; ic<m_longiClUCol.size(); ic++)
      map_showersXinlayer[layerindex[il]].push_back( m_longiClUCol[ic]->getBarShowersInLayer(layerindex[il]) );
    for(int ic=0; ic<m_longiClVCol.size(); ic++)
      map_showersYinlayer[layerindex[il]].push_back( m_longiClVCol[ic]->getBarShowersInLayer(layerindex[il]) );
  }
//cout<<"  XYClusterMatchingL2: map size X "<<map_showersXinlayer.size()<<", Y "<<map_showersYinlayer.size()<<endl;

  //Get the chi2 map for N*N
//cout<<"  XYClusterMatchingL2: Get chi2 maps"<<endl;
  const int Nclus = m_longiClUCol.size();
  double  **map_chi2[14] = {NULL};  //TODO: hardcoded layer number here! 
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

      if(m_longiClUCol[ic]->getCousinClusters().size()!=0 && m_longiClUCol[jc]->getCousinClusters().size()!=0) sumchi2[ic][jc] = 0.;
      else sumchi2[ic][jc] += map_chi2[il][ic][jc];
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
    const PandoraPlus::LongiCluster* m_clusX = m_longiClUCol[Index[ii].first];
    const PandoraPlus::LongiCluster* m_clusY = m_longiClVCol[Index[ii].second];
    PandoraPlus::Calo3DCluster* tmp_clus = new PandoraPlus::Calo3DCluster();

    XYClusterMatchingL0(m_clusX, m_clusY, tmp_clus);
    m_clusters.push_back(tmp_clus);
  }

  return StatusCode::SUCCESS;
}


//Longitudinal cluster: M*N
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL3( std::vector<const PandoraPlus::LongiCluster*>& m_longiClUCol, 
                                                       std::vector<const PandoraPlus::LongiCluster*>& m_longiClVCol, 
                                                       std::vector<PandoraPlus::Calo3DCluster*>& m_clusters )
{



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

  int _module = barShowerU->getBars()[0]->getModule();
  int _dlayer = barShowerU->getBars()[0]->getDlayer();
  int _part   = barShowerU->getBars()[0]->getPart();
  int _stave  = barShowerU->getBars()[0]->getStave();


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
  outsh->addShowerU( barShowerU );
  outsh->addShowerV( barShowerV );
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
    else{ cout<<"ERROR: Input shower has no seed! Check! Use the most energitic bar as seed. "<<endl; 
      double m_maxE = -99;
      int index = -1; 
      for(int ib=0; ib<shower1->getBars().size(); ib++){
        if(shower1->getBars()[ib]->getEnergy()>m_maxE) { m_maxE=shower1->getBars()[ib]->getEnergy(); index=ib; }
      }
      if(index>0) m_wiseed = shower1->getBars()[index]->Clone(); 
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
  TVector3 Cblock((barShowerUCol[0]->getBars())[0]->getPosition().x(), (barShowerUCol[0]->getBars())[0]->getPosition().y(), (barShowerVCol[0]->getBars())[0]->getPosition().z());
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
