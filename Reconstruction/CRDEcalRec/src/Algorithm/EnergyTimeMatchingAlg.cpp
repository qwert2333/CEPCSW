#ifndef _ENERGYTIMEMATCHING_ALG_C
#define _ENERGYTIMEMATCHING_ALG_C

#include "Algorithm/EnergyTimeMatchingAlg.h"

void EnergyTimeMatchingAlg::Settings::SetInitialValue(){
  chi2Wi_E = 1;
  chi2Wi_T = 10;
  th_chi2  = -1;
  sigmaE   = 0.10;
  fl_UseChi2 = true;
  Debug = 0;
}


EnergyTimeMatchingAlg::EnergyTimeMatchingAlg(){

}

StatusCode Initialize(){
	return StatusCode::SUCCESS;
}

StatusCode EnergyTimeMatchingAlg::RunAlgorithm( EnergyTimeMatchingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol ){
  settings = m_settings;

  std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_2DshowerCol; m_2DshowerCol.clear();
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_ClusterCol;  m_ClusterCol.clear();  

  //Basic unit: tower. 
  std::vector<CRDEcalEDM::CRDCaloTower> m_towerCol = m_datacol.TowerCol;
  if(m_towerCol.size()==0){ std::cout<<"Warning: Empty input in EnergyTimeMatchingAlg! Please check previous algorithm!"<<std::endl; return StatusCode::SUCCESS; }

  for(int it=0; it<m_towerCol.size(); it++){
    std::vector<CRDEcalEDM::CRDCaloBlock> m_blocks = m_towerCol[it].getBlocks(); 

    //Case1: Don't use shadow cluster or chi2 matching, Save all combinations for 3D clustering.
    if(!settings.fl_UseChi2){
      for(int ib=0;ib<m_blocks.size();ib++){
        std::vector<CRDEcalEDM::CRDCaloBarShower> showerXCol = m_blocks[ib].getShowerXCol();
        std::vector<CRDEcalEDM::CRDCaloBarShower> showerYCol = m_blocks[ib].getShowerYCol();
   
        std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_showerinlayer; m_showerinlayer.clear();
        CRDEcalEDM::CRDCaloHit3DCluster tmp_clus; tmp_clus.Clear();

        GetFullMatchedShowers( showerXCol, showerYCol, m_showerinlayer );
        tmp_clus.setShowerVec(m_showerinlayer);
        
        m_2DshowerCol.insert( m_2DshowerCol.end(), m_showerinlayer.begin(), m_showerinlayer.end() );   
        m_ClusterCol.push_back(tmp_clus);
      }
    }
    

    //Case2: chi2 matching (combined chi2 in all layers)
    else{
      std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_longiClXCol = m_towerCol[it].getLongiClusterXCol();
      std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_longiClYCol = m_towerCol[it].getLongiClusterYCol();

      const int NclusX = m_longiClXCol.size(); 
      const int NclusY = m_longiClYCol.size(); 
      if(NclusX==0 || NclusY==0) continue; 

      std::vector<CRDEcalEDM::CRDCaloHit3DCluster> tmp_clusters; tmp_clusters.clear(); 
      CRDEcalEDM::CRDCaloHit3DCluster tmp_clus; tmp_clus.Clear(); 
      
      //Case 2.1
      if(NclusX==1 && NclusY==1){
        chi2=0; sumchi2=0; 
        XYClusterMatchingL0(m_longiClXCol[0], m_longiClYCol[0], tmp_clus); 
        tmp_clusters.push_back(m_3Dclus);
        continue;
      }
      //Case 2.2
      else if(NclusX==1){ XYClusterMatchingL1(m_longiClXCol[0], m_longiClYCol, tmp_clusters); continue; }
      else if(NclusY==1){ XYClusterMatchingL1(m_longiClYCol[0], m_longiClXCol, tmp_clusters); continue; }

      //Case 2.3
      else if( NclusX==NclusY ){ 
        XYClusterMatchingL2(m_longiClXCol, m_longiClYCol, tmp_clusters);
      }
      //Case 2.4
      else{ 
        XYClusterMatchingL3(m_longiClXCol, m_longiClYCol, tmp_clusters);
      }

      m_ClusterCol.insert(m_ClusterCol.end(), tmp_clusters.begin(), tmp_clusters.end());
      for(int ic=0; ic<tmp_clusters.size(); ic++){ 
        std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_showersinclus = tmp_clusters[ic].get2DShowers(); 
        m_2DshowerCol.insert(m_2DshowerCol.end(), m_showersinclus.begin(), m_showersinclus.end());
      }

    }


  }
  m_datacol.Shower2DCol = m_2DshowerCol;
  m_datacol.Clus3DCol = m_ClusterCol; 


  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}

//Longitudinal cluster: 1*1
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL0( CRDEcalEDM::CRDCaloHitLongiCluster& m_longiClX, 
                                                      CRDEcalEDM::CRDCaloHitLongiCluster& m_longiClY, 
                                                      CRDEcalEDM::CRDCaloHit3DCluster& m_clus )
{
  m_clus.Clear(); 

  std::vector<int> layerindex; layerindex.clear(); 
  std::map<int, std::vector<CRDEcalEDM::CRDCaloBarShower> > map_showersXinlayer; map_showersXinlayer.clear();
  std::map<int, std::vector<CRDEcalEDM::CRDCaloBarShower> > map_showersYinlayer; map_showersYinlayer.clear();
  for(int is=0; is<m_longiClX.getBarShowers().size(); is++){
    int m_layer = m_longiClX.getBarShowers()[is].getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);  
    map_showersXinlayer[m_layer].push_back(m_longiClX.getBarShowers()[is]);
  }
  for(int is=0; is<m_longiClY.getBarShowers().size(); is++){
    int m_layer = m_longiClY.getBarShowers()[is].getDlayer(); 
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
    map_showersYinlayer[m_layer].push_back(m_longiClY.getBarShowers()[is]);
  }

  for(int il=0; il<layerindex.size(); il++){
    std::vector<CRDEcalEDM::CRDCaloBarShower> m_showerXcol = map_showersXinlayer[layerindex[il]];
    std::vector<CRDEcalEDM::CRDCaloBarShower> m_showerYcol = map_showersYinlayer[layerindex[il]];

    std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_showerinlayer; m_showerinlayer.clear();
    CRDEcalEDM::CRDCaloHitTransShower tmp_shower; tmp_shower.Clear();

    if(m_showerXcol.size()==0 || m_showerYcol.size()==0) continue;
    else if(m_showerXcol.size()==1 && m_showerYcol.size()==1){ GetMatchedShowersL0(m_showerXcol[0], m_showerYcol[0], tmp_shower); m_showerinlayer.push_back(tmp_shower); }
    else if(m_showerXcol.size()==1) GetMatchedShowersL1(m_showerXcol[0], m_showerYcol, m_showerinlayer );
    else if(m_showerYcol.size()==1) GetMatchedShowersL1(m_showerYcol[0], m_showerXcol, m_showerinlayer );
    else GetFullMatchedShowers(m_showerXcol, m_showerYcol, m_showerinlayer);
    //else if(m_showerXcol.size()== m_showerYcol.size()) GetMatchedShowersL2(m_showerXcol, m_showerYcol, m_showerinlayer );
    //else GetMatchedShowersL3(m_showerXcol, m_showerYcol, m_showerinlayer);

    for(int is=0; is<m_showerinlayer.size(); is++) m_clus.AddShower(m_showerinlayer[is]);
  }

  return StatusCode::SUCCESS;
}


//Longitudinal cluster: 1*N
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL1( CRDEcalEDM::CRDCaloHitLongiCluster& m_longiCl1, 
                                                       std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClN, 
                                                       std::vector<CRDEcalEDM::CRDCaloHit3DCluster>& m_clusters )
{
  m_clusters.clear(); 
  m_clusters.resize(m_longiClN.size());

  int slayer = m_longiCl1.getSlayer(); 
  std::vector<int> layerindex; layerindex.clear();
  std::map<int, std::vector<CRDEcalEDM::CRDCaloBarShower> > map_showersXinlayer; map_showersXinlayer.clear();
  std::map<int, std::vector<CRDEcalEDM::CRDCaloBarShower> > map_showersYinlayer; map_showersYinlayer.clear();

  for(int is=0; is<m_longiCl1.getBarShowers().size(); is++){
    int m_layer = m_longiCl1.getBarShowers()[is].getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);

    if(slayer==0) map_showersXinlayer[m_layer].push_back(m_longiCl1.getBarShowers()[is]);
    else map_showersYinlayer[m_layer].push_back(m_longiCl1.getBarShowers()[is]);
  }
  for(int ic=0; ic<m_longiClN.size(); ic++){
  for(int is=0; is<m_longiClN[ic].getBarShowers().size(); is++){
    int m_layer = m_longiClN[ic].getBarShowers()[is].getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);

    if(slayer==0) map_showersYinlayer[m_layer].push_back(m_longiClN[ic].getBarShowers()[is]);
    else map_showersXinlayer[m_layer].push_back(m_longiClN[ic].getBarShowers()[is]);
  }}
  
  std::map<int, std::vector<CRDEcalEDM::CRDCaloHitTransShower> > map_2Dshowersinlayer; map_2Dshowersinlayer.clear(); 
  for(int il=0; il<layerindex.size(); il++){
    std::vector<CRDEcalEDM::CRDCaloBarShower> m_showerXcol = map_showersXinlayer[layerindex[il]];
    std::vector<CRDEcalEDM::CRDCaloBarShower> m_showerYcol = map_showersYinlayer[layerindex[il]];

    std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_showerinlayer; m_showerinlayer.clear();
    CRDEcalEDM::CRDCaloHitTransShower tmp_shower; tmp_shower.Clear();

    if(m_showerXcol.size()==0 || m_showerYcol.size()==0) continue;
    else if(m_showerXcol.size()==1 && m_showerYcol.size()==1){ GetMatchedShowersL0(m_showerXcol[0], m_showerYcol[0], tmp_shower); m_showerinlayer.push_back(tmp_shower); }
    else if(m_showerXcol.size()==1) GetMatchedShowersL1(m_showerXcol[0], m_showerYcol, m_showerinlayer );
    else if(m_showerYcol.size()==1) GetMatchedShowersL1(m_showerYcol[0], m_showerXcol, m_showerinlayer );
    else GetFullMatchedShowers(m_showerXcol, m_showerYcol, m_showerinlayer);
    //else if(m_showerXcol.size()== m_showerYcol.size()) GetMatchedShowersL2(m_showerXcol, m_showerYcol, m_showerinlayer );
    //else GetMatchedShowersL3(m_showerXcol, m_showerYcol, m_showerinlayer);

    map_2Dshowersinlayer[layerindex[il]] = m_showerinlayer; 
  }

  for(int il=0; il<layerindex.size(); il++){
  for(int is=0; is<map_2Dshowersinlayer[layerindex[il]].size(); is++){
    CRDEcalEDM::CRDCaloBarShower m_barshower; 
    if(slayer==0) m_barshower=map_2Dshowersinlayer[layerindex[il]][is].getShowerY();
    else m_barshower=map_2Dshowersinlayer[layerindex[il]][is].getShowerX();

    int index_longiclus=-1; 
    bool fl_foundsh = false; 
    for(int ic=0; ic<m_longiClN.size() && !fl_foundsh; ic++){
    for(int jc=0; jc<m_longiClN[ic].getBarShowers().size() && !fl_foundsh; jc++){
      if(m_barshower==m_longiClN[ic].getBarShowers()[jc]){index_longiclus=ic; fl_foundsh=true; break; }
    }}

    if(index_longiclus<0){ std::cout<<"WARNING: did not find properate longitudinal cluster! "<<std::endl; continue; }
    m_clusters[index_longiclus].AddShower(map_2Dshowersinlayer[layerindex[il]][is]);
  }}

  return StatusCode::SUCCESS;
}


//Longitudinal cluster: N*N
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL2( std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClXCol, 
                                                       std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClYCol, 
                                                       std::vector<CRDEcalEDM::CRDCaloHit3DCluster>& m_clusters  )
{
  if(m_longiClXCol.size()==0 || m_longiClYCol.size()==0 || m_longiClXCol.size()!=m_longiClYCol.size()) return StatusCode::SUCCESS;

  std::vector<int> layerindex; layerindex.clear();
  std::map<int, std::vector<std::vector<CRDEcalEDM::CRDCaloBarShower>> > map_showersXinlayer; map_showersXinlayer.clear();
  std::map<int, std::vector<std::vector<CRDEcalEDM::CRDCaloBarShower>> > map_showersYinlayer; map_showersYinlayer.clear();

  for(int ic=0; ic<m_longiClXCol.size(); ic++){
  for(int is=0; is<m_longiClXCol[ic].getBarShowers().size(); is++){
    int m_layer = m_longiClXCol[ic].getBarShowers()[is].getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
  }}
  for(int ic=0; ic<m_longiClYCol.size(); ic++){
  for(int is=0; is<m_longiClYCol[ic].getBarShowers().size(); is++){
    int m_layer = m_longiClYCol[ic].getBarShowers()[is].getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
  }}

  for(int il=0; il<layerindex.size(); il++){
    for(int ic=0; ic<m_longiClXCol.size(); ic++)
      map_showersXinlayer[layerindex[il]].push_back( m_longiClXCol[ic].getBarShowersInLayer(layerindex[il]) );
    for(int ic=0; ic<m_longiClYCol.size(); ic++)
      map_showersYinlayer[layerindex[il]].push_back( m_longiClYCol[ic].getBarShowersInLayer(layerindex[il]) );
  }

  //Get the chi2 map for N*N
  const int Nclus = m_longiClXCol.size();
  double map_chi2[14][Nclus][Nclus] = {0};
  double sumchi2[Nclus][Nclus] = {0};
  std::map<int, std::vector<CRDEcalEDM::CRDCaloHitTransShower> > map_2Dshowersinlayer; map_2Dshowersinlayer.clear();  

  for(int il=0; il<layerindex.size(); il++){
    std::vector<std::vector<CRDEcalEDM::CRDCaloBarShower>> m_showerXcol = map_showersXinlayer[layerindex[il]];
    std::vector<std::vector<CRDEcalEDM::CRDCaloBarShower>> m_showerYcol = map_showersYinlayer[layerindex[il]];

    map_chi2[layerindex[il]] = GetClusterChi2Map(m_showerXcol, m_showerYcol);
  }

  for(int ic=0; ic<Nclus; ic++){
  for(int jc=0; jc<Nclus; jc++){
    for(int il=0; il<14; il++) sumchi2[ic][jc] += map_chi2[il][ic][jc];
  }}


  //Get the chi2 of N! combinations
  int Ncomb=1;
  for(int i=Nclus; i>0; i--) Ncomb = Ncomb*i;

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
    for(int i=0;i<Index.size();i++) chi2_tot += sumchi2[Index[i].first][Index[i].second];
    matchingMap[chi2_tot] = Index;

    Index.clear();
    if(!next_permutation(num, num+Nshower)) break;
  }

  //map is ordered with [double] value, first element has the smallest chi2 value. 
  map<double, vector<pair<int, int>> >::iterator iter = matchingMap.begin();
  vector<pair<int, int>> Index = iter->second;

  for(int ii=0; ii<Index.size(); ii++){
    CRDEcalEDM::CRDCaloHitLongiCluster m_clusX = m_longiClXCol[Index[ii].first];
    CRDEcalEDM::CRDCaloHitLongiCluster m_clusY = m_longiClYCol[Index[ii].first];
    CRDEcalEDM::CRDCaloHit3DCluster tmp_clus; tmp_clus.Clear();

    XYClusterMatchingL0(m_clusX, m_clusY, tmp_clus);
    m_clusters.push_back(tmp_clus);
  }

  return StatusCode::SUCCESS;
}


//Longitudinal cluster: M*N
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL3( std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClXCol, 
                                                       std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClYCol, 
                                                       std::vector<CRDEcalEDM::CRDCaloHit3DCluster>& m_clusters )
{



  return StatusCode::SUCCESS;
}



StatusCode EnergyTimeMatchingAlg::GetFullMatchedShowers(
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol,
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol,
           std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol )
{
  outshCol.clear();

  for(int is=0;is<barShowerXCol.size();is++){
  for(int js=0;js<barShowerYCol.size();js++){
    if(barShowerXCol[is].getBars().size()==0 || barShowerYCol[js].getBars().size()==0) continue;
    CRDEcalEDM::CRDCaloHitTransShower tmp_shower; tmp_shower.Clear();
    GetMatchedShowersL0(barShowerXCol[is], barShowerYCol[js], tmp_shower);
    outshCol.push_back(tmp_shower);
  }}

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL0( CRDEcalEDM::CRDCaloBarShower& barShowerX,   
                                                       CRDEcalEDM::CRDCaloBarShower& barShowerY,
                                                       CRDEcalEDM::CRDCaloHitTransShower& outsh )
{

  std::vector<edm4hep::ConstCalorimeterHit> m_digiCol; m_digiCol.clear();
  int NbarsX = barShowerX.getBars().size();
  int NbarsY = barShowerY.getBars().size();
  if(NbarsX==0 || NbarsY==0){ std::cout<<"WARNING: empty DigiHitsCol returned!"<<std::endl; return StatusCode::SUCCESS; }

  int _module = barShowerX.getBars()[0].getModule();
  int _dlayer = barShowerX.getBars()[0].getDlayer();
  int _part   = barShowerX.getBars()[0].getPart();
  int _stave  = barShowerX.getBars()[0].getStave();

  float rotAngle = -_module*PI/4.;

  TVector3 m_vec(0,0,0);
  for(int ibar=0;ibar<NbarsX;ibar++){
    CRDEcalEDM::CRDCaloBar barx = barShowerX.getBars()[ibar];

    m_vec.SetXYZ(barx.getPosition().x(), barx.getPosition().y(), barx.getPosition().z());
    m_vec.RotateZ(rotAngle);
    barx.setPosition( m_vec );

    for(int jbar=0;jbar<NbarsY;jbar++){
      CRDEcalEDM::CRDCaloBar bary = barShowerY.getBars()[jbar];
      m_vec.SetXYZ(bary.getPosition().x(), bary.getPosition().y(), bary.getPosition().z());
      m_vec.RotateZ(rotAngle);
      bary.setPosition( m_vec );

      TVector3 p_hit(bary.getPosition().x(), (barx.getPosition().y()+bary.getPosition().y())/2., barx.getPosition().z() );
      p_hit.RotateZ(-rotAngle);
      edm4hep::Vector3f m_vec3f(p_hit.x(), p_hit.y(), p_hit.z());
      float m_En = barx.getEnergy()*bary.getEnergy()/barShowerY.getE() + barx.getEnergy()*bary.getEnergy()/barShowerX.getE();
      edm4hep::CalorimeterHit hit;
      hit.setCellID(0);
      hit.setPosition(m_vec3f);
      hit.setEnergy(m_En);
      m_digiCol.push_back(hit);
    }
  }

  outsh.setBarShowers( barShowerX, barShowerY );
  outsh.setCaloHits( m_digiCol );
  outsh.setIDInfo( _module, _stave, _dlayer, _part);

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL1( CRDEcalEDM::CRDCaloBarShower& shower1,
                                                       std::vector<CRDEcalEDM::CRDCaloBarShower>& showerNCol, 
                                                       std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol )
{

  outshCol.clear();

  int _slayer = shower1.getBars()[0].getSlayer();

  const int NshY = showerNCol.size();
  double totE_shY = 0;
  double EshY[NshY] = {0};
  for(int is=0;is<NshY;is++){ EshY[is] = showerNCol[is].getE(); totE_shY += EshY[is]; }
  for(int is=0;is<NshY;is++){
    double wi_E = EshY[is]/totE_shY;
    CRDEcalEDM::CRDCaloBarShower m_splitshower1;

    CRDEcalEDM::CRDCaloBar m_wiseed = shower1.getSeed();
    m_wiseed.setQ( wi_E*m_wiseed.getQ1(), wi_E*m_wiseed.getQ2() );

    std::vector<CRDEcalEDM::CRDCaloBar> m_wibars;
    for(int ib=0;ib<shower1.getBars().size();ib++){
      CRDEcalEDM::CRDCaloBar m_wibar = shower1.getBars()[ib];
      m_wibar.setQ(wi_E*m_wibar.getQ1(), wi_E*m_wibar.getQ2());
      m_wibars.push_back(m_wibar);
    }

    m_splitshower1.setBars( m_wibars );
    m_splitshower1.setSeed( m_wiseed );

    CRDEcalEDM::CRDCaloHitTransShower m_shower; m_shower.Clear();
    if(_slayer==0 ) GetMatchedShowersL0( m_splitshower1, showerNCol[is], m_shower);
    else            GetMatchedShowersL0( showerNCol[is], m_splitshower1, m_shower);
    outshCol.push_back( m_shower );
  }

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL2( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol,
                                                       std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol,
                                                       std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol )
{
  outshCol.clear();
  if(barShowerXCol.size() != barShowerYCol.size() ) return StatusCode::FAILURE;

  const int Nshower = barShowerXCol.size();
  double chi2[Nshower][Nshower];
  double chi2_E[Nshower][Nshower];
  double chi2_tx[Nshower][Nshower];
  double chi2_ty[Nshower][Nshower];

  double wi_E = settings.chi2Wi_E/(settings.chi2Wi_E + settings.chi2Wi_T);
  double wi_T = settings.chi2Wi_T/(settings.chi2Wi_E + settings.chi2Wi_T);

  TVector3 m_vec(0,0,0);
  double rotAngle = -(barShowerXCol[0].getBars())[0].getModule()*PI/4.;
  TVector3 Cblock((barShowerXCol[0].getBars())[0].getPosition().x(), (barShowerXCol[0].getBars())[0].getPosition().y(), (barShowerYCol[0].getBars())[0].getPosition().z());
  Cblock.RotateZ(rotAngle);

  for(int ix=0;ix<Nshower;ix++){
  for(int iy=0;iy<Nshower;iy++){
    CRDEcalEDM::CRDCaloBarShower showerX = barShowerXCol[ix];
    CRDEcalEDM::CRDCaloBarShower showerY = barShowerYCol[iy];

    double Ex = showerX.getE();
    double Ey = showerY.getE();
    chi2_E[ix][iy] = pow(fabs(Ex-Ey)/settings.sigmaE, 2);
    double PosTx = C*(showerY.getT1()-showerY.getT2())/(2*settings.nMat) + showerY.getPos().z();
    chi2_tx[ix][iy] = pow( fabs(PosTx-showerX.getPos().z())/settings.sigmaPos, 2 );

    double PosTy = C*(showerX.getT1()-showerX.getT2())/(2*settings.nMat);
    m_vec.SetXYZ(showerY.getPos().x(), showerY.getPos().y(), showerY.getPos().z());
    m_vec.RotateZ(rotAngle);
    chi2_ty[ix][iy] = pow( fabs(PosTy - (m_vec-Cblock).x() )/settings.sigmaPos, 2);

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
    CRDEcalEDM::CRDCaloBarShower showerX = barShowerXCol[Index[i].first];
    CRDEcalEDM::CRDCaloBarShower showerY = barShowerYCol[Index[i].second];

    CRDEcalEDM::CRDCaloHitTransShower tmp_shower; tmp_shower.Clear();
    GetMatchedShowersL0(showerX, showerY, tmp_shower);
    outshCol.push_back(tmp_shower);
  }

  return StatusCode::SUCCESS;
}


double** EnergyTimeMatchingAlg::GetClusterChi2Map( std::vector<std::vector<CRDEcalEDM::CRDCaloBarShower>>& barShowerXCol, 
                                                   std::vector<std::vector<CRDEcalEDM::CRDCaloBarShower>>& barShowerYCol )
{
  if(NclusX==0 || NclusY==0) return nullptr; 

  const int NclusX = barShowerXCol.size(); 
  const int NclusY = barShowerYCol.size(); 
  double chi2map[NclusX][NclusY]={0};
  double chi2map_E[NclusX][NclusY]={0};
  double chi2map_tx[NclusX][NclusY]={0};
  double chi2map_ty[NclusX][NclusY]={0};

  double wi_E = settings.chi2Wi_E/(settings.chi2Wi_E + settings.chi2Wi_T);
  double wi_T = settings.chi2Wi_T/(settings.chi2Wi_E + settings.chi2Wi_T);

  TVector3 m_vec(0,0,0);
  double rotAngle = -(barShowerXCol[0][0].getBars())[0].getModule()*PI/4.;
  TVector3 Ctower((barShowerXCol[0][0].getBars())[0].getPosition().x(), (barShowerXCol[0][0].getBars())[0].getPosition().y(), (barShowerYCol[0][0].getBars())[0].getPosition().z());
  Ctower.RotateZ(rotAngle);

  for(int ix=0;ix<NclusX;ix++){
  for(int iy=0;iy<NclusY;iy++){
    std::vector<CRDEcalEDM::CRDCaloBarShower> clusterX = barShowerXCol[ix];
    std::vector<CRDEcalEDM::CRDCaloBarShower> clusterY = barShowerYCol[iy];

    double min_chi2E = 999;
    double min_chi2tx = 999;
    double min_chi2ty = 999;

    for(int icx=0; icx<clusterX.size(); icx++){
    for(int icy=0; icy<clusterY.size(); icy++){
      CRDEcalEDM::CRDCaloBarShower showerX = clusterX[icx];      
      CRDEcalEDM::CRDCaloBarShower showerY = clusterY[icy];      


      double Ex = showerX.getE();
      double Ey = showerY.getE();
      double chi2_E = pow(fabs(Ex-Ey)/settings.sigmaE, 2);
      double PosTx = C*(showerY.getT1()-showerY.getT2())/(2*settings.nMat) + showerY.getPos().z();
      double chi2_tx = pow( fabs(PosTx-showerX.getPos().z())/settings.sigmaPos, 2 );
   
      double PosTy = C*(showerX.getT1()-showerX.getT2())/(2*settings.nMat);
      m_vec.SetXYZ(showerY.getPos().x(), showerY.getPos().y(), showerY.getPos().z());
      m_vec.RotateZ(rotAngle);
      double chi2_ty = pow( fabs(PosTy - (m_vec-Ctower).x() )/settings.sigmaPos, 2);

      if(chi2_E<min_chi2E) min_chi2E=chi2_E; 
      if(chi2_tx<min_chi2tx) min_chi2tx=chi2_tx; 
      if(chi2_ty<min_chi2ty) min_chi2ty=chi2_ty;
    }}

    chi2map_E[ix][iy] = min_chi2E; 
    chi2map_tx[ix][iy] = min_chi2tx;
    chi2map_ty[ix][iy] = min_chi2ty;
    chi2map[ix][iy] = chi2_E[ix][iy]*wi_E + (chi2_tx[ix][iy]+chi2_ty[ix][iy])*wi_T ;

  }}

  return chi2map; 

}


#endif
