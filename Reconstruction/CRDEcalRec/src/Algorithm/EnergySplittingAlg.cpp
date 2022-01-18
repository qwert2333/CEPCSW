#ifndef _ENERGYSPLITTING_ALG_C
#define _ENERGYSPLITTING_ALG_C

#include "Algorithm/EnergySplittingAlg.h"
void EnergySplittingAlg::Settings::SetInitialValue(){

  th_GoodLayer = 4;
  th_split = -1; 
  Debug = 0;
}

StatusCode EnergySplittingAlg::Initialize(){

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::RunAlgorithm( EnergySplittingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){
  settings = m_settings;

  for(int it=0; it<m_datacol.TowerCol.size(); it++){
    std::vector<CRDEcalEDM::CRDCaloBlock> m_blocksInTower = m_datacol.TowerCol[it].getBlocks(); 
    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_clusXCol = m_datacol.TowerCol[it].getLongiClusterXCol();
    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_clusYCol = m_datacol.TowerCol[it].getLongiClusterYCol();


    //Identify longitudinal cluster quality
    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_goodClusXCol; m_goodClusXCol.clear();
    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_badClusXCol; m_badClusXCol.clear();
    for(int ic=0; ic<m_clusXCol.size(); ic++){
      if( m_clusXCol[ic].getEndDlayer()-m_clusXCol[ic].getBeginningDlayer()>=settings.th_GoodLayer ) m_goodClusXCol.push_back(m_clusXCol[ic]);
      else m_badClusXCol.push_back(m_clusXCol[ic]);
    }

    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_goodClusYCol; m_goodClusYCol.clear();
    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_badClusYCol; m_badClusYCol.clear();
    for(int ic=0; ic<m_clusYCol.size(); ic++){
      if( m_clusYCol[ic].getEndDlayer()-m_clusYCol[ic].getBeginningDlayer()>=settings.th_GoodLayer ) m_goodClusYCol.push_back(m_clusYCol[ic]);
      else m_badClusYCol.push_back(m_clusYCol[ic]);
    }
    

    //Make clusters and shwoers in blocks
    for(int ib=0; ib<m_blocksInTower.size(); ib++ ){
      m_blocksInTower[ib].ClearShower();
      ClusteringInBlock(m_blocksInTower[ib], m_goodClusXCol, m_goodClusYCol); 
    }

    std::vector<CRDCaloHitLongiCluster> m_longiClusXCol; m_longiClusXCol.clear();
    std::vector<CRDCaloHitLongiCluster> m_longiClusYCol; m_longiClusYCol.clear();
    LongiClusterLinking(m_blocksInTower, m_clusXCol, m_longiClusXCol);
    LongiClusterLinking(m_blocksInTower, m_clusYCol, m_longiClusYCol);

    m_datacol.TowerCol[it].SetLongiClusters( m_longiClusXCol, m_longiClusYCol );
    m_datacol.TowerCol[it].SetBlocks( m_blocksInTower );
  }

  return StatusCode::SUCCESS;
}

StatusCode EnergySplittingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}

StatusCode EnergySplittingAlg::ClusteringInBlock( CRDEcalEDM::CRDCaloBlock& m_block, 
                                                  std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClusXCol,
                                                  std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClusYCol  )
{

  std::vector<CRDEcalEDM::CRDCaloBarCluster> m_barClusXCol; m_barClusXCol.clear();
  std::vector<CRDEcalEDM::CRDCaloBarCluster> m_barClusYCol; m_barClusYCol.clear();
  
  std::vector<CRDEcalEDM::CRDCaloBar> m_barColX = m_block.getBarXCol(); 
  std::vector<CRDEcalEDM::CRDCaloBar> m_barColY = m_block.getBarYCol(); 
  Clustering( m_barColX, m_barClusXCol, m_longiClusXCol );
  Clustering( m_barColY, m_barClusYCol, m_longiClusYCol );

  m_block.setClusterXCol(m_barClusXCol);
  m_block.setClusterYCol(m_barClusYCol);

  std::vector<CRDEcalEDM::CRDCaloBarShower> m_barShowerXCol; m_barShowerXCol.clear();
  for(int ic=0; ic<m_barClusXCol.size(); ic++){
    std::vector<CRDEcalEDM::CRDCaloBarShower> m_showers; m_showers.clear();
    ClusterSplitting( m_barClusXCol[ic], m_showers );

    if(m_showers.size()==0) continue;
    m_barShowerXCol.insert( m_barShowerXCol.end(), m_showers.begin(), m_showers.end() );
  }

  std::vector<CRDEcalEDM::CRDCaloBarShower> m_barShowerYCol; m_barShowerYCol.clear();
  for(int ic=0; ic<m_barClusYCol.size(); ic++){
    std::vector<CRDEcalEDM::CRDCaloBarShower> m_showers; m_showers.clear();
    ClusterSplitting( m_barClusYCol[ic], m_showers );

    if(m_showers.size()==0) continue;
    m_barShowerYCol.insert( m_barShowerYCol.end(), m_showers.begin(), m_showers.end() );
  }
  m_block.setShowerXCol( m_barShowerXCol );
  m_block.setShowerYCol( m_barShowerYCol );

  return StatusCode::SUCCESS;
}


StatusCode LongiClusterLinking( std::vector<CRDEcalEDM::CRDCaloBlock>& m_blocks, 
                                std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_oldClusCol, 
                                std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_outClusCol )
{
  if(m_blocks.size()==0 || m_oldClusCol.size()==0) return StatusCode::SUCCESS;
  m_outClusCol.clear();

  bool fl_isXclus = (m_oldClusCol[0].getSlayer()==0);
  std::vector<CRDEcalEDM::CRDCaloBarShower> m_showers; m_showers.clear(); 
  for(int ib=0; ib<m_blocks.size(); ib++){
    std::vector<CRDEcalEDM::CRDCaloBarShower> tmp_showersinblock; tmp_showersinblock.clear(); 
    if(fl_isXclus) tmp_showersinblock = m_blocks[ib].getShowerXCol();
    else tmp_showersinblock = m_blocks[ib].getShowerYCol();
    m_showers.insert(m_showers.end(), tmp_showersinblock.begin(), tmp_showersinblock.end());
  }

  for(int ic=0; ic<m_oldClusCol.size(); ic++){
    CRDEcalEDM::CRDCaloHitLongiCluster m_newClus; m_newClus.Clear();
    
    for(int is=0; is<m_oldClusCol[ic].getBarShowers().size(); is++){
      CRDEcalEDM::CRDCaloBarShower m_shower = m_oldClusCol[ic].getBarShowers()[is];
      
      CRDEcalEDM::CRDCaloBarShower m_selshower; m_selshower.Clear();
      for(int js=0; js<m_showers.size(); js++){
        if( m_showers[js].getModule()==m_shower.getModule() && 
            m_showers[js].getStave()==m_shower.getStave() && 
            m_showers[js].getPart()==m_shower.getPart() && 
            m_showers[js].getDlayer()==m_shower.getDlayer() && 
            (m_showers[js].getPosV3()-m_shower.getPosV3()).Mag()<10 ) {m_selshower = m_showers[js]; break; }
      }
      m_newClus.AddBarShower( m_selshower );
    }
    m_newClus.FitAxis();
    m_outClusCol.push_back( m_newClus );
  }

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::Clustering( std::vector<CRDEcalEDM::CRDCaloBar>& barCol, 
                                           std::vector<CRDEcalEDM::CRDCaloBarCluster>& outClus,  
                                           std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClusCol )
{
  if(barCol.size()==0) return StatusCode::SUCCESS;

  //Neighbor clustering
  std::sort(barCol.begin(), barCol.end());
  std::vector<CRDEcalEDM::CRDCaloBarCluster> m_clusCol; m_clusCol.clear();
  for(int i=0;i<barCol.size();i++){
    CRDEcalEDM::CRDCaloBar iBar = barCol[i];
    if((m_clusCol.size()!=0) && (m_clusCol[m_clusCol.size()-1].isNeighbor(iBar)) ){
      m_clusCol[m_clusCol.size()-1].addBar( iBar );
      continue;
    }

    CRDEcalEDM::CRDCaloBarCluster clus;
    clus.addBar(iBar);
    m_clusCol.push_back(clus);
  }

  //Save out CaloBars in longiClusCol as seed. 
  std::vector<CRDEcalEDM::CRDCaloBar> m_seedbars; m_seedbars.clear(); 
  for(int il=0; il<m_longiClusCol.size(); il++){
  for(int is=0; is<m_longiClusCol[il].getBarShowers().size(); is++){
    std::vector<CRDEcalEDM::CRDCaloBar> m_bars = m_longiClusCol[il].getBarShowers()[is].getBars();
    m_seedbars.insert(m_seedbars.end(), m_bars.begin(), m_bars.end());
  }}

  //Find seed with LongiCluster
  for(int ic=0; ic<m_clusCol.size(); ic++){
  for(int ib=0; ib<m_clusCol[ic].getBars().size(); ib++){
      std::vector<CRDEcalEDM::CRDCaloBar>::iterator iter = find(m_seedbars.begin(), m_seedbars.end(), m_clusCol[ic].getBars()[ib]);
      if(iter==m_seedbars.end()) continue; 
      m_clusCol[ic].addSeed( *iter );
    
  }}

  //Merge clusters without seed
  for(int ic=0; ic<m_clusCol.size(); ic++){
    if(m_clusCol[ic].getNseeds()==0){
      MergeToClosestCluster( m_clusCol[ic], m_clusCol );
      m_clusCol.erase(m_clusCol.begin()+ic);
      ic--;
    }
  }
  outClus = m_clusCol; 

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::ClusterSplitting( CRDEcalEDM::CRDCaloBarCluster& m_cluster, std::vector<CRDEcalEDM::CRDCaloBarShower>& outshCol ){
  std::vector<CRDEcalEDM::CRDCaloBarShower> m_showers; m_showers.clear();

  //No seed in cluster: return empty vector.
  if(m_cluster.getNseeds()==0) { std::cout<<"WARNING: Still have no-seed cluster!!"<<std::endl; outshCol = m_showers; return StatusCode::SUCCESS; }

  //1 seed or second moment less than threshold: Not split. Turn cluster to shower and return
  else if(m_cluster.getNseeds()<2 || m_cluster.getScndMoment()<settings.th_split){
    CRDEcalEDM::CRDCaloBarShower shower; shower.Clear();
    shower.setBars(m_cluster.getBars());
    if(m_cluster.getNseeds()!=0) shower.setSeed(m_cluster.getSeeds()[0]);
    m_showers.push_back(shower);
    outshCol = m_showers;
    return StatusCode::SUCCESS;
  }
  
  //2 or more seeds, large second moment: Split
  int Nshower = m_cluster.getNseeds();
  int Nbars = m_cluster.getBars().size();
  double Eseed[Nshower] = {0};
  double weight[Nbars][Nshower] = {0};
  dd4hep::Position SeedPos[Nshower];
  for(int is=0;is<Nshower;is++) SeedPos[is] = m_cluster.getSeeds()[is].getPosition();
  CalculateInitialEseed(m_cluster.getSeeds(), SeedPos, Eseed);

  bool isConverge = false;
  int iter=0;
  dd4hep::Position SeedPos_prev[Nshower];
  do{

    for(int ibar=0;ibar<m_cluster.getBars().size();ibar++){
      double Eexp[Nshower];
      double Eexp_tot=0;
      for(int is=0;is<Nshower;is++){ Eexp[is] = Eseed[is]*GetShowerProfile(m_cluster.getBars()[ibar].getPosition(), SeedPos[is] ); Eexp_tot+= Eexp[is];}
      for(int is=0;is<Nshower;is++) weight[ibar][is] = Eexp[is]/Eexp_tot;
    }
    for(int is=0;is<Nshower;is++){
      SeedPos_prev[is]=SeedPos[is];

      std::vector<CRDEcalEDM::CRDCaloBar> Bars = m_cluster.getBars();
      for(int ib=0;ib<Bars.size();ib++)  Bars[ib].setQ( Bars[ib].getQ1()*weight[ib][is], Bars[ib].getQ2()*weight[ib][is]  );

      CRDEcalEDM::CRDCaloBarShower shower;
      shower.setBars(Bars);
      SeedPos[is] = shower.getPos();

      double Emax = -99;
      for(int ib=0;ib<Bars.size();ib++) if(Bars[ib].getEnergy()>Emax) Emax=Bars[ib].getEnergy();
      Eseed[is] = Emax;
    }

    isConverge=true;
    for(int is=0;is<Nshower;is++) if( (SeedPos_prev[is]-SeedPos[is]).Mag2()>2.89 ){ isConverge=false; break;}
    iter++;
  }
  while(iter<20 && !isConverge);
  if(iter>=20){
    std::cout<<"WARNING: Iteration time larger than 20! Might not converge!"<<std::endl;
    std::cout<<"  For Check: NBars: "<<m_cluster.getBars().size()<<"  Nseeds: "<<Nshower<<std::endl;
  }

  for(int is=0;is<Nshower;is++){
    CRDEcalEDM::CRDCaloBarShower shower;

    std::vector<CRDEcalEDM::CRDCaloBar> Bars = m_cluster.getBars();
    int iseed=-1;
    double _Emax = -99;
    for(int ib=0;ib<Bars.size();ib++){
      Bars[ib].setQ( Bars[ib].getQ1()*weight[ib][is], Bars[ib].getQ2()*weight[ib][is]  );
      if( Bars[ib].getEnergy()>_Emax ) { _Emax=Bars[ib].getEnergy(); iseed=ib; }
    }
    if(iseed<0) { std::cout<<"ERROR: Can not find seed(max energy bar) in this shower! Please Check!"<<std::endl; iseed=0;}

    shower.setBars(Bars);
    shower.setSeed(Bars[iseed]);
    m_showers.push_back(shower);
  }

  outshCol = m_showers;
  
  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::MergeToClosestCluster( CRDEcalEDM::CRDCaloBarCluster& iclus, std::vector<CRDEcalEDM::CRDCaloBarCluster>& clusvec ){

  int cLedge = iclus.getLeftEdge();
  int cRedge = iclus.getRightEdge();

  //Find the closest cluster with iclus. 
  int minD = 99;
  int index = -1;
  for(int icl=0; icl<clusvec.size(); icl++){
    if(clusvec[icl].getNseeds()==0) continue;
    int iLedge = clusvec[icl].getLeftEdge();
    int iRedge = clusvec[icl].getRightEdge();

    int dis = (cLedge-iRedge>0 ? cLedge-iRedge : iLedge-cRedge );
    if(dis>5) continue; //Don't merge to a too far cluster. 
    if(dis<minD && clusvec[icl].getE()>2.*iclus.getE() ){ minD = dis; index=icl; } //Don't merge to a too small cluster. 
  }
  if(index<0) return StatusCode::FAILURE;

  //Merge to the selected cluster
  for(int icl=0; icl<iclus.getBars().size(); icl++)  clusvec[index].addBar(iclus.getBars()[icl]);

  return StatusCode::SUCCESS;
}


void EnergySplittingAlg::CalculateInitialEseed( const std::vector<CRDEcalEDM::CRDCaloBar>& Seeds, const dd4hep::Position* pos, double* Eseed){
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
    vecE[i] = Seeds[i].getEnergy();
    for(int j=0;j<Nele;j++) matrixR[i][j] = GetShowerProfile(Seeds[i].getPosition(), pos[j]);
  }

  matrixR.Invert();
  TVector sol = matrixR*vecE;

  for(int i=0;i<Nele;i++) Eseed[i] = sol[i];

}


double EnergySplittingAlg::GetShowerProfile(const dd4hep::Position& p_bar, const dd4hep::Position& p_seed ){
  dd4hep::Position rpos = p_bar-p_seed;
  double dis = sqrt(rpos.Mag2());
  double a1, a2, b1, b2, Rm;
  a1=0.037; a2=0.265; b1=0.101; b2=0.437; Rm=1.868;

  return a1*exp(-b1*dis/Rm)+a2*exp(-b2*dis/Rm);
}
#endif
