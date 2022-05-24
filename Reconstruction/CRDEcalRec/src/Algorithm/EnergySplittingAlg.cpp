#ifndef _ENERGYSPLITTING_ALG_C
#define _ENERGYSPLITTING_ALG_C

#include "Algorithm/EnergySplittingAlg.h"
void EnergySplittingAlg::Settings::SetInitialValue(){

  th_split = -1; 
  Eth_SeedAbs = 0.005; 
  Debug = 0;
}

StatusCode EnergySplittingAlg::Initialize(){

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::RunAlgorithm( EnergySplittingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){
  settings = m_settings;
cout<<"EnergySplittingAlg: TowerCol size = "<<m_datacol.TowerCol.size()<<endl;

  for(int it=0; it<m_datacol.TowerCol.size(); it++){
    std::vector<CRDEcalEDM::CRDCaloBlock> m_blocksInTower = m_datacol.TowerCol[it].getBlocks(); 
    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_clusXCol = m_datacol.TowerCol[it].getLongiClusterXCol();
    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_clusYCol = m_datacol.TowerCol[it].getLongiClusterYCol();
    if(m_clusXCol.size()==0 || m_clusYCol.size()==0) continue;

    int minLayerX = 99;
    int minLayerY = 99;
    int maxLayerX = -99;
    int maxLayerY = -99;
    for(int ic=0; ic<m_clusXCol.size(); ic++){
      if(minLayerX>m_clusXCol[ic].getBeginningDlayer()) minLayerX = m_clusXCol[ic].getBeginningDlayer();
      if(maxLayerX<m_clusXCol[ic].getEndDlayer()) maxLayerX = m_clusXCol[ic].getEndDlayer();
    }
    for(int ic=0; ic<m_clusYCol.size(); ic++){
      if(minLayerY>m_clusYCol[ic].getBeginningDlayer()) minLayerY = m_clusYCol[ic].getBeginningDlayer();
      if(maxLayerY<m_clusYCol[ic].getEndDlayer()) maxLayerY = m_clusYCol[ic].getEndDlayer();
    }

printf("  In Tower #%d: HoughCLX range: (%d, %d), HoughCLY range: (%d, %d) \n", it, minLayerX, maxLayerX, minLayerY, maxLayerY);

    //Make clusters and shwoers in blocks
    for(int ib=0; ib<m_blocksInTower.size(); ib++ ){
      m_blocksInTower[ib].ClearShower();

      std::vector<CRDEcalEDM::CRDCaloBarCluster> m_barClusXCol; m_barClusXCol.clear();
      std::vector<CRDEcalEDM::CRDCaloBarCluster> m_barClusYCol; m_barClusYCol.clear();
      std::vector<CRDEcalEDM::CRDCaloBar> m_barColX = m_blocksInTower[ib].getBarXCol();
      std::vector<CRDEcalEDM::CRDCaloBar> m_barColY = m_blocksInTower[ib].getBarYCol();

      //Neighbor clustering
      int dlayer = m_blocksInTower[ib].getDlayer();


      if( dlayer>maxLayerX || dlayer<minLayerX ) Clustering( m_barColX, m_barClusXCol );
      else Clustering( m_barColX, m_barClusXCol, m_clusXCol );

      if( dlayer>maxLayerY || dlayer<minLayerY ) Clustering( m_barColY, m_barClusYCol );
      else Clustering( m_barColY, m_barClusYCol, m_clusYCol );

      m_blocksInTower[ib].setClusterXCol(m_barClusXCol);
      m_blocksInTower[ib].setClusterYCol(m_barClusYCol);


      //Split cluster to showers
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
     
      m_blocksInTower[ib].setShowerXCol( m_barShowerXCol );
      m_blocksInTower[ib].setShowerYCol( m_barShowerYCol );     

printf("    Clustering in Block #%d: Layer = %d, bar size (%d, %d), cluster size (%d, %d), shower size (%d, %d) \n", ib, dlayer, m_barColX.size(), m_barColY.size(), m_barClusXCol.size(), m_barClusYCol.size(), m_barShowerXCol.size(), m_barShowerYCol.size());

    }

cout<<"  EnergySplittingAlg: Update clusters energy"<<endl;
    std::vector<CRDCaloHitLongiCluster> m_longiClusXCol; m_longiClusXCol.clear();
    std::vector<CRDCaloHitLongiCluster> m_longiClusYCol; m_longiClusYCol.clear();
    LongiClusterLinking(m_blocksInTower, m_clusXCol, m_longiClusXCol);
    LongiClusterLinking(m_blocksInTower, m_clusYCol, m_longiClusYCol);

    m_datacol.TowerCol[it].SetLongiClusters( m_longiClusXCol, m_longiClusYCol );
    m_datacol.TowerCol[it].SetBlocks( m_blocksInTower );
  }
cout<<"End EnergySplittingAlg"<<endl;

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::LongiClusterLinking( std::vector<CRDEcalEDM::CRDCaloBlock>& m_blocks, 
                                                    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_oldClusCol, 
                                                    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_outClusCol )
{
  if(m_blocks.size()==0 || m_oldClusCol.size()==0) return StatusCode::SUCCESS;
  m_outClusCol.clear();

  bool fl_isXclus = (m_oldClusCol[0].getSlayer()==0);

  std::vector<CRDEcalEDM::CRDCaloBarShower> m_showers; m_showers.clear(); 
  std::vector<CRDEcalEDM::CRDCaloBarShower> m_showerswotrk; m_showerswotrk.clear(); 
cout<<"  LongiClusterLinking: input block Dlayer: ";
  for(int ib=0; ib<m_blocks.size(); ib++){
cout<<m_blocks[ib].getDlayer()<<"(";

    std::vector<CRDEcalEDM::CRDCaloBarShower> tmp_showersinblock; tmp_showersinblock.clear(); 
    if(fl_isXclus) tmp_showersinblock = m_blocks[ib].getShowerXCol();
    else tmp_showersinblock = m_blocks[ib].getShowerYCol();

    int dlayer = m_blocks[ib].getDlayer();
    bool fl_hastrk = false; 
    for(int icl=0; icl<m_oldClusCol.size(); icl++) 
      if(dlayer<=m_oldClusCol[icl].getEndDlayer() && dlayer>=m_oldClusCol[icl].getBeginningDlayer()) { fl_hastrk=true; break; }
cout<<fl_hastrk<<", shower size="<<tmp_showersinblock.size()<<")"<<'\t';

    if(fl_hastrk) m_showers.insert(m_showers.end(), tmp_showersinblock.begin(), tmp_showersinblock.end());
    else m_showerswotrk.insert(m_showerswotrk.end(), tmp_showersinblock.begin(), tmp_showersinblock.end());
  }

cout<<endl;
cout<<"  LongiClusterLinking: #normal shower = "<<m_showers.size()<<", #shower wo trk = "<<m_showerswotrk.size()<<endl;

cout<<"  Normal shower position: "<<endl;
for(int is=0; is<m_showers.size(); is++) printf("    (%.2f, %.2f, %.2f) \n", m_showers[is].getPos().x(), m_showers[is].getPos().y(), m_showers[is].getPos().z());

  //Update old Hough clusters
  for(int ic=0; ic<m_oldClusCol.size(); ic++){
    CRDEcalEDM::CRDCaloHitLongiCluster m_newClus; m_newClus.Clear();
cout<<"  Showers in cluster "<<ic<<endl;
    for(int is=0; is<m_oldClusCol[ic].getBarShowers().size(); is++){
      CRDEcalEDM::CRDCaloBarShower m_shower = m_oldClusCol[ic].getBarShowers()[is];
printf("    (%.2f, %.2f, %.2f) \n", m_shower.getPos().x(), m_shower.getPos().y(), m_shower.getPos().z());

      CRDEcalEDM::CRDCaloBarShower m_selshower; m_selshower.Clear();
      bool fl_foundshower = false; 
      for(int js=0; js<m_showers.size(); js++){
        if( m_showers[js].getModule()==m_shower.getModule() && 
            m_showers[js].getStave()==m_shower.getStave() && 
            m_showers[js].getPart()==m_shower.getPart() && 
            m_showers[js].getDlayer()==m_shower.getDlayer() && 
            (m_showers[js].getPos()-m_shower.getPos()).Mag()<10 ) {m_selshower = m_showers[js]; fl_foundshower=true; break; }
      }
      if(fl_foundshower) m_newClus.AddBarShower( m_selshower );
else{ printf(" In Cluster #%d: shower #%d does not find the new shower. Layer = %d. \n ", ic, is, m_shower.getDlayer()); }
    }

    m_newClus.SetHoughPars(m_oldClusCol[ic].getHoughAlpha(), m_oldClusCol[ic].getHoughRho() );
    m_newClus.SetIntercept(m_oldClusCol[ic].getHoughIntercept());
    m_newClus.FitAxis();
    m_outClusCol.push_back( m_newClus );
  }

cout<<"    Before merging shower: "<<endl;
cout<<"    HoughCluster Nhit: ";
for(int icl=0; icl<m_outClusCol.size(); icl++) cout<<m_outClusCol[icl].getBarShowers().size()<<'\t';
cout<<endl;

  //Merge showers into closest Hough cluster
  for(int is=0; is<m_showerswotrk.size(); is++) 
    MergeToClosestCluster( m_showerswotrk[is], m_outClusCol );

cout<<"    After merging shower: "<<endl;
cout<<"    HoughCluster Nhit: ";
for(int icl=0; icl<m_outClusCol.size(); icl++) cout<<m_outClusCol[icl].getBarShowers().size()<<'\t';
cout<<endl;

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
//cout<<"Clustering: After Neighbor clustering. Cluster size: "<<m_clusCol.size()<<endl;
//cout<<"  Bar size in cluster: ";
//for(int j=0; j<m_clusCol.size(); j++) cout<<m_clusCol[j].getBars().size()<<'\t';
//cout<<endl;

  //Save out CaloBars in longiClusCol as seed. 
  std::vector<CRDEcalEDM::CRDCaloBar> m_seedbars; m_seedbars.clear(); 
  for(int il=0; il<m_longiClusCol.size(); il++){
  for(int is=0; is<m_longiClusCol[il].getBarShowers().size(); is++){
    std::vector<CRDEcalEDM::CRDCaloBar> m_bars = m_longiClusCol[il].getBarShowers()[is].getBars();
    m_seedbars.insert(m_seedbars.end(), m_bars.begin(), m_bars.end());
  }}
//cout<<"Clustering: seed size from LongiCluster: "<<m_seedbars.size()<<endl;

  //Find seed with LongiCluster
  for(int ic=0; ic<m_clusCol.size(); ic++){
  for(int ib=0; ib<m_clusCol[ic].getBars().size(); ib++){
      std::vector<CRDEcalEDM::CRDCaloBar>::iterator iter = find(m_seedbars.begin(), m_seedbars.end(), m_clusCol[ic].getBars()[ib]);
      if(iter==m_seedbars.end()) continue; 
      m_clusCol[ic].addSeed( *iter );
    
  }}
//cout<<"Clustering: After seed finding. Cluster size: "<<m_clusCol.size()<<endl;
//cout<<"  Bar size in cluster: ";
//for(int j=0; j<m_clusCol.size(); j++) cout<<m_clusCol[j].getBars().size()<<'\t';
//cout<<endl;
//cout<<"  Seed size in cluster: ";
//for(int j=0; j<m_clusCol.size(); j++) cout<<m_clusCol[j].getSeeds().size()<<'\t';
//cout<<endl;

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


StatusCode EnergySplittingAlg::Clustering( std::vector<CRDEcalEDM::CRDCaloBar>& barCol,
                                           std::vector<CRDEcalEDM::CRDCaloBarCluster>& outClus){
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

  //Find seed in cluster
  for(int i=0;i<m_clusCol.size();i++){
    CRDEcalEDM::CRDCaloBarCluster iclus = m_clusCol[i];

    std::vector<CRDEcalEDM::CRDCaloBar> m_seedVec; m_seedVec.clear();
    findSeeds(iclus, m_seedVec);

    m_clusCol[i].setSeeds(m_seedVec);
    m_clusCol[i].getScndMoment();
  }

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
//cout<<"ClusterSplitting: input cluster seed size = "<<m_cluster.getNseeds()<<", second moment = "<<m_cluster.getScndMoment()<<endl;

  //No seed in cluster: return empty vector.
  if(m_cluster.getNseeds()==0) { std::cout<<"WARNING: Still have no-seed cluster!!"<<std::endl; outshCol = m_showers; return StatusCode::SUCCESS; }

  //1 seed or second moment less than threshold: Not split. Turn cluster to shower and return
  else if(m_cluster.getNseeds()<2 || m_cluster.getScndMoment()<settings.th_split){
    CRDEcalEDM::CRDCaloBarShower shower; shower.Clear();
    std::vector<CRDEcalEDM::CRDCaloBar> tmp_bars = m_cluster.getBars(); 
    shower.setBars(tmp_bars);
    if(m_cluster.getNseeds()!=0) shower.setSeed(m_cluster.getSeeds()[0]);
    shower.setIDInfo();
    m_showers.push_back(shower);
    outshCol = m_showers;
    return StatusCode::SUCCESS;
  }
  
  //2 or more seeds, large second moment: Split
  int Nshower = m_cluster.getNseeds();
  int Nbars = m_cluster.getBars().size();
  double Eseed[Nshower] = {0};
  double weight[Nbars][Nshower] = {0};
  TVector3 SeedPos[Nshower];
  for(int is=0;is<Nshower;is++) SeedPos[is] = m_cluster.getSeeds()[is].getPosition();
  CalculateInitialEseed(m_cluster.getSeeds(), SeedPos, Eseed);

  bool isConverge = false;
  int iter=0;
  TVector3 SeedPos_prev[Nshower];
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
    shower.setIDInfo(); 
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


StatusCode EnergySplittingAlg::MergeToClosestCluster( CRDEcalEDM::CRDCaloBarShower& m_shower, 
                                                      std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_clusters ){
  if(m_clusters.size()==0) return StatusCode::SUCCESS;


  int minLayer = 99;
  int maxLayer = -99;
  for(int ic=0; ic<m_clusters.size(); ic++){
    if(minLayer>m_clusters[ic].getBeginningDlayer()) minLayer = m_clusters[ic].getBeginningDlayer();
    if(maxLayer<m_clusters[ic].getEndDlayer()) maxLayer = m_clusters[ic].getEndDlayer();
  }
  if(minLayer==99 || maxLayer<0) return StatusCode::SUCCESS;

  int dlayer = m_shower.getDlayer();
  TVector3 sh_pos = m_shower.getPos();

cout<<"  Cluster range: ("<<minLayer<<", "<<maxLayer<<") "<<endl;
cout<<"  Merging shower into cluster: Input shower ";
printf(" (%.3f, %.3f, %.3f, %.3f), layer #%d \n", sh_pos.X(), sh_pos.Y(), sh_pos.Z(), m_shower.getE(), dlayer);

  double minR = 999;
  int index_cluster = -1;
  double m_distance = 0;
  for(int ic=0; ic<m_clusters.size(); ic++ ){
    if(dlayer<minLayer)
      m_distance = (sh_pos-m_clusters[ic].getBarShowersInLayer(m_clusters[ic].getBeginningDlayer())[0].getPos()).Mag();     
    else if(dlayer>maxLayer)
      m_distance = (sh_pos-m_clusters[ic].getBarShowersInLayer(m_clusters[ic].getEndDlayer())[0].getPos()).Mag();    
    else
      m_distance = (sh_pos-m_clusters[ic].getPos()).Mag();

cout<<"    Cluster #"<<ic<<": distance with shower = "<<m_distance<<endl;
    
    if( m_distance<minR )  { minR=m_distance; index_cluster=ic; }
  }

printf("  minR = %.3f, in #cl %d", minR, index_cluster);

  if(index_cluster>=0) m_clusters[index_cluster].AddBarShower( m_shower );

  return StatusCode::SUCCESS;
}

StatusCode EnergySplittingAlg::findSeeds( CRDEcalEDM::CRDCaloBarCluster& m_cluster, std::vector<CRDEcalEDM::CRDCaloBar>& seedCol){

  std::vector<CRDEcalEDM::CRDCaloBar> m_seeds; m_seeds.clear();
  for(int i=0;i<m_cluster.getBars().size();i++){
    CRDEcalEDM::CRDCaloBar ibar = m_cluster.getBars()[i];
    std::vector<CRDEcalEDM::CRDCaloBar> m_neighbor = getNeighbors(m_cluster, ibar);
    if(m_neighbor.size()==0 && ibar.getEnergy()>settings.Eth_SeedAbs)
      m_seeds.push_back(ibar); 

    else {
      bool isLocalMax=true; 
      for(int j=0;j<m_neighbor.size();j++){
        if(m_neighbor[j].getEnergy()>ibar.getEnergy()) isLocalMax=false;
      }
      if(ibar.getEnergy()>settings.Eth_SeedAbs && isLocalMax) m_seeds.push_back(ibar);
    }
  }

  seedCol = m_seeds;
  return StatusCode::SUCCESS;
}

std::vector<CRDEcalEDM::CRDCaloBar>  EnergySplittingAlg::getNeighbors(CRDEcalEDM::CRDCaloBarCluster& m_cluster, CRDEcalEDM::CRDCaloBar& seed){

  std::vector<CRDEcalEDM::CRDCaloBar> m_neighbor;
  for(int i=0;i<m_cluster.getBars().size();i++){

    if( seed.isNeighbor(m_cluster.getBars()[i]) ) m_neighbor.push_back(m_cluster.getBars()[i]);
  }
  if(m_neighbor.size()>2) std::cout<<"WARNING: more than 2 hits in neighborCol!!"<<std::endl;

  return m_neighbor;
}


void EnergySplittingAlg::CalculateInitialEseed( const std::vector<CRDEcalEDM::CRDCaloBar>& Seeds, const TVector3* pos, double* Eseed){
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


double EnergySplittingAlg::GetShowerProfile(const TVector3& p_bar, const TVector3& p_seed ){
  TVector3 rpos = p_bar-p_seed;
  double dis = sqrt(rpos.Mag2());
  double a1, a2, b1, b2, Rm;
  a1=0.037; a2=0.265; b1=0.101; b2=0.437; Rm=1.868;

  return a1*exp(-b1*dis/Rm)+a2*exp(-b2*dis/Rm);
}
#endif
