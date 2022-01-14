#ifndef _ENERGYSPLITTING_ALG_C
#define _ENERGYSPLITTING_ALG_C

#include "Algorithm/EnergySplittingAlg.h"
void EnergySplittingAlg::Settings::SetInitialValue(){
  fl_findMaxOnly = true; 
  Eth_localMax = 0.003;
  Eth_MaxWithNeigh = 0.;

  Sth_split = -1;
  Eth_SeedAbs = 0.007; //0.1MeV
  Eth_ShowerAbs = 0.0;
  Eth_ClusAbs = 0.0;
  Eth_SeedWithNeigh = 0.6;
  Eth_SeedWithTot = 0.;
  Eth_ShowerWithTot = 0.;
  Eth_ClusterWithTot = 0.;
  Etot=0;
  Debug = 0;
}

StatusCode EnergySplittingAlg::Initialize(){

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::RunAlgorithm( EnergySplittingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){
  settings = m_settings;

  std::vector<CRDEcalEDM::CRDCaloLayer> m_LayerCol; m_LayerCol.clear();
  std::vector<CRDEcalEDM::CRDCaloTower> m_TowerCol; m_TowerCol.clear();
  int Nblocks = m_datacol.BlockVec.size();
  for(int ib=0;ib<Nblocks;ib++){
    CRDEcalEDM::CRDCaloLayer m_layer; m_layer.Clear();
    if( settings.fl_findMaxOnly ) GetLocalMax( m_datacol.BlockVec[ib] );
    else ClusteringinLayer( m_datacol.BlockVec[ib], m_layer );
    m_LayerCol.push_back(m_layer);

    CRDEcalEDM::CRDCaloTower m_tower; m_tower.Clear();
    m_tower.SetTowerID( m_datacol.BlockVec[ib].getModule(), 
                        m_datacol.BlockVec[ib].getStave(),
                        m_datacol.BlockVec[ib].getPart() );
    std::vector<CRDEcalEDM::CRDCaloTower>::iterator iter = find(m_TowerCol.begin(), m_TowerCol.end(), m_tower);
    if(iter==m_TowerCol.end()){ m_tower.AddBlock(m_datacol.BlockVec[ib]); m_TowerCol.push_back(m_tower); }
    else iter->AddBlock(m_datacol.BlockVec[ib]);
    
  }
  m_datacol.LayerCol = m_LayerCol;
  m_datacol.TowerCol = m_TowerCol;
  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}


StatusCode EnergySplittingAlg::GetLocalMax( CRDEcalEDM::CRDCaloBlock& m_block ){
  if(m_block.getBarXCol().size()==0 && m_block.getBarYCol().size()==0) return StatusCode::SUCCESS;

  std::vector<CRDEcalEDM::CRDCaloBar> m_barXCol = m_block.getBarXCol();
  std::vector<CRDEcalEDM::CRDCaloBar> m_barYCol = m_block.getBarYCol();

  std::vector<CRDEcalEDM::CRDCaloBar> localMaxXCol; localMaxXCol.clear();
  std::vector<CRDEcalEDM::CRDCaloBar> localMaxYCol; localMaxYCol.clear();
  GetLocalMaxBar( m_barXCol, localMaxXCol );
  GetLocalMaxBar( m_barYCol, localMaxYCol );


  std::vector<CRDEcalEDM::CRDCaloBarShower> m_showerColX; m_showerColX.clear();  
  std::vector<CRDEcalEDM::CRDCaloBarShower> m_showerColY; m_showerColY.clear();  
  for(int j=0; j<localMaxXCol.size(); j++){
    CRDEcalEDM::CRDCaloBarShower m_shower; m_shower.Clear();
    m_shower.addBar( localMaxXCol[j] );
    m_shower.setSeed( localMaxXCol[j] );
    m_shower.setIDInfo( localMaxXCol[j].getModule(),    
                        localMaxXCol[j].getStave(),
                        localMaxXCol[j].getDlayer(),
                        localMaxXCol[j].getPart(),
                        localMaxXCol[j].getSlayer() );
    m_showerColX.push_back(m_shower);
  }
  for(int j=0; j<localMaxYCol.size(); j++){
    CRDEcalEDM::CRDCaloBarShower m_shower; m_shower.Clear();
    m_shower.addBar( localMaxYCol[j] );
    m_shower.setSeed( localMaxYCol[j] );
    m_shower.setIDInfo( localMaxYCol[j].getModule(),
                        localMaxYCol[j].getStave(),
                        localMaxYCol[j].getDlayer(),
                        localMaxYCol[j].getPart(),
                        localMaxYCol[j].getSlayer() );
    m_showerColY.push_back(m_shower);
  }

  m_block.setShowerXCol(m_showerColX);
  m_block.setShowerYCol(m_showerColY);
  return StatusCode::SUCCESS;
}

StatusCode EnergySplittingAlg::GetLocalMaxBar( std::vector<CRDEcalEDM::CRDCaloBar>& barCol, std::vector<CRDEcalEDM::CRDCaloBar>& localMaxCol ){
  std::sort( barCol.begin(), barCol.end() );

  for(int ib=0; ib<barCol.size(); ib++){
    std::vector<CRDEcalEDM::CRDCaloBar> m_neighbors = getNeighbors( barCol[ib], barCol );
    if( m_neighbors.size()==0 && barCol[ib].getEnergy()>settings.Eth_localMax ) { localMaxCol.push_back( barCol[ib] ); continue; }

    bool isLocalMax=true;
    bool isIso;
    double Eneigh=0;
    if(barCol[ib].getEnergy()<settings.Eth_localMax) isLocalMax = false;
    for(int j=0;j<m_neighbors.size();j++){
      if(m_neighbors[j].getEnergy()>barCol[ib].getEnergy()) isLocalMax=false;
      Eneigh += m_neighbors[j].getEnergy();
    }
    isIso = (barCol[ib].getEnergy()/(barCol[ib].getEnergy()+Eneigh))>settings.Eth_MaxWithNeigh;
    if(!isIso) continue;

    if(isLocalMax) localMaxCol.push_back( barCol[ib] );

  }

  return StatusCode::SUCCESS;
}




StatusCode EnergySplittingAlg::ClusteringinLayer( CRDEcalEDM::CRDCaloBlock& m_blocks, CRDEcalEDM::CRDCaloLayer& m_layer ){
  m_layer.Clear(); 

  m_layer.barXCol = m_blocks.getBarXCol();
  m_layer.barYCol = m_blocks.getBarYCol();

  bool f_succX = true;  bool f_succY = true;
  if( !Clustering( m_layer.barXCol, m_layer.barClusXCol, m_blocks.getAllShadowClusCol()) ) {
       Clustering( m_layer.barXCol, m_layer.barClusXCol ); f_succX = false;
  }
  if( !Clustering( m_layer.barYCol, m_layer.barClusYCol, m_blocks.getAllShadowClusCol()) ) {
       Clustering( m_layer.barYCol, m_layer.barClusYCol ); f_succY = false;
  }

  if(settings.Debug>0){
    printf("EnergySplittingAlg::ClusterinLayer  BlockID(%d, %d, %d, %d), Bar size(X/Y): (%d / %d) \n",
           m_blocks.getModule(), m_blocks.getStave(), m_blocks.getPart(),  m_blocks.getDlayer(), m_layer.barXCol.size(), m_layer.barYCol.size() );
    std::cout<<"  Cluster with candidates in X: "<<f_succX<<"  Cluster size: "<<m_layer.barClusXCol.size()<<std::endl;
    std::cout<<"  Cluster with candidates in Y: "<<f_succY<<"  Cluster size: "<<m_layer.barClusYCol.size()<<std::endl;
  }

  for(int i=0;i<m_layer.barClusXCol.size();i++){
    std::vector<CRDEcalEDM::CRDCaloBarShower> showers; showers.clear();
    //---WIP: Split the cluster with expected showers (position & energy).
//    if(f_succX)
//      ClusterSplitting( m_layer.barClusXCol[i], showers,  m_blocks.getNeuShadowClusCol(), m_blocks.getTrkShadowClusCol() );
//    else
      ClusterSplitting(m_layer.barClusXCol[i], showers);

    if(showers.size()==0) continue;
    m_layer.barShowerXCol.insert(m_layer.barShowerXCol.end(), showers.begin(), showers.end());
  }

  for(int i=0;i<m_layer.barClusYCol.size();i++){
    std::vector<CRDEcalEDM::CRDCaloBarShower> showers; showers.clear();
//    if(f_succY)
//      ClusterSplitting( m_layer.barClusYCol[i], showers, m_blocks.getNeuShadowClusCol(), m_blocks.getTrkShadowClusCol());
//    else
      ClusterSplitting(m_layer.barClusYCol[i], showers);

    if(showers.size()==0) continue;
    m_layer.barShowerYCol.insert(m_layer.barShowerYCol.end(), showers.begin(), showers.end());
  }

  return StatusCode::SUCCESS;
}


bool EnergySplittingAlg::Clustering( std::vector<CRDEcalEDM::CRDCaloBar>& barCol,
                                     std::vector<CRDEcalEDM::CRDCaloBarCluster>& outClus)
{
  if(barCol.size()==0) return false;
  const std::vector<CRDEcalEDM::CRDShadowCluster> emptyCol(0);
  return Clustering( barCol, outClus, emptyCol );
}


bool EnergySplittingAlg::Clustering( std::vector<CRDEcalEDM::CRDCaloBar>& barCol,
                                     std::vector<CRDEcalEDM::CRDCaloBarCluster>& outClus,
                                     const std::vector<CRDEcalEDM::CRDShadowCluster>& m_CandiCol )
{

  if(settings.Debug>1) std::cout<<"Bar Clustering: barCol size = "<<barCol.size()<<std::endl;
  if(barCol.size()==0) return false;

  if(settings.Debug>2){
    for(int i=0; i<barCol.size(); i++)
      printf("  Bar Pos: (%.2f, %.2f, %.2f, %.5f) \t", barCol[i].getPosition().x(), barCol[i].getPosition().y(), barCol[i].getPosition().z(), barCol[i].getEnergy());
    cout<<endl;
  }

  //Neighbor clustering
  for(int i=0;i<barCol.size();i++) settings.Etot+=barCol[i].getEnergy();
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

  //Find seeds in cluster.
  std::vector<CRDEcalEDM::CRDCaloBarCluster> out_clusCol; out_clusCol.clear();
  for(int i=0;i<m_clusCol.size();i++){
    CRDEcalEDM::CRDCaloBarCluster iclus = m_clusCol[i];
    if( (iclus.getE()/settings.Etot < settings.Eth_ClusterWithTot) || iclus.getE()<settings.Eth_ClusAbs ) continue;
    std::vector<CRDEcalEDM::CRDCaloBar> m_seedVec; m_seedVec.clear();

    m_seedVec = findSeeds(iclus, m_CandiCol);
    iclus.setSeeds(m_seedVec);
    iclus.getScndMoment();
    out_clusCol.push_back(iclus);
  }

  if(settings.Debug>1){
    std::cout<<"  N neiCluster: "<<m_clusCol.size()<<std::endl;
    std::cout<<"  N outCluster: "<<out_clusCol.size()<<std::endl;
    std::cout<<"  N bars in outClus: "; for(int i=0; i<out_clusCol.size(); i++) std::cout<<out_clusCol[i].getBars().size()<<"  "; std::cout<<endl;
    std::cout<<"  N seed in outClus: "; for(int i=0; i<out_clusCol.size(); i++) std::cout<<out_clusCol[i].getNseeds()<<"  "; std::cout<<endl;
  }

  //Merge Clusters without seed
  for(int i=0; i<out_clusCol.size(); i++){
    if(out_clusCol[i].getNseeds()==0){
      MergeToClosestCluster( out_clusCol[i], out_clusCol );
      out_clusCol.erase(out_clusCol.begin()+i);
      i--;
    }
  }

  outClus = out_clusCol;
  return true;
}



bool EnergySplittingAlg::ClusterSplitting( CRDEcalEDM::CRDCaloBarCluster& m_cluster, std::vector<CRDEcalEDM::CRDCaloBarShower>& outshCol ){

  std::vector<CRDEcalEDM::CRDCaloBarShower> m_showers; m_showers.clear();

  if(settings.Debug>1)
    printf("Cluster for Splitting: (%.2f, %.2f, %.2f), scndM=%.2f \n", m_cluster.getPos().x(), m_cluster.getPos().y(), m_cluster.getPos().z(), m_cluster.getScndMoment());

  //No seed in cluster: return empty vector.
  if(m_cluster.getNseeds()==0) { std::cout<<"WARNING: Still have no-seed cluster!!"<<std::endl; outshCol = m_showers; return false; }

  //1 seed or second moment less than threshold: Not split. Turn cluster to shower and return
  else if(m_cluster.getNseeds()<2 || m_cluster.getScndMoment()<settings.Sth_split){
    CRDEcalEDM::CRDCaloBarShower shower; shower.Clear();
    shower.setBars(m_cluster.getBars());
    if(m_cluster.getNseeds()!=0) shower.setSeed(m_cluster.getSeeds()[0]);
    if( (shower.getE()/settings.Etot > settings.Eth_ShowerWithTot) && shower.getE()>settings.Eth_ShowerAbs  ) m_showers.push_back(shower);
    outshCol = m_showers;
    return true;
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
    if(settings.Debug>1){
      std::cout<<"  Seed position in last iter: "; for(int is=0; is<Nshower; is++) printf("(%.2f, %.2f, %.2f)  ",SeedPos[is].x(), SeedPos[is].y(),SeedPos[is].z()); std::cout<<std::endl;
      std::cout<<"  Seed position in 2nd last iter: "; for(int is=0; is<Nshower; is++) printf("(%.2f, %.2f, %.2f)  ",SeedPos[is].x(), SeedPos[is].y(),SeedPos[is].z()); std::cout<<std::endl;
    }
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
  return true;

}

bool EnergySplittingAlg::ClusterSplitting( CRDEcalEDM::CRDCaloBarCluster& m_cluster,
                                           std::vector<CRDEcalEDM::CRDCaloBarShower>& outshCol,
                                           const std::vector<CRDEcalEDM::CRDShadowCluster>& m_Neucandidates,
                                           const std::vector<CRDEcalEDM::CRDShadowCluster>& m_Trkcandidates  )
{
  std::vector<CRDEcalEDM::CRDCaloBarShower> m_showers; m_showers.clear();
  bool f_succ = ClusterSplitting( m_cluster, m_showers );
  if(!f_succ){ outshCol=m_showers;  return false;}

  int slayer = m_showers[0].getBars()[0].getSlayer();
  int module = m_showers[0].getBars()[0].getModule();
  int f_isXclus = (slayer==0 ? 1 : 0);
  double rotAngle = -module*PI/4.;

  if(settings.Debug>1) std::cout<<"ClusterSplitting Adding ShadowClus: isXclus "<<f_isXclus<<"  shower size: "<<m_showers.size()<<std::endl;

  for(int is=0; is<m_showers.size(); is++){
    TVector3 p_shower(0, 0, 0);
    p_shower.SetXYZ( m_showers[is].getPos().x(), m_showers[is].getPos().y(), m_showers[is].getPos().z() );

    if(settings.Debug>2) printf("  Shower #%d: (%.2f, %.2f, %.2f) \n", is, p_shower.x(), p_shower.y(), p_shower.z() );

    //Neutral candidate
    for(int ic=0; ic<m_Neucandidates.size(); ic++){
      TVector3 p_candi = m_Neucandidates[ic].ExpPos;

      if(settings.Debug>2) printf("    Neutral candidate: ((%.2f, %.2f, %.2f) \n", p_candi.x(), p_candi.y(), p_candi.z());

      if(f_isXclus==1 && fabs(p_candi.z()-p_shower.z())<10 ){  //10mm distance
        m_showers[is].addNeuShadowClus( m_Neucandidates[ic] );
        if(settings.Debug>2) cout<<"    !!Added!"<<endl;
      }
      if(f_isXclus==0){
        p_shower.RotateZ(rotAngle);
        p_candi.RotateZ(rotAngle);
        if( fabs(p_shower.x()-p_candi.x())<10. ){
          m_showers[is].addNeuShadowClus( m_Neucandidates[ic] );
          if(settings.Debug>2) cout<<"    !!Added!"<<endl;
        }
        p_shower.RotateZ(-rotAngle);
        p_candi.RotateZ(-rotAngle);
      }
    }

    //Charged candidate
    p_shower.SetXYZ( m_showers[is].getPos().x(), m_showers[is].getPos().y(), m_showers[is].getPos().z() );
    for(int ic=0; ic<m_Trkcandidates.size(); ic++){
      TVector3 p_candi = m_Trkcandidates[ic].ExpPos;

      if(settings.Debug>2) printf("    Track candidate: ((%.2f, %.2f, %.2f) \n", p_candi.x(), p_candi.y(), p_candi.z());

      if(f_isXclus==1 && fabs(p_candi.z()-p_shower.z())<10 ){
        m_showers[is].addTrkShadowClus( m_Trkcandidates[ic] );
        if(settings.Debug>2) cout<<"    !!Added!"<<endl;
      }
      if(f_isXclus==0){
        p_shower.RotateZ(rotAngle);
        p_candi.RotateZ(rotAngle);
        if( fabs(p_shower.x()-p_candi.x())<10. ){
          m_showers[is].addTrkShadowClus( m_Trkcandidates[ic] );
          if(settings.Debug>2) cout<<"    !!Added!"<<endl;
        }
        p_shower.RotateZ(-rotAngle);
        p_candi.RotateZ(-rotAngle);
      }
    }

  } //End loop showers

  outshCol=m_showers;
  return true;
}

bool EnergySplittingAlg::GetMatchedSeeds( std::vector<CRDEcalEDM::CRDCaloBar>& seeds, 
                                          const std::vector<CRDEcalEDM::CRDShadowCluster>& m_CandiCol, 
                                          std::vector<CRDEcalEDM::CRDCaloBar>& m_matchedCol,
                                          std::vector<CRDEcalEDM::CRDCaloBar>& m_unmatchedCol ) 
{
  m_matchedCol.clear(); m_unmatchedCol.clear(); 

  for(int is=0;is<seeds.size(); is++){
  for(int j=0; j<m_CandiCol.size(); j++){
    if( seeds[is].getSlayer()==0){

    if(settings.Debug>2)
      printf("  Slayer0: Seed(%.2f, %.2f, %.2f).  Candidate(%.2f, %.2f, %.2f) \n",
                seeds[is].getPosition().x(), seeds[is].getPosition().y(), seeds[is].getPosition().z(),
                m_CandiCol[j].ExpPos.x(), m_CandiCol[j].ExpPos.y(), m_CandiCol[j].ExpPos.z());


      if( fabs( seeds[is].getPosition().z()-m_CandiCol[j].ExpPos.z() ) < 12 ) // Position difference within 12mm
        { m_matchedCol.push_back(seeds[is]); break; }
      else 
        { m_unmatchedCol.push_back(seeds[is]); break; }
    }
    else if (seeds[is].getSlayer()==1){

      if(settings.Debug>2)
      printf("  Slayer1: Seed(%.2f, %.2f, %.2f).  Candidate(%.2f, %.2f, %.2f) \n",
              seeds[is].getPosition().x(), seeds[is].getPosition().y(), seeds[is].getPosition().z(),
              m_CandiCol[j].ExpPos.x(), m_CandiCol[j].ExpPos.y(), m_CandiCol[j].ExpPos.z());


      double rotAngle = -seeds[is].getModule()*PI/4.;
      TVector3 p_exps = m_CandiCol[j].ExpPos;
      TVector3 p_seed(0,0,0);
      p_seed.SetXYZ( seeds[is].getPosition().x(), seeds[is].getPosition().y(), seeds[is].getPosition().z() );

      p_exps.RotateZ(rotAngle);
      p_seed.RotateZ(rotAngle);

      if( fabs(p_exps.x()-p_seed.x()) < 12 )
        { m_matchedCol.push_back(seeds[is]); break; }
      else
        { m_unmatchedCol.push_back(seeds[is]); break; }
    }
  }}
  return true;
}


std::vector<CRDEcalEDM::CRDCaloBar>  EnergySplittingAlg::findSeeds( CRDEcalEDM::CRDCaloBarCluster& m_cluster){
  const std::vector<CRDEcalEDM::CRDShadowCluster> emptyCol(0);
  return findSeeds(m_cluster, emptyCol);
}


std::vector<CRDEcalEDM::CRDCaloBar>  EnergySplittingAlg::findSeeds( CRDEcalEDM::CRDCaloBarCluster& m_cluster, const std::vector<CRDEcalEDM::CRDShadowCluster>& m_CandiCol){

  std::vector<CRDEcalEDM::CRDCaloBar> m_seeds; m_seeds.clear();
  std::vector<CRDEcalEDM::CRDCaloBar> m_localMax; m_localMax.clear();
  for(int i=0;i<m_cluster.getBars().size();i++){
    CRDEcalEDM::CRDCaloBar ibar = m_cluster.getBars()[i];
    std::vector<CRDEcalEDM::CRDCaloBar> m_neighbor = getNeighbors(m_cluster, ibar);
    if(m_neighbor.size()==0){
      //std::cout<<"WARNING: An isolated bar without neighbors! "<<std::endl;
      if(ibar.getEnergy()>settings.Eth_localMax) { m_localMax.push_back(ibar); }
      continue;
    }

    bool isLocalMax=true;
    double Eneigh=0;
    for(int j=0;j<m_neighbor.size();j++){
      if(m_neighbor[j].getEnergy()>ibar.getEnergy()) isLocalMax=false;
      Eneigh += m_neighbor[j].getEnergy();
    }
    if(isLocalMax && ibar.getEnergy()>settings.Eth_localMax ) m_localMax.push_back(ibar);
  }

  if(m_CandiCol.size()==0){
    for(int is=0; is<m_localMax.size(); is++)
      if(m_localMax[is].getEnergy()>settings.Eth_SeedAbs) m_seeds.push_back(m_localMax[is]);
    
  }
  else{
    std::vector<CRDEcalEDM::CRDCaloBar> m_matched; m_matched.clear(); 
    std::vector<CRDEcalEDM::CRDCaloBar> m_unmatched; m_unmatched.clear(); 
    GetMatchedSeeds(m_localMax, m_CandiCol, m_matched, m_unmatched);

    m_seeds.insert(m_seeds.end(), m_matched.begin(), m_matched.end());

    for(int is=0; is<m_unmatched.size(); is++){
      std::vector<CRDEcalEDM::CRDCaloBar> m_neighbor = getNeighbors(m_cluster, m_unmatched[is]);
      double Eneigh=0;
      for(int j=0;j<m_neighbor.size();j++){
        Eneigh += m_neighbor[j].getEnergy();
      }
      bool isIso = ( m_unmatched.size()>0 && (m_unmatched[is].getEnergy()/(m_unmatched[is].getEnergy()+Eneigh))>settings.Eth_SeedWithNeigh );
      if(m_unmatched[is].getEnergy()>settings.Eth_SeedAbs && isIso) m_seeds.push_back(m_unmatched[is]);
    }
  }
  return m_seeds; 

}


std::vector<CRDEcalEDM::CRDCaloBar> EnergySplittingAlg::getNeighbors( CRDEcalEDM::CRDCaloBar& seed, std::vector<CRDEcalEDM::CRDCaloBar>& barCol){
  std::vector<CRDEcalEDM::CRDCaloBar> m_neighbor;
  for(int i=0;i<barCol.size();i++){

    if( seed.isNeighbor(barCol[i]) ) m_neighbor.push_back(barCol[i]);
  }
  if(m_neighbor.size()>2) std::cout<<"WARNING: more than 2 hits in neighborCol!!"<<std::endl;

  return m_neighbor;
}

std::vector<CRDEcalEDM::CRDCaloBar>  EnergySplittingAlg::getNeighbors(CRDEcalEDM::CRDCaloBarCluster& m_cluster, CRDEcalEDM::CRDCaloBar& seed){

  std::vector<CRDEcalEDM::CRDCaloBar> m_neighbor;
  for(int i=0;i<m_cluster.getBars().size();i++){

    if( seed.isNeighbor(m_cluster.getBars()[i]) ) m_neighbor.push_back(m_cluster.getBars()[i]);
  }
  if(m_neighbor.size()>2) std::cout<<"WARNING: more than 2 hits in neighborCol!!"<<std::endl;

  return m_neighbor;
}


bool EnergySplittingAlg::MergeToClosestCluster( CRDEcalEDM::CRDCaloBarCluster& iclus, std::vector<CRDEcalEDM::CRDCaloBarCluster>& clusvec ){

  int cLedge = iclus.getLeftEdge();
  int cRedge = iclus.getRightEdge();

  int minD = 99;
  int index = -1;
  for(int icl=0; icl<clusvec.size(); icl++){
    if(clusvec[icl].getNseeds()==0) continue;
    int iLedge = clusvec[icl].getLeftEdge();
    int iRedge = clusvec[icl].getRightEdge();

    int dis = (cLedge-iRedge>0 ? cLedge-iRedge : iLedge-cRedge );
    if(dis>3) continue;
    if(dis<minD){ minD = dis; index=icl; }
  }
  if(index<0) return false;
  if(clusvec[index].getE()/iclus.getE()<2.) return false;

  for(int icl=0; icl<iclus.getBars().size(); icl++)  clusvec[index].addBar(iclus.getBars()[icl]);
  clusvec[index].getScndMoment();
  return true;

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
