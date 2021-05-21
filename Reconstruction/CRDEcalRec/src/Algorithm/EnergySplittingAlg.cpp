#ifndef _ENERGYSPLITTING_ALG_C
#define _ENERGYSPLITTING_ALG_C

#include "Algorithm/EnergySplittingAlg.h"

EnergySplittingAlg::EnergySplittingAlg(){

} 

void EnergySplittingAlg::Settings::SetInitialValue(){
  Sth_split = -1;
  Eth_SeedAbs = 0.005; //0.1MeV
  Eth_ShowerAbs = 0.0;
  Eth_ClusAbs = 0.0;
  Eth_SeedWithNeigh = 0.;
  Eth_SeedWithTot = 0.;
  Eth_ShowerWithTot = 0.;
  Eth_ClusterWithTot = 0.;
  Etot=0;
  isForced = false;
}

void EnergySplittingAlg::Settings::SetValues( double _Sth_split, double _Eth_seedAbs, double _Eth_showerAbs, double _Eth_clusAbs,
                                              double _Eth_seedwnei, double _Eth_seedwtot, double _Eth_showerwtot, double _Eth_cluswtot, bool _isforced ){
  Sth_split          = _Sth_split; 
  Eth_SeedAbs        = _Eth_seedAbs; 
  Eth_ShowerAbs      = _Eth_showerAbs; 
  Eth_ClusAbs        = _Eth_clusAbs;
  Eth_SeedWithNeigh  = _Eth_seedwnei;
  Eth_SeedWithTot    = _Eth_seedwtot; 
  Eth_ShowerWithTot  = _Eth_showerwtot;
  Eth_ClusterWithTot = _Eth_cluswtot;
  Etot               = 0; 
  isForced           = _isforced;

}

StatusCode EnergySplittingAlg::Initialize(){

  return StatusCode::SUCCESS;
}

StatusCode EnergySplittingAlg::RunAlgorithm( EnergySplittingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol ){

  settings = m_settings;
  std::vector<CRDEcalEDM::CRDCaloLayer> m_LayerCol; m_LayerCol.clear();
  int Nblocks = m_datacol.BlockVec.size();
  for(int ib=0;ib<Nblocks;ib++){
    CRDEcalEDM::CRDCaloLayer m_layer = ClusterinLayer( m_datacol.BlockVec[ib] );
    m_LayerCol.push_back(m_layer);
  }
  m_datacol.LayerCol = m_LayerCol;
  return StatusCode::SUCCESS;

}


CRDEcalEDM::CRDCaloLayer EnergySplittingAlg::ClusterinLayer( CRDEcalEDM::CRDCaloBlock& m_blocks){
  CRDEcalEDM::CRDCaloLayer m_layer; m_layer.Clear();

  m_layer.barXCol = m_blocks.getBarXCol();
  m_layer.barYCol = m_blocks.getBarYCol();

  bool f_succX = true;  bool f_succY = true;
  if( !ClusteringWithCandidate( m_layer.barXCol, m_layer.barClusXCol, m_blocks.getCandidateCol()) ) {  m_layer.barClusXCol = Clustering( m_layer.barXCol ); f_succX = false; }
  if( !ClusteringWithCandidate( m_layer.barYCol, m_layer.barClusYCol, m_blocks.getCandidateCol()) ) {  m_layer.barClusYCol = Clustering( m_layer.barYCol ); f_succY = false; }


  for(int i=0;i<m_layer.barClusXCol.size();i++){
    std::vector<CRDEcalEDM::CRDCaloBarShower> showers; showers.clear();
    //---WIP: Split the cluster with expected showers (position & energy). 
    //if(settings.isForced && f_succX) showers = ClusterSplitting( m_layer.barClusXCol[i], m_blocks.getCandidateCol() );
    //else showers = ClusterSplitting(m_layer.barClusXCol[i]);
    showers = ClusterSplitting(m_layer.barClusXCol[i]);
    if(showers.size()==0) continue;
    m_layer.barShowerXCol.insert(m_layer.barShowerXCol.end(), showers.begin(), showers.end());
  }


  for(int i=0;i<m_layer.barClusYCol.size();i++){
    std::vector<CRDEcalEDM::CRDCaloBarShower> showers; showers.clear(); 
    //if(settings.isForced && f_succY) showers = ClusterSplitting( m_layer.barClusYCol[i], m_blocks.getCandidateCol() );
    //else showers = ClusterSplitting(m_layer.barClusYCol[i]);
    showers = ClusterSplitting(m_layer.barClusYCol[i]);
    if(showers.size()==0) continue;
    m_layer.barShowerYCol.insert(m_layer.barShowerYCol.end(), showers.begin(), showers.end());
  }

  return m_layer;
}


std::vector<CRDEcalEDM::CRDCaloBarCluster> EnergySplittingAlg::Clustering(std::vector<CRDEcalEDM::CRDCaloBar>& barCol){

  for(int i=0;i<barCol.size();i++) settings.Etot+=barCol[i].getEnergy();
  std::sort(barCol.begin(), barCol.end());

  std::vector<CRDEcalEDM::CRDCaloBarCluster> m_clusCol; m_clusCol.clear();

  //Neighbor clustering
  for(int i=0;i<barCol.size();i++){
    CRDEcalEDM::CRDCaloBar iBar = barCol[i];

    if((m_clusCol.size()!=0) && (m_clusCol[m_clusCol.size()-1].isNeighbor(iBar)) ){
      m_clusCol[m_clusCol.size()-1].addBar( iBar );
      continue;
    }

    CRDEcalEDM::CRDCaloBarCluster clus; clus.Clear();
    clus.addBar(iBar);
    m_clusCol.push_back(clus);
  }


  //Find seeds in cluster. 
  std::vector<CRDEcalEDM::CRDCaloBarCluster> out_clusCol; out_clusCol.clear();
  for(int i=0;i<m_clusCol.size();i++){
    CRDEcalEDM::CRDCaloBarCluster iclus = m_clusCol[i];
    if( (iclus.getE()/settings.Etot < settings.Eth_ClusterWithTot) || iclus.getE()<settings.Eth_ClusAbs ) continue;
    std::vector<CRDEcalEDM::CRDCaloBar> m_seedVec = findSeeds(iclus);
    iclus.setSeeds(m_seedVec);
    iclus.getScndMoment();
    out_clusCol.push_back(iclus);
  }

  //Merge Clusters without seed
  for(int i=0; i<out_clusCol.size(); i++){
    if(out_clusCol[i].getNseeds()==0){  
      MergeToClosestCluster( out_clusCol[i], out_clusCol );
      out_clusCol.erase(out_clusCol.begin()+i);
      i--;
    }
  }

  return out_clusCol;
}



bool EnergySplittingAlg::ClusteringWithCandidate( std::vector<CRDEcalEDM::CRDCaloBar>& barCol, 
                                                  std::vector<CRDEcalEDM::CRDCaloBarCluster>& outClus,  
                                                  const std::vector<CRDEcalEDM::CRDShowerCandidate>& m_CandiCol )
{

  //If no expected EM shower, return false, and clustering bars with normal way. 

  if(m_CandiCol.size()==0) return false;
  if(barCol.size()==0) return false; 

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

  std::vector<CRDEcalEDM::CRDCaloBarCluster> out_clusCol; out_clusCol.clear();
  for(int i=0;i<m_clusCol.size();i++){
    CRDEcalEDM::CRDCaloBarCluster iclus = m_clusCol[i];
    if( (iclus.getE()/settings.Etot < settings.Eth_ClusterWithTot) || iclus.getE()<settings.Eth_ClusAbs ) continue;
    std::vector<CRDEcalEDM::CRDCaloBar> m_seedVec; m_seedVec.clear();

    m_seedVec = findSeeds(iclus, m_CandiCol, true);
    iclus.setSeeds(m_seedVec);
    iclus.getScndMoment();
    out_clusCol.push_back(iclus);
  }

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


std::vector<CRDEcalEDM::CRDCaloBar> EnergySplittingAlg::GetMatchedSeeds( std::vector<CRDEcalEDM::CRDCaloBar>& seeds, const std::vector<CRDEcalEDM::CRDShowerCandidate>& m_CandiCol ) {

  std::vector<CRDEcalEDM::CRDCaloBar> m_matchedseeds; m_matchedseeds.clear();

  for(int is=0;is<seeds.size(); is++){
  for(int j=0; j<m_CandiCol.size(); j++){
    if( seeds[is].getSlayer()==0){
      if(    fabs( seeds[is].getPosition().z()-m_CandiCol[j].ExpPos.z() ) < 20  // Position difference within 20mm
//          && fabs( seeds[is].getEnergy()-m_CandiCol[j].ExpEseed )/m_CandiCol[j].ExpEseed < 0.2
      ){ m_matchedseeds.push_back(seeds[is]); break; }
    }
    else if (seeds[is].getSlayer()==1){
      double rotAngle = -seeds[is].getModule()*PI/4.;
      TVector3 p_exps = m_CandiCol[j].ExpPos; 
      TVector3 p_seed(0,0,0);
      p_seed.SetXYZ( seeds[is].getPosition().x(), seeds[is].getPosition().y(), seeds[is].getPosition().z() );

      p_exps.RotateZ(rotAngle); 
      p_seed.RotateZ(rotAngle);
      if( (p_exps-p_seed).Mag() < 20  
//          && fabs( seeds[is].getEnergy()-m_CandiCol[j].ExpEseed )/m_CandiCol[j].ExpEseed < 0.2  
      ){ m_matchedseeds.push_back(seeds[is]); break; }
    }
  }}
  return m_matchedseeds;
}


std::vector<CRDEcalEDM::CRDShowerCandidate> EnergySplittingAlg::CandidateInCluster( CRDEcalEDM::CRDCaloBarCluster& cluster, const std::vector<CRDEcalEDM::CRDShowerCandidate>& m_CandiCol){

  std::vector<CRDEcalEDM::CRDShowerCandidate> out_expsh; out_expsh.clear(); 
  double xmin, xmax, ymin, ymax, zmin, zmax; 
  cluster.getGlobalRange(xmin, ymin, zmin, xmax, ymax, zmax);
//printf("cluster range: (%.2f, %.2f), (%.2f, %.2f), (%.2f, %.2f) \n", xmin, xmax, ymin, ymax, zmin, zmax );

  for(int i=0;i<m_CandiCol.size();i++){
    TVector3 pos = m_CandiCol[i].ExpPos; 
//printf("Exp shower position: (%.2f, %.2f, %.2f) \n", pos.x(), pos.y(), pos.z());
    if( pos.x()>=xmin && pos.x()<=xmax && 
        pos.y()>=ymin && pos.y()<=ymax &&
        pos.z()>=zmin && pos.z()<=zmax
      )
    out_expsh.push_back(m_CandiCol[i]);
  }

  return out_expsh;
}


std::vector<CRDEcalEDM::CRDCaloBar>  EnergySplittingAlg::findSeeds( CRDEcalEDM::CRDCaloBarCluster& m_cluster){
  const std::vector<CRDEcalEDM::CRDShowerCandidate> emptyCol; 
  return findSeeds(m_cluster, emptyCol, false);
}

std::vector<CRDEcalEDM::CRDCaloBar>  EnergySplittingAlg::findSeeds( CRDEcalEDM::CRDCaloBarCluster& m_cluster, const std::vector<CRDEcalEDM::CRDShowerCandidate>& m_CandiCol, bool f_force ){


  std::vector<CRDEcalEDM::CRDCaloBar> m_seeds; m_seeds.clear();
  std::vector<CRDEcalEDM::CRDCaloBar> m_localMax; m_localMax.clear();
  for(int i=0;i<m_cluster.getBars().size();i++){
    CRDEcalEDM::CRDCaloBar ibar = m_cluster.getBars()[i];
    std::vector<CRDEcalEDM::CRDCaloBar> m_neighbor = getNeighbors(m_cluster, ibar);
    if(m_neighbor.size()==0){
      //std::cout<<"WARNING: An isolated bar without neighbors! "<<std::endl;
      if(ibar.getEnergy()>settings.Eth_SeedAbs) m_seeds.push_back(ibar);
      continue;
    }

    bool isLocalMax=true;
    bool isIso;
    double Eneigh=0;
    for(int j=0;j<m_neighbor.size();j++){
      if(m_neighbor[j].getEnergy()>ibar.getEnergy()) isLocalMax=false;
      Eneigh += m_neighbor[j].getEnergy();
    }
    if(isLocalMax) m_localMax.push_back(ibar);
    isIso = (ibar.getEnergy()/(ibar.getEnergy()+Eneigh))>settings.Eth_SeedWithNeigh;
    if(ibar.getEnergy()>settings.Eth_SeedAbs && isLocalMax && isIso) m_seeds.push_back(ibar);
  }

  //Run without expected shower: return seeds.
  if(!f_force) return m_seeds;

  //Run with expected shower: 
  std::vector<CRDEcalEDM::CRDCaloBar> m_forcedSeeds; m_forcedSeeds.clear();

  if( m_CandiCol.size()==0 || m_localMax.size()==0 ) return m_forcedSeeds;

  m_forcedSeeds = GetMatchedSeeds(m_localMax, m_CandiCol);
  return m_forcedSeeds;

}



std::vector<CRDEcalEDM::CRDCaloBar>  EnergySplittingAlg::getNeighbors(CRDEcalEDM::CRDCaloBarCluster& m_cluster, CRDEcalEDM::CRDCaloBar& seed){

  std::vector<CRDEcalEDM::CRDCaloBar> m_neighbor;
  for(int i=0;i<m_cluster.getBars().size();i++){

    if( seed.isNeighbor(m_cluster.getBars()[i]) ) m_neighbor.push_back(m_cluster.getBars()[i]);
  }
  if(m_neighbor.size()>2) std::cout<<"WARNING: more than 2 hits in neighborCol!!"<<std::endl;
  //if(m_neighbor.size()==0) std::cout<<"WARNING: Can not find neighborCol!!"<<std::endl;

  return m_neighbor;
}



std::vector<CRDEcalEDM::CRDCaloBarShower>  EnergySplittingAlg::ClusterSplitting( CRDEcalEDM::CRDCaloBarCluster& m_cluster ){

  std::vector<CRDEcalEDM::CRDCaloBarShower> m_showers; m_showers.clear();

  //No seed in cluster: return empty vector. 
  if(m_cluster.getNseeds()==0) { std::cout<<"WARNING: Still have no-seed cluster!!"<<std::endl; return m_showers; }

  //1 seed or second moment less than threshold: Not split. Turn cluster to shower and return
  else if(m_cluster.getNseeds()<2 || m_cluster.getScndMoment()<settings.Sth_split){
    CRDEcalEDM::CRDCaloBarShower shower; shower.Clear();
    shower.setBars(m_cluster.getBars());
    if(m_cluster.getNseeds()!=0) shower.setSeed(m_cluster.getSeeds()[0]);
    if( (shower.getE()/settings.Etot > settings.Eth_ShowerWithTot) && shower.getE()>settings.Eth_ShowerAbs  ) m_showers.push_back(shower);
    return m_showers;
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
  do{

    for(int ibar=0;ibar<m_cluster.getBars().size();ibar++){
      double Eexp[Nshower];
      double Eexp_tot=0;
      for(int is=0;is<Nshower;is++){ Eexp[is] = Eseed[is]*GetShowerProfile(m_cluster.getBars()[ibar].getPosition(), SeedPos[is] ); Eexp_tot+= Eexp[is];}
      for(int is=0;is<Nshower;is++) weight[ibar][is] = Eexp[is]/Eexp_tot;
    }

    dd4hep::Position SeedPos_prev[Nshower];
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
  if(iter>=20) std::cout<<"WARNING: Iteration time larger than 20! Might not converge!"<<std::endl;

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

  return m_showers;

}

std::vector<CRDEcalEDM::CRDCaloBarShower>  EnergySplittingAlg::ClusterSplitting( CRDEcalEDM::CRDCaloBarCluster& m_cluster, 
                                                                                              const std::vector<CRDEcalEDM::CRDShowerCandidate>& m_candidates )
{
  std::vector<CRDEcalEDM::CRDCaloBarShower> m_showers; m_showers.clear();
  return m_showers;
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
    if(dis<minD){ minD = dis; index=icl; }
  }
  if(index<0) return false; 

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

