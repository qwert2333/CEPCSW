#ifndef _ENERGYSPLITTING_ALG_C
#define _ENERGYSPLITTING_ALG_C

#include "Algorithm/EnergySplittingAlg.h"

EnergySplittingAlg::EnergySplittingAlg(){

} 

void EnergySplittingAlg::Settings::SetInitialValue(){
      Sth_split = -1;
      Eth_SeedAbs = 0.02; //10MeV
      Eth_ShowerAbs = 0.03;
      Eth_ClusAbs = 0.03;
      Eth_SeedWithNeigh = 0.33;
      Eth_SeedWithTot = 0.;
      Eth_ShowerWithTot = 0.01;
      Eth_ClusterWithTot = 0.01;
      Etot=0;
}


StatusCode EnergySplittingAlg::Initialize(){

  return StatusCode::SUCCESS;
}

StatusCode EnergySplittingAlg::RunAlgorithm( EnergySplittingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol ){

  settings = m_settings;
  std::vector<CRDEcalEDM::CRDCaloLayer> m_LayerCol; m_LayerCol.clear();
  int Nblocks = m_datacol.blockVec.size();
  for(int ib=0;ib<Nblocks;ib++){
    CRDEcalEDM::CRDCaloLayer m_layer = ClusterinLayer( m_datacol.blockVec[ib]);
    m_LayerCol.push_back(m_layer);
  }
  m_datacol.LayerCol = m_LayerCol;
  return StatusCode::SUCCESS;

}


CRDEcalEDM::CRDCaloLayer EnergySplittingAlg::ClusterinLayer( CRDEcalEDM::DigiBlock& m_blocks){
  CRDEcalEDM::CRDCaloLayer m_layer;

  for(int i=0;i<m_blocks.size();i++){
    if(m_blocks[i].getSlayer()==0) m_layer.barXCol.push_back(m_blocks[i]);
    else m_layer.barYCol.push_back(m_blocks[i]);
  }

  m_layer.barClusXCol = Clustering(m_layer.barXCol);
  m_layer.barClusYCol = Clustering(m_layer.barYCol);

  for(int i=0;i<m_layer.barClusXCol.size();i++){
    double _Etot=0;
    for(int i=0;i<m_layer.barXCol.size();i++) _Etot += m_layer.barXCol[i].getEnergy();
    if(m_layer.barClusXCol[i].getE()/_Etot < settings.Eth_ClusterWithTot || m_layer.barClusXCol[i].getE()<settings.Eth_ClusAbs ) continue;
    std::vector<CRDEcalEDM::CRDCaloBarShower> showers = ClusterSplitting(m_layer.barClusXCol[i]);
    if(showers.size()==0) continue;
    m_layer.barShowerXCol.insert(m_layer.barShowerXCol.end(), showers.begin(), showers.end());
  }
  for(int i=0;i<m_layer.barClusYCol.size();i++){
    double _Etot=0;
    for(int i=0;i<m_layer.barYCol.size();i++) _Etot += m_layer.barYCol[i].getEnergy();
    if(m_layer.barClusYCol[i].getE()/_Etot < settings.Eth_ClusterWithTot || m_layer.barClusYCol[i].getE()<settings.Eth_ClusAbs ) continue;
    std::vector<CRDEcalEDM::CRDCaloBarShower> showers = ClusterSplitting(m_layer.barClusYCol[i]);
    if(showers.size()==0) continue;
    m_layer.barShowerYCol.insert(m_layer.barShowerYCol.end(), showers.begin(), showers.end());
  }

  return m_layer;
}


std::vector<CRDEcalEDM::CRDCaloBarCluster> EnergySplittingAlg::Clustering(std::vector<CRDEcalEDM::CRDCaloBar>& barCol){

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
    iclus.CRDEcalEDM::CRDCaloBarCluster::sortByPos();
    std::vector<CRDEcalEDM::CRDCaloBar> m_seedVec = findSeeds(iclus);
    iclus.setSeeds(m_seedVec);
    iclus.getScndMoment();
    out_clusCol.push_back(iclus);
  }

  return out_clusCol;

}


std::vector<CRDEcalEDM::CRDCaloBar>  EnergySplittingAlg::findSeeds( CRDEcalEDM::CRDCaloBarCluster& m_cluster ){

  std::vector<CRDEcalEDM::CRDCaloBar> m_seeds;
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
    isIso = (ibar.getEnergy()/(ibar.getEnergy()+Eneigh))>settings.Eth_SeedWithNeigh;
    if(ibar.getEnergy()>settings.Eth_SeedAbs && isLocalMax && isIso) m_seeds.push_back(ibar);
  }

  return m_seeds;

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



std::vector<CRDEcalEDM::CRDCaloBarShower>  EnergySplittingAlg::ClusterSplitting( CRDEcalEDM::CRDCaloBarCluster& m_cluster){

  std::vector<CRDEcalEDM::CRDCaloBarShower> m_showers; m_showers.clear();

  //Not split. Turn cluster to shower and return
  if(m_cluster.getNseeds()<2 || m_cluster.getScndMoment()<settings.Sth_split){
    CRDEcalEDM::CRDCaloBarShower shower;
    shower.setBars(m_cluster.getBars());
    if(m_cluster.getNseeds()!=0) shower.setSeed(m_cluster.getSeeds()[0]);
    if( (shower.getE()/settings.Etot > settings.Eth_ShowerWithTot) && shower.getE()>settings.Eth_ShowerAbs  ) m_showers.push_back(shower);
    return m_showers;
  }

  //Split
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

