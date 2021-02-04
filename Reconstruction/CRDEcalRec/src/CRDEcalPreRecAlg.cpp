#include <iostream>
#include <algorithm>
#include <cmath>
#include "CRDEcalPreRecAlg.h"
#include "CRDEcalDigiEDM.h"
using namespace std;

CRDEcalPreRecAlg::CRDEcalPreRecAlg(){
	Sth_split = -1;
	Eth_SeedAbs = 0.05; //50MeV
	Eth_ShowerAbs = 0.01; 
	Eth_ClusAbs = 0.01;
	Eth_SeedWithNeigh = 0.4;
	Eth_SeedWithTot = 0.15;
	Eth_ShowerWithTot = 0.03;
	Eth_ClusterWithTot = 0.03;
	Etot=0;
}

CRDEcalPreRecAlg::~CRDEcalPreRecAlg(){
}

std::vector<CRDEcalDigiEDM::BarCollection>  CRDEcalPreRecAlg::Bars2Shower( std::vector<CRDEcalDigiEDM::DigiBar>& barCol){

	for(int i=0;i<barCol.size();i++) Etot+=barCol[i].getEnergy();

	std::vector<CRDEcalDigiEDM::BarCollection> m_showers; m_showers.clear();
	std::vector<CRDEcalDigiEDM::BarCluster> m_clusCol; m_clusCol.clear();
	m_clusCol = CRDEcalPreRecAlg::Clustering(barCol); 

	for(int i=0;i<m_clusCol.size();i++){
		std::vector<CRDEcalDigiEDM::BarCollection> showers = ClusterSplitting(m_clusCol[i]);
		if(showers.size()==0) continue;
		m_showers.insert(m_showers.end(), showers.begin(), showers.end());
	}

	return m_showers;
}

std::vector<CRDEcalDigiEDM::BarCluster> CRDEcalPreRecAlg::Clustering(std::vector<CRDEcalDigiEDM::DigiBar>& barCol){

	for(int i=0;i<barCol.size();i++) Etot+=barCol[i].getEnergy();
	std::sort(barCol.begin(), barCol.end());

	std::vector<CRDEcalDigiEDM::BarCluster> m_clusCol; m_clusCol.clear();
	
	for(int i=0;i<barCol.size();i++){
		CRDEcalDigiEDM::DigiBar iBar = barCol[i];

		if((m_clusCol.size()!=0) && (m_clusCol[m_clusCol.size()-1].isNeighbor(iBar)) ){
			m_clusCol[m_clusCol.size()-1].Bars.push_back(iBar);
			continue;
		}

		CRDEcalDigiEDM::BarCluster clus;
		clus.Bars.push_back(iBar);
		m_clusCol.push_back(clus);
	}


	std::vector<CRDEcalDigiEDM::BarCluster> out_clusCol; out_clusCol.clear();
	for(int i=0;i<m_clusCol.size();i++){
		CRDEcalDigiEDM::BarCluster iclus = m_clusCol[i];
		if( (iclus.getE()/Etot < Eth_ClusterWithTot) || iclus.getE()<Eth_ClusAbs ) continue;
		iclus.CRDEcalDigiEDM::BarCluster::sortByPos();
		std::vector<CRDEcalDigiEDM::DigiBar> m_seedVec = findSeeds(iclus);
		iclus.Seeds = m_seedVec;
		iclus.Nseeds = m_seedVec.size();
		iclus.ScndMoment = iclus.getScndMoment();
		out_clusCol.push_back(iclus);
	}

	return out_clusCol;
}


std::vector<CRDEcalDigiEDM::BarCollection> CRDEcalPreRecAlg::ClusterSplitting( CRDEcalDigiEDM::BarCluster& m_cluster){

//cout<<"Nseed: "<<m_cluster.Nseeds<<endl;
//if(m_cluster.Nseeds>=2){
//m_cluster.PrintBars();
//m_cluster.PrintSeeds();
//}

	std::vector<CRDEcalDigiEDM::BarCollection> m_showers; m_showers.clear();

	//Not split. Turn cluster to shower and return
	if(m_cluster.Nseeds<2 || m_cluster.ScndMoment<Sth_split){ 
		CRDEcalDigiEDM::BarCollection shower;
		shower.Bars = m_cluster.Bars;
		if(m_cluster.Nseeds!=0) shower.Seed = m_cluster.Seeds[0];
		if( (shower.getE()/Etot > Eth_ShowerWithTot) && shower.getE()>Eth_ShowerAbs  ) m_showers.push_back(shower);
		return m_showers;
	}

	//Split
	int Nshower = m_cluster.Nseeds;
	int Nbars = m_cluster.Bars.size();
	double Eseed[Nshower] = {0};
	double weight[Nbars][Nshower] = {0};
	dd4hep::Position SeedPos[Nshower];
	for(int is=0;is<Nshower;is++) SeedPos[is] = m_cluster.Seeds[is].position;
	CalculateInitialEseed(m_cluster.Seeds, SeedPos, Eseed);

	bool isConverge = false;
	int iter=0;
	do{
//cout<<"Iterate time: "<<iter<<endl;
//cout<<"x      y      z      E"<<endl;
//for(int i=0;i<Nshower;i++){
//cout<<SeedPos[i].x()<<'\t'<<SeedPos[i].y()<<'\t'<<SeedPos[i].z()<<'\t'<<Eseed[i]<<endl;
//}

		for(int ibar=0;ibar<m_cluster.Bars.size();ibar++){
			double Eexp[Nshower]; 
			double Eexp_tot=0;
			for(int is=0;is<Nshower;is++){ Eexp[is] = Eseed[is]*GetShowerProfile(m_cluster.Bars[ibar].position, SeedPos[is] ); Eexp_tot+= Eexp[is];}
			for(int is=0;is<Nshower;is++) weight[ibar][is] = Eexp[is]/Eexp_tot;
		}

		dd4hep::Position SeedPos_prev[Nshower];
	   for(int is=0;is<Nshower;is++){
			SeedPos_prev[is]=SeedPos[is];
   	   CRDEcalDigiEDM::BarCollection shower;
      	shower.Bars = m_cluster.Bars;
	      for(int ib=0;ib<shower.Bars.size();ib++){
   	      shower.Bars[ib].Q1 = shower.Bars[ib].Q1*weight[ib][is];
      	   shower.Bars[ib].Q2 = shower.Bars[ib].Q2*weight[ib][is];
	      }
			SeedPos[is] = shower.getPos();
			double Emax = -99;
			for(int ib=0;ib<shower.Bars.size();ib++) if(shower.Bars[ib].getEnergy()>Emax) Emax=shower.Bars[ib].getEnergy();
			Eseed[is] = Emax;
   	}

		isConverge=true;
		for(int is=0;is<Nshower;is++) if((SeedPos_prev[is]-SeedPos[is]).Mag2()>2.89){ isConverge=false; break;}
		iter++;
	}
	while(iter<20 && !isConverge);
	if(iter>=20) std::cout<<"WARNING: Iteration time larger than 20! Might not converge!"<<std::endl;

	for(int is=0;is<Nshower;is++){
		CRDEcalDigiEDM::BarCollection shower;
		shower.Bars = m_cluster.Bars;
		for(int ib=0;ib<shower.Bars.size();ib++){
			shower.Bars[ib].Q1 = shower.Bars[ib].Q1*weight[ib][is];
			shower.Bars[ib].Q2 = shower.Bars[ib].Q2*weight[ib][is];
		}
		shower.Energy = shower.getE();
		shower.pos = shower.getPos();
		m_showers.push_back(shower);
	}

	return m_showers;

}


std::vector<CRDEcalDigiEDM::DigiBar> CRDEcalPreRecAlg::findSeeds(CRDEcalDigiEDM::BarCluster& m_cluster){

	std::vector<CRDEcalDigiEDM::DigiBar> m_seeds;
	for(int i=0;i<m_cluster.Bars.size();i++){
		CRDEcalDigiEDM::DigiBar ibar = m_cluster.Bars[i]; 
		std::vector<CRDEcalDigiEDM::DigiBar> m_neighbor = getNeighbors(m_cluster, ibar);
		if(m_neighbor.size()==0){ 
			//std::cout<<"WARNING: An isolated bar without neighbors! "<<std::endl;
			if(ibar.getEnergy()>Eth_SeedAbs) m_seeds.push_back(ibar);
			continue;
		}

		bool isLocalMax=true;
		bool isIso;
		double Eneigh=0;
		for(int j=0;j<m_neighbor.size();j++){
			if(m_neighbor[j].getEnergy()>ibar.getEnergy()) isLocalMax=false;
			Eneigh += m_neighbor[j].getEnergy();
		}
		isIso = (ibar.getEnergy()/(ibar.getEnergy()+Eneigh))>Eth_SeedWithNeigh;
		if(ibar.getEnergy()>Eth_SeedAbs && isLocalMax && isIso) m_seeds.push_back(ibar);
	}

	return m_seeds;
}


std::vector<CRDEcalDigiEDM::DigiBar> CRDEcalPreRecAlg::getNeighbors(CRDEcalDigiEDM::BarCluster& m_cluster, CRDEcalDigiEDM::DigiBar& seed){

	std::vector<CRDEcalDigiEDM::DigiBar> m_neighbor;
	for(int i=0;i<m_cluster.Bars.size();i++){

		if( seed.isNeighbor(m_cluster.Bars[i]) ) m_neighbor.push_back(m_cluster.Bars[i]);
	}
	if(m_neighbor.size()>2) std::cout<<"WARNING: more than 2 hits in neighborCol!!"<<std::endl;
	//if(m_neighbor.size()==0) std::cout<<"WARNING: Can not find neighborCol!!"<<std::endl;


	return m_neighbor;
}

void CRDEcalPreRecAlg::RemoveBars(CRDEcalDigiEDM::BarCluster& m_cluster, CRDEcalDigiEDM::BarCollection& m_barCol){
	for(int i=0;i<m_barCol.Bars.size();i++){
		auto iter = std::remove(m_cluster.Bars.begin(), m_cluster.Bars.end(), m_barCol.Bars[i]);
		m_cluster.Bars.erase(iter, m_cluster.Bars.end());
	}
}

void CRDEcalPreRecAlg::RemoveBars(std::vector<CRDEcalDigiEDM::DigiBar>& m_bars, std::vector<CRDEcalDigiEDM::DigiBar>& m_subbars){
	for(int i=0;i<m_subbars.size();i++){
		auto iter = std::remove(m_bars.begin(), m_bars.end(), m_subbars[i]);
		m_bars.erase(iter, m_bars.end());
	}
}

void CRDEcalPreRecAlg::CalculateInitialEseed(std::vector<CRDEcalDigiEDM::DigiBar>& Seeds, dd4hep::Position* pos, double* Eseed){
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
		for(int j=0;j<Nele;j++) matrixR[i][j] = GetShowerProfile(Seeds[i].position, pos[j]);
	}

	matrixR.Invert();
	TVector sol = matrixR*vecE;

	for(int i=0;i<Nele;i++) Eseed[i] = sol[i];

}

double CRDEcalPreRecAlg::GetShowerProfile(dd4hep::Position& p_bar, dd4hep::Position& p_seed ){
	dd4hep::Position rpos = p_bar-p_seed;
	double dis = sqrt(rpos.Mag2());
	double a1, a2, b1, b2, Rm;
	a1=0.037; a2=0.265; b1=0.101; b2=0.437; Rm=1.868;

	return a1*exp(-b1*dis/Rm)+a2*exp(-b2*dis/Rm);
}

