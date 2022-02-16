/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
// Unit in code: mm, ns. 
// Digitization for Fan crystal ECAL design from Huaqiao Zhang.

#ifndef FAN_ECAL_DIGI_ALG_C
#define FAN_ECAL_DIGI_ALG_C

#include "FanEcalDigiAlg.h"

#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Vector3f.h"
#include "edm4hep/Cluster.h"

#include "DD4hep/Detector.h"
#include <DD4hep/Objects.h>
#include <DDRec/CellIDPositionConverter.h>

#include "TVector3.h"
#include <math.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <map>

#define C 299.79  // unit: mm/ns
#define PI 3.141592653
using namespace std;
using namespace dd4hep;

DECLARE_COMPONENT( FanEcalDigiAlg )

FanEcalDigiAlg::FanEcalDigiAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{
  
	// Input collections
	declareProperty("SimCaloHitCollection", r_SimCaloCol, "Handle of the Input SimCaloHit collection");
 	declareProperty("MCParticleCol", m_mcParCol, "MCParticle collection (input)");
	
	// Output collections
	declareProperty("CaloHitCollection", w_DigiCaloCol, "Handle of Digi CaloHit collection");
	declareProperty("CaloAssociationCollection", w_CaloAssociationCol, "Handle of CaloAssociation collection");
	
   
}

StatusCode FanEcalDigiAlg::initialize()
{

	std::string s_outfile = _filename;
	m_wfile = new TFile(s_outfile.c_str(), "recreate");
	t_SimStep = new TTree("SimStep", "SimStep");
	t_SimBar = new TTree("SimBarHit", "SimBarHit");
	t_ClusBar = new TTree("ClusBar", "ClusBar");
        t_RecClus = new TTree("RecClus", "RecClus");
//	t_MCdata = new TTree("MCdata", "MCdata");
//	t_MCep = new TTree("MCep", "MCep");
	t_SimStep->Branch("step_x", &m_step_x);
	t_SimStep->Branch("step_y", &m_step_y);
	t_SimStep->Branch("step_z", &m_step_z);
	t_SimStep->Branch("stepBar_x", &m_stepBar_x);
	t_SimStep->Branch("stepBar_y", &m_stepBar_y);
	t_SimStep->Branch("stepBar_z", &m_stepBar_z);
	t_SimStep->Branch("step_E", &m_step_E);
	t_SimStep->Branch("step_T1", &m_step_T1);
	t_SimStep->Branch("step_T2", &m_step_T2);
	t_SimBar->Branch("simBar_x", &m_simBar_x);
	t_SimBar->Branch("simBar_y", &m_simBar_y);
	t_SimBar->Branch("simBar_z", &m_simBar_z);
	t_SimBar->Branch("simBar_T1", &m_simBar_T1);
	t_SimBar->Branch("simBar_T2", &m_simBar_T2);
	t_SimBar->Branch("simBar_Q1", &m_simBar_Q1);
	t_SimBar->Branch("simBar_Q2", &m_simBar_Q2);
	t_SimBar->Branch("simBar_module", &m_simBar_module);
	t_SimBar->Branch("simBar_crystal", &m_simBar_crystal);
	t_SimBar->Branch("simBar_id", &m_simBar_id);
//=================bars in cluster=============================
	t_ClusBar->Branch("clusBar_x", &m_clusBar_x);
	t_ClusBar->Branch("clusBar_y", &m_clusBar_y);
	t_ClusBar->Branch("clusBar_z", &m_clusBar_z);
	t_ClusBar->Branch("clusBar_T1", &m_clusBar_T1);
	t_ClusBar->Branch("clusBar_T2", &m_clusBar_T2);
	t_ClusBar->Branch("clusBar_Q1", &m_clusBar_Q1);
	t_ClusBar->Branch("clusBar_Q2", &m_clusBar_Q2);
	t_ClusBar->Branch("clusBar_module", &m_clusBar_module);
	t_ClusBar->Branch("clusBar_crystal", &m_clusBar_crystal);
	t_ClusBar->Branch("nclusBar", &m_nclusBar);

  t_RecClus->Branch("Ncluster", &m_Ncluster);
  t_RecClus->Branch("clus_id", &m_clus_id);
  t_RecClus->Branch("clus_phi", &m_clus_phi);
  t_RecClus->Branch("clus_phi_start", &m_clus_phi_start);
  t_RecClus->Branch("clus_E", &m_clus_E);
  t_RecClus->Branch("clus_Z", &m_clus_Z);
  t_RecClus->Branch("clus_aphi", &m_clus_aphi);
  t_RecClus->Branch("clus_Nbars", &m_clus_Nbars);
  t_RecClus->Branch("clus_chi2", &m_clus_chi2);
  t_RecClus->Branch("clus_alpha", &m_clus_alpha);
  t_RecClus->Branch("clus_beta", &m_clus_beta);
//===================MC data====================================
  	t_RecClus->Branch("MCendpoint_x", &m_MCendpoint_x);
	t_RecClus->Branch("MCendpoint_y", &m_MCendpoint_y);
	t_RecClus->Branch("MCendpoint_z", &m_MCendpoint_z);
	t_RecClus->Branch("MCendpoint_phi", &m_MCendpoint_phi);
        t_RecClus->Branch("MCendpoint_theta", &m_MCendpoint_theta);
	t_RecClus->Branch("MCendpoint_R", &m_MCendpoint_R);
	t_RecClus->Branch("MCmomentum", &m_MCmomentum);
	t_RecClus->Branch("MCmomentum_x", &m_MCmomentum_x);
	t_RecClus->Branch("MCmomentum_y", &m_MCmomentum_y);
	t_RecClus->Branch("MCmomentum_z", &m_MCmomentum_z);

	std::cout<<"FanEcalDigiAlg::m_scale="<<m_scale<<std::endl;
	m_geosvc = service<IGeomSvc>("GeomSvc");
	if ( !m_geosvc )  throw "FanEcalDigiAlg :Failed to find GeoSvc ...";
	dd4hep::Detector* m_dd4hep = m_geosvc->lcdd();
	if ( !m_dd4hep )  throw "FanEcalDigiAlg :Failed to get dd4hep::Detector ...";
	m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep);
   m_decoder = m_geosvc->getDecoder(_readout);
	if (!m_decoder) {
		error() << "Failed to get the decoder. " << endmsg;
		return StatusCode::FAILURE;
	}
	
	rndm.SetSeed(_seed);
	std::cout<<"FanEcalDigiAlg::initialize"<<std::endl;
	return GaudiAlgorithm::initialize();
}

StatusCode FanEcalDigiAlg::execute()
{
	if(_nEvt==0) std::cout<<"FanEcalDigiAlg::execute Start"<<std::endl;
	std::cout<<"Processing event: "<<_nEvt<<std::endl;
   if(_nEvt<_Nskip){ _nEvt++; return StatusCode::SUCCESS; }

	Clear();

  //=========Readin SimHit collections and create output collections============//
 	const edm4hep::SimCalorimeterHitCollection* SimHitCol =  r_SimCaloCol.get();

	edm4hep::CalorimeterHitCollection* caloVec = w_DigiCaloCol.createAndPut(); //Not used beacuse present edm4hep::CalorimeterHitCollection can not handle double-side (Q,T) readout. 
	edm4hep::MCRecoCaloAssociationCollection* caloAssoVec = w_CaloAssociationCol.createAndPut();
 	std::vector<edm4hep::SimCalorimeterHit> m_simhitCol; m_simhitCol.clear();

//===========MC data===================
 
   auto mcCol = m_mcParCol.get();
   int m_nParticles = 0;
   std::cout<<mcCol->size()<<std::endl; 
   for ( auto p : *mcCol ) {

	m_MCmomentum.push_back(p.getEnergy());		
     	m_MCmomentum_x.push_back(p.getMomentum().x);
	m_MCmomentum_y.push_back(p.getMomentum().y);
	m_MCmomentum_z.push_back(p.getMomentum().z);
	m_MCendpoint_x.push_back(p.getEndpoint().x);
	m_MCendpoint_y.push_back(p.getEndpoint().y);
	m_MCendpoint_z.push_back(p.getEndpoint().z);
    	m_MCendpoint_phi.push_back(getPhi(p.getEndpoint().x,p.getEndpoint().y));
	m_MCendpoint_theta.push_back(getTheta(p.getEndpoint().x,p.getEndpoint().y,p.getEndpoint().z));
	m_MCendpoint_R.push_back(getR(p.getEndpoint().x,p.getEndpoint().y,p.getEndpoint().z));
    }
   // t_MCdata->Fill();
  //================Digitization=================
  if(SimHitCol == 0) 
  {
     std::cout<<"not found SimCalorimeterHitCollection"<< std::endl;
     return StatusCode::SUCCESS;
  }
  if(_Debug>=1) std::cout<<"digi, input sim hit size="<< SimHitCol->size() <<std::endl;

	double totE_bar=0;

	//Merge input simhit(steps) to real simhit(bar).
	m_simhitCol = MergeHits(*SimHitCol);
	if(_Debug>=1) std::cout<<"Finish Hit Merge, with Nhit: "<<m_simhitCol.size()<<std::endl;

  std::vector<CaloBar> HitBarVec; HitBarVec.clear(); 
	//Loop in SimHit, digitalize SimHit to DigiBar
	for(int i=0;i<m_simhitCol.size();i++){

		edm4hep::SimCalorimeterHit SimHit = m_simhitCol.at(i);
		if(SimHit.getEnergy()<_Eth) continue;


		unsigned long long id = SimHit.getCellID();
		CaloBar hitbar;
		hitbar.setcellID( id);
		hitbar.setcellID(	m_decoder->get(id, "system"), 
											m_decoder->get(id, "module"), 
											m_decoder->get(id, "crystal"));

		dd4hep::Position hitpos = m_cellIDConverter->position(id);
		dd4hep::Position barpos(10*hitpos.x(), 10*hitpos.y(), 10*hitpos.z());	//cm to mm.
		hitbar.setPosition(barpos);
		//hitbar.T1 = 99999; hitbar.T2 = 99999;
		//if(_Debug>=2) std::cout<<"SimHit contribution size: "<<SimHit.contributions_size()<<std::endl;


		std::vector<CaloStep> DigiLvec; DigiLvec.clear();
		std::vector<CaloStep> DigiRvec; DigiRvec.clear();
		double totQ1 = 0;
		double totQ2 = 0;

		//Loop in all SimHitContribution(G4Step). 
		for(int iCont=0; iCont < SimHit.contributions_size(); ++iCont){
			edm4hep::ConstCaloHitContribution conb = SimHit.getContributions(iCont);
			if( !conb.isAvailable() ) { std::cout<<"FanEcalDigiAlg  Can not get SimHitContribution: "<<iCont<<std::endl; continue;}

			double en = conb.getEnergy();
			if(en == 0) continue;
      
			dd4hep::Position steppos(conb.getStepPosition().x, conb.getStepPosition().y, conb.getStepPosition().z);
			dd4hep::Position rpos = steppos-hitbar.getPosition();
      
			m_step_x.push_back(steppos.x());
			m_step_y.push_back(steppos.y());
			m_step_z.push_back(steppos.z());
			m_step_E.push_back(en);
			m_stepBar_x.push_back(hitbar.getPosition().x());
			m_stepBar_y.push_back(hitbar.getPosition().y());
			m_stepBar_z.push_back(hitbar.getPosition().z());

			if(_Debug>=3){
				cout<<"Cell Pos: "<<hitbar.getPosition().x()<<'\t'<<hitbar.getPosition().y()<<'\t'<<hitbar.getPosition().z()<<endl;
				cout<<"step pos: "<<steppos.x()<<'\t'<<steppos.y()<<'\t'<<steppos.z()<<endl;
				cout<<"Relative pos: "<<rpos.x()<<'\t'<<rpos.y()<<'\t'<<rpos.z()<<endl;
				cout<<"Cell: "<<hitbar.getModule()<<"  "<<hitbar.getCrystal()<<endl;
			}
      
			//Get digitalized signal(Q1, Q2, T1, T2) from step
			//Define: 1 is left, 2 is right, clockwise direction in phi. 
      int sign = (steppos.Mag2()-hitbar.getPosition().Mag2()>0) ? -1 : 1; // Inner: sign=-1; outer: sign=1; 
      double Qi_inner = en*exp( -( Lbar/2 - sign*sqrt(rpos.Mag2()) )/Latt );
      double Qi_outer = en*exp( -( Lbar/2 + sign*sqrt(rpos.Mag2()) )/Latt );



			double Ti_inner = -1; int looptime=0;
			while(Ti_inner<0){ 
				Ti_inner = Tinit + rndm.Gaus(nMat*(Lbar/2 - sign*sqrt(rpos.Mag2()))/C, Tres); 
				looptime++;
				if(looptime>500){ std::cout<<"ERROR: Step "<<iCont<<" can not get a positive left-side time!"<<std::endl; break;}
			}
			if(looptime>500) continue;		
			double Ti_outer = -1; looptime=0;
			while(Ti_outer<0){ 
				Ti_outer = Tinit + rndm.Gaus(nMat*(Lbar/2 + sign*sqrt(rpos.Mag2()))/C, Tres); 
				looptime++;
            if(looptime>500){ std::cout<<"ERROR: Step "<<iCont<<" can not get a positive right-side time!"<<std::endl; break;}
			}
			if(looptime>500) continue;		


			m_step_T1.push_back(Ti_inner);
			m_step_T2.push_back(Ti_outer);
			totQ1 += Qi_inner;
			totQ2 += Qi_outer;
	
			CaloStep stepoutL, stepoutR;
			stepoutL.setQ(Qi_inner); stepoutL.setT(Ti_inner);
			stepoutR.setQ(Qi_outer); stepoutR.setT(Ti_outer);
			DigiLvec.push_back(stepoutL);
			DigiRvec.push_back(stepoutR);


		}

		//Time digitalization
		//if(_Debug>=2) std::cout<<"Time Digitalize: time at Q >"<<_Qthfrac<<"*totQ"<<std::endl;
		std::sort(DigiLvec.begin(), DigiLvec.end());
		std::sort(DigiRvec.begin(), DigiRvec.end());
		double thQ1=0;
		double thQ2=0;
		double thT1, thT2; 
		for(int iCont=0;iCont<DigiLvec.size();iCont++){
			thQ1 += DigiLvec[iCont].getQ();
			if(thQ1>totQ1*_Qthfrac){ 
				thT1 = DigiLvec[iCont].getT(); 
				if(_Debug>=3) std::cout<<"Get T1 at index: "<<iCont<<std::endl;
				break;
			}
		}
		for(int iCont=0;iCont<DigiRvec.size();iCont++){
			thQ2 += DigiRvec[iCont].getQ();
			if(thQ2>totQ2*_Qthfrac){ 
				thT2 = DigiRvec[iCont].getT(); 
				if(_Debug>=3) std::cout<<"Get T2 at index: "<<iCont<<std::endl;
				break;
			}
		}
		hitbar.setT(thT1, thT2);
		hitbar.setQ(totQ1, totQ2);
    HitBarVec.push_back(hitbar);

		totE_bar+=hitbar.getEnergy();
		//unsigned long int blockID = coder(hitbar);
		//DigiBlocks[blockID].push_back(hitbar);
		
		m_simBar_x.push_back(hitbar.getPosition().x());
		m_simBar_y.push_back(hitbar.getPosition().y());
		m_simBar_z.push_back(hitbar.getPosition().z());
		m_simBar_Q1.push_back(hitbar.getQ1());
		m_simBar_Q2.push_back(hitbar.getQ2());
		m_simBar_T1.push_back(hitbar.getT1());
		m_simBar_T2.push_back(hitbar.getT2());
		m_simBar_module.push_back(hitbar.getModule());
		m_simBar_crystal.push_back(hitbar.getCrystal());
		m_simBar_id.push_back(_nEvt);
	}
	t_SimStep->Fill();
	t_SimBar->Fill();
	if(_Debug>=1) std::cout<<"End Loop: Bar Digitalization!"<<std::endl;
	std::cout<<"Total Bar Energy: "<<totE_bar<<std::endl;
  //=======Digitization end, results are stored in HitBarVec=======


  //================Reconstruction================
  std::vector<CaloCluster> ClusVec; ClusVec.clear();
  NeighborClustering(HitBarVec, ClusVec);

  double tot_clus_E = 0;
  double Eth_bar = totE_bar/14.;
  m_Ncluster = ClusVec.size();
  cout<<"m_Ncluster is:"<<m_Ncluster<<endl;
  for(int ic=0; ic<m_Ncluster; ic++){
    m_clus_phi_start.push_back(ClusVec[ic].getPhiStart());
    m_clus_phi.push_back(ClusVec[ic].getPhiStartFit());
    m_clus_Z.push_back(ClusVec[ic].getZ());
    m_clus_aphi.push_back(ClusVec[ic].getAvePhi());
    m_clus_E.push_back(ClusVec[ic].getEnergy());
    m_clus_chi2.push_back(ClusVec[ic].getFitChi2());
    m_clus_alpha.push_back(ClusVec[ic].getFitAlpha());
    m_clus_beta.push_back(ClusVec[ic].getFitBeta());
    m_clus_Nbars.push_back(ClusVec[ic].getBars().size());
    m_clus_id.push_back(ic+100*_nEvt); 
    tot_clus_E+=ClusVec[ic].getEnergy();
  //  cout<<"This is cluster ic: "<<ic<<" energy: "<<ClusVec[ic].getEnergy()<<" ; total energy is: "<<tot_clus_E<<endl;
//=====================ADD============================
   int nseedBar = 0;
   vector<CaloBar> n_hitBars;n_hitBars.clear();
    n_hitBars = ClusVec[ic].getBars();
    if(n_hitBars.size()==1 && n_hitBars[0].getEnergy()>Eth_bar){
	vector<CaloBar> neighborBar; neighborBar.clear();
	neighborBar.push_back(n_hitBars[0]);	
	continue;
    }/////
    for(int ibar=0;ibar<n_hitBars.size();ibar++){
	if(ClusVec[ic].getEnergy()>totE_bar*0.4){
		m_clusBar_x.push_back(n_hitBars[ibar].getPosition().x());
                m_clusBar_y.push_back(n_hitBars[ibar].getPosition().y());
                m_clusBar_z.push_back(n_hitBars[ibar].getPosition().z());
                m_clusBar_Q1.push_back(n_hitBars[ibar].getQ1());
                m_clusBar_Q2.push_back(n_hitBars[ibar].getQ2());
                m_clusBar_T1.push_back(n_hitBars[ibar].getT1());
                m_clusBar_T2.push_back(n_hitBars[ibar].getT2());
                m_clusBar_module.push_back(n_hitBars[ibar].getModule());
                m_clusBar_crystal.push_back(n_hitBars[ibar].getCrystal());
                m_nclusBar.push_back(ic+100*_nEvt);
	}
	if(n_hitBars[ibar].getEnergy()>Eth_bar){
	   vector<CaloBar> neighborBar; neighborBar.clear();
	   neighborBar = ClusVec[ic].getNeighbor(n_hitBars[ibar]);
	   if(isMaxEnergy(neighborBar,n_hitBars[ibar])){
	    // ClusVec[ic].AddBar(n_hitBars[ibar]);	     
	     nseedBar+=1;
	     cout<<"here is seed bar: "<<n_hitBars[ibar].getModule()<<" "<<n_hitBars[ibar].getCrystal()<<" "<<n_hitBars[ibar].getEnergy()<<endl;
	   }
        } 
//cout<<"nseedbar : "<<nseedBar<<endl;
    }
    if(nseedBar==2){
	ClusVec[ic].Fit2Profile();
    }
//====================================================
  }
  t_RecClus->Fill();
  t_ClusBar->Fill();
  
  _nEvt ++ ;
  return StatusCode::SUCCESS;
}

StatusCode FanEcalDigiAlg::finalize()
{
	m_wfile->cd();
	t_SimStep->Write();
	t_SimBar->Write();
	t_ClusBar->Write();
  t_RecClus->Write();
//	t_MCdata->Write(); 
	m_wfile->Close();

  info() << "Processed " << _nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}

std::vector<edm4hep::SimCalorimeterHit> FanEcalDigiAlg::MergeHits(const edm4hep::SimCalorimeterHitCollection& m_col){
	std::vector<edm4hep::SimCalorimeterHit> m_mergedhit;
	m_mergedhit.clear();

	for(int iter=0;iter<m_col.size();iter++){
		edm4hep::SimCalorimeterHit m_step = m_col[iter];
		if(!m_step.isAvailable()){ cout<<"ERROR HIT!"<<endl; continue;}
		if(m_step.getEnergy()==0) continue;
		unsigned long long cellid = m_step.getCellID();
		dd4hep::Position hitpos = m_cellIDConverter->position(cellid);
		edm4hep::Vector3f pos(hitpos.x()*10, hitpos.y()*10, hitpos.z()*10);

		edm4hep::CaloHitContribution conb;
		conb.setEnergy(m_step.getEnergy());
		conb.setStepPosition(m_step.getPosition());
		edm4hep::SimCalorimeterHit m_hit = find(m_mergedhit, cellid);
		if(m_hit.getCellID()==0){
			//m_hit = new edm4hep::SimCalorimeterHit();
			m_hit.setCellID(cellid);
			m_hit.setPosition(pos);
			m_mergedhit.push_back(m_hit);
		}
		m_hit.addToContributions(conb);
		m_hit.setEnergy(m_hit.getEnergy()+m_step.getEnergy());
	}
	return m_mergedhit;
}


edm4hep::SimCalorimeterHit FanEcalDigiAlg::find(std::vector<edm4hep::SimCalorimeterHit>& m_col, unsigned long long& cellid){
   for(int i=0;i<m_col.size();i++){
		edm4hep::SimCalorimeterHit hit=m_col.at(i);
		if(hit.getCellID() == cellid) return hit;
	}
	edm4hep::SimCalorimeterHit hit;
	hit.setCellID(0);
	return hit;
}


StatusCode FanEcalDigiAlg::NeighborClustering(std::vector<CaloBar>& m_barVec, std::vector<CaloCluster>& m_clusVec ){

  if(m_barVec.size()==0) return StatusCode::SUCCESS;

  m_clusVec.clear();
  for(int ib=0; ib<m_barVec.size(); ib++ ){

    bool isNew = true;
    for(int ic=0; ic<m_clusVec.size(); ic++){
      if( m_clusVec[ic].isNeighbor(m_barVec[ib]) ){
        m_clusVec[ic].AddBar(m_barVec[ib]); 
        isNew = false; 
        break; 
    }}
    if(isNew){
      CaloCluster m_clus; m_clus.Clear();
      m_clus.AddBar( m_barVec[ib] );
      m_clusVec.push_back( m_clus );
    }
  }

  return StatusCode::SUCCESS;
}
//===========add=================
bool FanEcalDigiAlg::isMaxEnergy(std::vector<CaloBar>& m_barVec, CaloBar& m_bar){
  double Max_E = 0;
  for(int i=0;i<m_barVec.size();i++){
    if(m_barVec[i].getEnergy()>Max_E){
	Max_E = m_barVec[i].getEnergy();
    }
    Max_E = Max_E;
  } 
  if(m_bar.getEnergy()>=Max_E){
    return true;
  }
  return false;
}

double FanEcalDigiAlg::getPhi(double x, double y){
  double phi=0;
  if(x>0&&y>0) {phi = atan2(y,x)*180./PI;}
  if(x<0&&y>0) {phi = 180.-atan2(abs(y),abs(x))*180./PI;}
  if(x<0&&y<0) {phi = 180.+atan2(abs(y),abs(x))*180./PI;}
  if(x>0&&y<0) {phi = 360.-atan2(abs(y),abs(x))*180./PI;}
  return phi;  
}

double FanEcalDigiAlg::getR(double x, double y, double z){
      
      return sqrt( x*x + y*y + z*z) ; 
}

double FanEcalDigiAlg::getTheta(double x, double y, double z){
  double r = getR(x,y,z);
  if(x == 0.0 && y == 0.0 && z == 0.0){return 0.0;}
  return acos(z/r)*180./PI;
}
//=================================

void FanEcalDigiAlg::Clear(){
	m_step_x.clear();
	m_step_y.clear();
	m_step_z.clear();
	m_step_E.clear();
	m_stepBar_x.clear();
	m_stepBar_y.clear();
	m_stepBar_z.clear();
	m_step_T1.clear();
	m_step_T2.clear();
	m_simBar_x.clear();
	m_simBar_y.clear();
	m_simBar_z.clear();
	m_simBar_T1.clear();
	m_simBar_T2.clear();
	m_simBar_Q1.clear();
	m_simBar_Q2.clear();
	m_simBar_module.clear();
	m_simBar_crystal.clear();
	m_simBar_id.clear();

  m_Ncluster = -99;
  m_clus_id.clear();
  m_clus_phi.clear();
  m_clus_phi_start.clear();
  m_clus_Z.clear();
  m_clus_aphi.clear();
  m_clus_E.clear();
  m_clus_chi2.clear();
  m_clus_alpha.clear();
  m_clus_beta.clear();
  m_clus_Nbars.clear(); 
//==============MC data====================
  m_MCmomentum.clear();
  m_MCmomentum_x.clear();
  m_MCmomentum_y.clear(); 
  m_MCmomentum_z.clear();
  m_MCendpoint_x.clear();
  m_MCendpoint_y.clear();
  m_MCendpoint_z.clear();
  m_MCendpoint_phi.clear();
  m_MCendpoint_theta.clear();
  m_MCendpoint_R.clear();

//==============bars in cluster=============
	m_clusBar_x.clear();
        m_clusBar_y.clear();
        m_clusBar_z.clear();
        m_clusBar_T1.clear();
        m_clusBar_T2.clear();
        m_clusBar_Q1.clear();
        m_clusBar_Q2.clear();
        m_clusBar_module.clear();
        m_clusBar_crystal.clear();
        m_nclusBar.clear();
//==========================================
}

#endif
