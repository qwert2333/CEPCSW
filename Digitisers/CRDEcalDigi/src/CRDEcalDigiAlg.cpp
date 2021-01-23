/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
// Unit in code: mm, ns. 
// NOTE: This digitialization highly matches detector geometry CRDEcalBarrel_v01. 
#include "CRDEcalDigiEDM.h"
#include "CRDEcalDigiAlg.h"
#include "CRDEcalPreRecAlg.h"
#include "DigiHitsWithPos.h"
#include "DigiHitsWithMatching.h"
#include "DigiHitsWithMatchingL2.h"

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
DECLARE_COMPONENT( CRDEcalDigiAlg )

CRDEcalDigiAlg::CRDEcalDigiAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{
  
	// Input collections
	declareProperty("SimCaloHitCollection", r_SimCaloCol, "Handle of the Input SimCaloHit collection");
  
	// Output collections
	declareProperty("CaloHitCollection", w_DigiCaloCol, "Handle of Digi CaloHit collection");
	declareProperty("TruthSimCaloHitCollection", w_SimCaloTruth, "Handle of Truth SimHit collection");
	declareProperty("CaloAssociationCollection", w_CaloAssociationCol, "Handle of CaloAssociation collection");
   
}

StatusCode CRDEcalDigiAlg::initialize()
{

	std::string s_outfile = _filename;
	m_wfile = new TFile(s_outfile.c_str(), "recreate");
	t_SimCont = new TTree("SimStep", "SimStep");
	t_SimBar = new TTree("SimBarHit", "SimBarHit");
	t_SimTruth = new TTree("TruthSimHit", "TruthSimHit");
	t_Rec = new TTree("RecHit", "RecHit");
	t_PreRec = new TTree("RecBlock", "RecBlock");
	t_SimCont->Branch("step_x", &m_step_x);
	t_SimCont->Branch("step_y", &m_step_y);
	t_SimCont->Branch("step_z", &m_step_z);
	t_SimCont->Branch("stepBar_x", &m_stepBar_x);
	t_SimCont->Branch("stepBar_y", &m_stepBar_y);
	t_SimCont->Branch("stepBar_z", &m_stepBar_z);
	t_SimCont->Branch("step_E", &m_step_E);
	t_SimCont->Branch("step_T1", &m_step_T1);
	t_SimCont->Branch("step_T2", &m_step_T2);
	t_SimBar->Branch("simBar_x", &m_simBar_x);
	t_SimBar->Branch("simBar_y", &m_simBar_y);
	t_SimBar->Branch("simBar_z", &m_simBar_z);
	t_SimBar->Branch("simBar_T1", &m_simBar_T1);
	t_SimBar->Branch("simBar_T2", &m_simBar_T2);
	t_SimBar->Branch("simBar_Q1", &m_simBar_Q1);
	t_SimBar->Branch("simBar_Q2", &m_simBar_Q2);
	t_SimBar->Branch("simBar_module", &m_simBar_module);
	t_SimBar->Branch("simBar_dlayer", &m_simBar_dlayer);
	t_SimBar->Branch("simBar_part", &m_simBar_part);
	t_SimBar->Branch("simBar_block", &m_simBar_block);
	t_SimBar->Branch("simBar_slayer", &m_simBar_slayer);
	t_SimTruth->Branch("simTruth_x", &m_simTruth_x);
	t_SimTruth->Branch("simTruth_y", &m_simTruth_y);
	t_SimTruth->Branch("simTruth_z", &m_simTruth_z);
	t_SimTruth->Branch("simTruth_E", &m_simTruth_E);
	t_SimTruth->Branch("simTruth_dlayer", &m_simTruth_dlayer);
	t_SimTruth->Branch("simTruth_slayer", &m_simTruth_slayer);
	t_Rec->Branch("Rec_x", &m_Rec_x);
	t_Rec->Branch("Rec_y", &m_Rec_y);
	t_Rec->Branch("Rec_z", &m_Rec_z);
	t_Rec->Branch("Rec_E", &m_Rec_E);

	t_PreRec->Branch("PreRec_Bar0x",&m_PreRec_Bar0x);
	t_PreRec->Branch("PreRec_Bar0y",&m_PreRec_Bar0y);
	t_PreRec->Branch("PreRec_Bar0z",&m_PreRec_Bar0z);
	t_PreRec->Branch("PreRec_Bar0E",&m_PreRec_Bar0E);
	t_PreRec->Branch("PreRec_Bar1x",&m_PreRec_Bar1x);
	t_PreRec->Branch("PreRec_Bar1y",&m_PreRec_Bar1y);
	t_PreRec->Branch("PreRec_Bar1z",&m_PreRec_Bar1z);
	t_PreRec->Branch("PreRec_Bar1E",&m_PreRec_Bar1E);
	t_PreRec->Branch("PreRec_shower0X", &m_PreRec_shower0X);
	t_PreRec->Branch("PreRec_shower0Y", &m_PreRec_shower0Y);
	t_PreRec->Branch("PreRec_shower0Z", &m_PreRec_shower0Z);
	t_PreRec->Branch("PreRec_shower0E", &m_PreRec_shower0E);
	t_PreRec->Branch("PreRec_shower0T1", &m_PreRec_shower0T1);
	t_PreRec->Branch("PreRec_shower0T2", &m_PreRec_shower0T2);
	t_PreRec->Branch("PreRec_shower1X", &m_PreRec_shower1X);
	t_PreRec->Branch("PreRec_shower1Y", &m_PreRec_shower1Y);
	t_PreRec->Branch("PreRec_shower1Z", &m_PreRec_shower1Z);
	t_PreRec->Branch("PreRec_shower1E", &m_PreRec_shower1E);
	t_PreRec->Branch("PreRec_shower1T1", &m_PreRec_shower1T1);
	t_PreRec->Branch("PreRec_shower1T2", &m_PreRec_shower1T2);
	t_PreRec->Branch("PreRec_NshowerX", &m_PreRec_NshowerX);
	t_PreRec->Branch("PreRec_NshowerY", &m_PreRec_NshowerY);
	t_PreRec->Branch("PreRec_NclusterX", &m_PreRec_NclusterX);
	t_PreRec->Branch("PreRec_NclusterY", &m_PreRec_NclusterY);
	t_PreRec->Branch("PreRec_ClusX_ScndM", &m_ClusX_ScndM);
	t_PreRec->Branch("PreRec_ClusY_ScndM", &m_ClusY_ScndM);
	t_PreRec->Branch("PreRec_chi2", &m_chi2);
	t_PreRec->Branch("PreRec_chi2E", &m_chi2E);
	t_PreRec->Branch("PreRec_chi2Tx", &m_chi2Tx);
	t_PreRec->Branch("PreRec_chi2Ty", &m_chi2Ty);
	t_PreRec->Branch("PreRec_chi2comb_C2", &m_chi2comb);
	t_PreRec->Branch("PreRec_showerX_Nbars", &m_showerX_Nbars);
	t_PreRec->Branch("PreRec_showerX_barx", &m_showerX_barx);
	t_PreRec->Branch("PreRec_showerX_bary", &m_showerX_bary);
	t_PreRec->Branch("PreRec_showerX_barz", &m_showerX_barz);
	t_PreRec->Branch("PreRec_showerX_barE", &m_showerX_barE);
	t_PreRec->Branch("PreRec_showerY_Nbars", &m_showerY_Nbars);
	t_PreRec->Branch("PreRec_showerY_barx", &m_showerY_barx);
	t_PreRec->Branch("PreRec_showerY_bary", &m_showerY_bary);
	t_PreRec->Branch("PreRec_showerY_barz", &m_showerY_barz);
	t_PreRec->Branch("PreRec_showerY_barE", &m_showerY_barE);

	std::cout<<"CRDEcalDigiAlg::m_scale="<<m_scale<<std::endl;
	m_geosvc = service<IGeomSvc>("GeoSvc");
	if ( !m_geosvc )  throw "CRDEcalDigiAlg :Failed to find GeoSvc ...";
	dd4hep::Detector* m_dd4hep = m_geosvc->lcdd();
	if ( !m_dd4hep )  throw "CRDEcalDigiAlg :Failed to get dd4hep::Detector ...";
	m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep);
	//const std::string name_readout = "CaloHitsCollection";
	m_decoder = m_geosvc->getDecoder(_readout);
	if (!m_decoder) {
		error() << "Failed to get the decoder. " << endmsg;
		return StatusCode::FAILURE;
	}
	
	rndm.SetSeed(_seed);
	std::cout<<"CRDEcalDigiAlg::initialize"<<std::endl;
	return GaudiAlgorithm::initialize();
}

StatusCode CRDEcalDigiAlg::execute()
{
	if(_nEvt==0) std::cout<<"CRDEcalDigiAlg::execute Start"<<std::endl;
	std::cout<<"Processing event: "<<_nEvt<<std::endl;
	Clear();

 	const edm4hep::SimCalorimeterHitCollection* SimHitCol =  r_SimCaloCol.get();

	edm4hep::CalorimeterHitCollection* caloVec = w_DigiCaloCol.createAndPut();
	edm4hep::MCRecoCaloAssociationCollection* caloAssoVec = w_CaloAssociationCol.createAndPut();
	edm4hep::SimCalorimeterHitCollection* SimHitCellCol = w_SimCaloTruth.createAndPut();
 	std::vector<edm4hep::SimCalorimeterHit> m_simhitCol; m_simhitCol.clear();


  if(SimHitCol == 0) 
  {
     std::cout<<"not found SimCalorimeterHitCollection"<< std::endl;
     return StatusCode::SUCCESS;
  }
  if(_Debug>=1) std::cout<<"digi, input sim hit size="<< SimHitCol->size() <<std::endl;

	double totE_bar=0;
	double totE_Digi=0;

	//Merge input simhit(steps) to real simhit(bar).
	m_simhitCol = MergeHits(SimHitCol);
	if(_Debug>=1) std::cout<<"Finish Hit Merge, with Nhit: "<<m_simhitCol.size()<<std::endl;



	std::map< unsigned long int, std::vector<CRDEcalDigiEDM::DigiBar> > DigiBlocks; DigiBlocks.clear();

	//Loop in SimHit, digitalize SimHit to DigiBar
	for(int i=0;i<m_simhitCol.size();i++){


		edm4hep::SimCalorimeterHit SimHit = m_simhitCol.at(i);
		if(SimHit.getEnergy()<_Eth) continue;


		CRDEcalDigiEDM::DigiBar hitbar;
		unsigned long long id = SimHit.getCellID();
		hitbar.cellID = id;
		hitbar.system = m_decoder->get(id, "system");
		hitbar.module = m_decoder->get(id, "module");
		hitbar.dlayer = m_decoder->get(id, "dlayer");
		hitbar.block  = m_decoder->get(id, "block");
		hitbar.part   = m_decoder->get(id, "part");
		hitbar.slayer = m_decoder->get(id, "slayer");
		hitbar.bar    = m_decoder->get(id, "bar");


		double Lbar = GetBarLength(hitbar);  //NOTE: Is fixed with geometry CRDEcalBarrel_v01. 
		dd4hep::Position hitpos = m_cellIDConverter->position(id);
		dd4hep::Position barpos(10*hitpos.x(), 10*hitpos.y(), 10*hitpos.z());	//cm to mm.
		hitbar.position = barpos;
		hitbar.T1 = 99999; hitbar.T2 = 99999;
		//if(_Debug>=2) std::cout<<"SimHit contribution size: "<<SimHit.contributions_size()<<std::endl;


		std::vector<CRDEcalDigiEDM::StepDigiOut> DigiLvec; DigiLvec.clear();
		std::vector<CRDEcalDigiEDM::StepDigiOut> DigiRvec; DigiRvec.clear();
		double totQ1 = 0;
		double totQ2 = 0;

		//Loop in all SimHitContribution(G4Step). 
		for(int iCont=0; iCont < SimHit.contributions_size(); ++iCont){
			edm4hep::ConstCaloHitContribution conb = SimHit.getContributions(iCont);
			if( !conb.isAvailable() ) { std::cout<<"CRDEcalDigiAlg  Can not get SimHitContribution: "<<iCont<<std::endl; continue;}

			double en = conb.getEnergy();
			if(en == 0) continue;

			dd4hep::Position steppos(conb.getStepPosition().x, conb.getStepPosition().y, conb.getStepPosition().z);
			dd4hep::Position rpos = steppos-hitbar.position;

			m_step_x.push_back(steppos.x());
			m_step_y.push_back(steppos.y());
			m_step_z.push_back(steppos.z());
			m_step_E.push_back(en);
			m_stepBar_x.push_back(hitbar.position.x());
			m_stepBar_y.push_back(hitbar.position.y());
			m_stepBar_z.push_back(hitbar.position.z());

			if(_Debug>=3){
				cout<<"Cell Pos: "<<hitbar.position.x()<<'\t'<<hitbar.position.y()<<'\t'<<hitbar.position.z()<<endl;
				cout<<"step pos: "<<steppos.x()<<'\t'<<steppos.y()<<'\t'<<steppos.z()<<endl;
				cout<<"Relative pos: "<<rpos.x()<<'\t'<<rpos.y()<<'\t'<<rpos.z()<<endl;
				cout<<"Cell: "<<hitbar.module<<"  "<<hitbar.dlayer<<"  "<<hitbar.slayer<<endl;
			}

			//Get digitalized signal(Q1, Q2, T1, T2) from step
			//Define: 1 is left, 2 is right, clockwise direction in phi. 

			int sign=-999;
			if(hitbar.slayer==1) sign = rpos.z()==0 ? 1 : rpos.z()/fabs(rpos.z());
			else{
				if(hitbar.module==0 || hitbar.module==1 || hitbar.module==7) sign = rpos.x()==0 ?  1: rpos.x()/fabs(rpos.x());
				if(hitbar.module==3 || hitbar.module==4 || hitbar.module==5) sign = rpos.x()==0 ? -1:-rpos.x()/fabs(rpos.x());
				else if(hitbar.module==2) sign = rpos.y()==0 ?  1: rpos.y()/fabs(rpos.y());
				else if(hitbar.module==6) sign = rpos.y()==0 ? -1:-rpos.y()/fabs(rpos.y());
			}
			if(!fabs(sign)) {std::cout<<"ERROR: Wrong bar direction/position!"<<std::endl; continue;}


			double Qi_left = en*exp(-(Lbar/2 + sign*sqrt(rpos.Mag2()))/Latt);	//FIXME: Need to use z rather than magnitude.
			double Qi_right = en*exp(-(Lbar/2 - sign*sqrt(rpos.Mag2()))/Latt);

			if(_Debug>=3){
				cout<<Qi_left<<'\t'<<Qi_right<<endl;
				cout<<Lbar<<'\t'<<sign*sqrt(rpos.Mag2())<<endl;
			}


			double Ti_left = -1; int looptime=0;
			while(Ti_left<0){ 
				Ti_left = Tinit + rndm.Gaus(nMat*(Lbar/2 + sign*sqrt(rpos.Mag2()))/C, Tres); //FIXME: same as Qi.
				looptime++;
				if(looptime>500){ std::cout<<"ERROR: Step "<<iCont<<" can not get a positive left-side time!"<<std::endl; break;}
			}
			if(looptime>500) continue;		
			double Ti_right = -1; looptime=0;
			while(Ti_right<0){ 
				Ti_right = Tinit + rndm.Gaus(nMat*(Lbar/2 - sign*sqrt(rpos.Mag2()))/C, Tres); //FIXME: same as Qi.
				looptime++;
            if(looptime>500){ std::cout<<"ERROR: Step "<<iCont<<" can not get a positive right-side time!"<<std::endl; break;}
			}
			if(looptime>500) continue;		


			m_step_T1.push_back(Ti_left);
			m_step_T2.push_back(Ti_right);

			hitbar.Q1 += Qi_left;
			hitbar.Q2 += Qi_right;
			//hitbar.T1 = hitbar.T1>Ti_left ? Ti_left : hitbar.T1;
			//hitbar.T2 = hitbar.T2>Ti_right ? Ti_right : hitbar.T2;
	
			CRDEcalDigiEDM::StepDigiOut stepoutL, stepoutR;
			stepoutL.T = Ti_left;
			stepoutL.Q = Qi_left;
			stepoutR.T = Ti_right;
			stepoutR.Q = Qi_right;
			DigiLvec.push_back(stepoutL);
			DigiRvec.push_back(stepoutR);
			totQ1 = hitbar.Q1;
			totQ2 = hitbar.Q2;

			//Create truth SimHit(1cm*1cm*1cm cellsize)
			//NOTE: NO cellID for this truth SimHit. 
			dd4hep::Position truthpos = GetCellPos(steppos, hitbar); 
			edm4hep::SimCalorimeterHit simhitTruth = find(SimHitCellCol, truthpos); 
			if(simhitTruth.getCellID()==0){ 
				simhitTruth = SimHitCellCol->create();
				edm4hep::Vector3f m_vec(truthpos.x(), truthpos.y(), truthpos.z());
				simhitTruth.setPosition(m_vec);
				simhitTruth.setCellID(1);
				m_simTruth_slayer.push_back(hitbar.slayer);
				m_simTruth_dlayer.push_back(hitbar.dlayer);
			}

			simhitTruth.addToContributions(conb);
			simhitTruth.setEnergy(simhitTruth.getEnergy()+en );

		}

		//Time digitalization
		//if(_Debug>=2) std::cout<<"Time Digitalize: time at Q >"<<_Qthfrac<<"*totQ"<<std::endl;
		std::sort(DigiLvec.begin(), DigiLvec.end());
		std::sort(DigiRvec.begin(), DigiRvec.end());
		double thQ1=0;
		double thQ2=0;
		for(int iCont=0;iCont<DigiLvec.size();iCont++){
			thQ1 += DigiLvec[iCont].Q;
			if(thQ1>totQ1*_Qthfrac){ 
				hitbar.T1 = DigiLvec[iCont].T; 
				if(_Debug>=3) std::cout<<"Get T1 at index: "<<iCont<<std::endl;
				break;
			}
		}
		for(int iCont=0;iCont<DigiRvec.size();iCont++){
			thQ2 += DigiRvec[iCont].Q;
			if(thQ2>totQ2*_Qthfrac){ 
				hitbar.T2 = DigiRvec[iCont].T; 
				if(_Debug>=3) std::cout<<"Get T2 at index: "<<iCont<<std::endl;
				break;
			}
		}


		totE_bar+=(hitbar.Q1+hitbar.Q2)/2;
		unsigned long int blockID = coder(hitbar);
		DigiBlocks[blockID].push_back(hitbar);
		
		m_simBar_x.push_back(hitbar.position.x());
		m_simBar_y.push_back(hitbar.position.y());
		m_simBar_z.push_back(hitbar.position.z());
		m_simBar_Q1.push_back(hitbar.Q1);
		m_simBar_Q2.push_back(hitbar.Q2);
		m_simBar_T1.push_back(hitbar.T1);
		m_simBar_T2.push_back(hitbar.T2);
		m_simBar_module.push_back(hitbar.module);
		m_simBar_dlayer.push_back(hitbar.dlayer);
		m_simBar_part.push_back(hitbar.part);
		m_simBar_block.push_back(hitbar.block);
		m_simBar_slayer.push_back(hitbar.slayer);

	}
	t_SimCont->Fill();
	t_SimBar->Fill();
	if(_Debug>=1) std::cout<<"End Loop: Bar Digitalization!"<<std::endl;

	
	//DigiHit reconstruction
	if(_Debug>=1) std::cout<<"Block number: "<<DigiBlocks.size()<<std::endl;
	for(auto iter=DigiBlocks.begin(); iter!=DigiBlocks.end();iter++){
		ClearPreRec();

		std::vector<CRDEcalDigiEDM::DigiBar> m_block = iter->second;
		if(m_block.size()==0){ std::cout<<"WARNING: No DigiBar in this block!"<<std::endl; continue;}

		std::vector<CRDEcalDigiEDM::DigiBar> barColX; barColX.clear();
		std::vector<CRDEcalDigiEDM::DigiBar> barColY; barColY.clear();
		for(int i=0;i<m_block.size();i++){
			if(m_block[i].slayer==0) barColX.push_back(m_block[i]);
			else if(m_block[i].slayer==1) barColY.push_back(m_block[i]);
		}

		CRDEcalPreRecAlg preRec; 
		preRec.Eth_SeedWithNeigh  = _Eth_SeedWithNeigh;
		preRec.Eth_SeedWithTot    = _Eth_SeedWithTot;
		preRec.Eth_ShowerWithTot  = _Eth_ShowerWithTot;
		preRec.Eth_ClusterWithTot = _Eth_ClusterWithTot;
		preRec.Sth_split			  	= _Sth_split;
		preRec.Eth_SeedAbs				=	_Eth_SeedAbs;
		preRec.Eth_ShowerAbs			= _Eth_ShowerAbs;
		preRec.Eth_ClusAbs				= _Eth_ClusAbs;
	
		//Step1: Cluster bars in 1D. 
		//std::vector<CRDEcalDigiEDM::BarCollection> barShowerXCol = preRec.Bars2Shower(barColX);
		//std::vector<CRDEcalDigiEDM::BarCollection> barShowerYCol = preRec.Bars2Shower(barColY);

		std::vector<CRDEcalDigiEDM::BarCollection> barShowerXCol; barShowerXCol.clear();
		std::vector<CRDEcalDigiEDM::BarCollection> barShowerYCol; barShowerYCol.clear();
		std::vector<CRDEcalDigiEDM::BarCluster> barClusterX = preRec.Clustering(barColX);
		std::vector<CRDEcalDigiEDM::BarCluster> barClusterY = preRec.Clustering(barColY);

		for(int i=0;i<barClusterX.size();i++){
			double _Etot=0;
			for(int i=0;i<barColX.size();i++) _Etot+=barColX[i].getEnergy();
			if(barClusterX[i].getE()/_Etot < _Eth_ClusterWithTot || barClusterX[i].getE()<_Eth_ClusAbs ) continue;
			std::vector<CRDEcalDigiEDM::BarCollection> showers = preRec.ClusterSplitting(barClusterX[i]);
			if(showers.size()==0) continue;
			barShowerXCol.insert(barShowerXCol.end(), showers.begin(), showers.end());
			m_ClusX_ScndM.push_back(barClusterX[i].getScndMoment());
		}
		for(int i=0;i<barClusterY.size();i++){
			double _Etot=0;
			for(int i=0;i<barColY.size();i++) _Etot+=barColY[i].getEnergy();
			if(barClusterY[i].getE()/_Etot < _Eth_ClusterWithTot || barClusterY[i].getE()<_Eth_ClusAbs) continue;
			std::vector<CRDEcalDigiEDM::BarCollection> showers = preRec.ClusterSplitting(barClusterY[i]);
			if(showers.size()==0) continue;
			barShowerYCol.insert(barShowerYCol.end(), showers.begin(), showers.end());
			m_ClusY_ScndM.push_back(barClusterY[i].getScndMoment());
		}


		//Record 1D clustering output for debugging.
		for(int i=0;i<barColX.size();i++){
			m_PreRec_Bar0x.push_back(barColX[i].position.x());
			m_PreRec_Bar0y.push_back(barColX[i].position.y());
			m_PreRec_Bar0z.push_back(barColX[i].position.z());
			m_PreRec_Bar0E.push_back(barColX[i].getEnergy());
		}
		for(int i=0;i<barColY.size();i++){
			m_PreRec_Bar1x.push_back(barColY[i].position.x());
			m_PreRec_Bar1y.push_back(barColY[i].position.y());
			m_PreRec_Bar1z.push_back(barColY[i].position.z());
			m_PreRec_Bar1E.push_back(barColY[i].getEnergy());
		}
		for(int i=0;i<barShowerXCol.size();i++){ 
			m_PreRec_shower0X.push_back(barShowerXCol[i].getPos().x());
			m_PreRec_shower0Y.push_back(barShowerXCol[i].getPos().y());
			m_PreRec_shower0Z.push_back(barShowerXCol[i].getPos().z());
			m_PreRec_shower0E.push_back(barShowerXCol[i].getE());
			m_PreRec_shower0T2.push_back(barShowerXCol[i].getT2());
			m_PreRec_shower0T1.push_back(barShowerXCol[i].getT1());
		
			m_showerX_Nbars.push_back(barShowerXCol[i].Bars.size());
			for(int ii=0;ii<barShowerXCol[i].Bars.size();ii++){
				m_showerX_barx.push_back(barShowerXCol[i].Bars[ii].position.x());
				m_showerX_bary.push_back(barShowerXCol[i].Bars[ii].position.y());
				m_showerX_barz.push_back(barShowerXCol[i].Bars[ii].position.z());
				m_showerX_barE.push_back(barShowerXCol[i].Bars[ii].getEnergy());
			}
		}
		for(int i=0;i<barShowerYCol.size();i++){ 
			m_PreRec_shower1X.push_back(barShowerYCol[i].getPos().x());
			m_PreRec_shower1Y.push_back(barShowerYCol[i].getPos().y());
			m_PreRec_shower1Z.push_back(barShowerYCol[i].getPos().z());
			m_PreRec_shower1E.push_back(barShowerYCol[i].getE());
			m_PreRec_shower1T1.push_back(barShowerYCol[i].getT1());
			m_PreRec_shower1T2.push_back(barShowerYCol[i].getT2());

			m_showerY_Nbars.push_back(barShowerYCol[i].Bars.size());
			for(int ii=0;ii<barShowerYCol[i].Bars.size();ii++){
				m_showerY_barx.push_back(barShowerYCol[i].Bars[ii].position.x());
				m_showerY_bary.push_back(barShowerYCol[i].Bars[ii].position.y());
				m_showerY_barz.push_back(barShowerYCol[i].Bars[ii].position.z());
				m_showerY_barE.push_back(barShowerYCol[i].Bars[ii].getEnergy());
			}
		}
		m_PreRec_NshowerX = barShowerXCol.size();
		m_PreRec_NshowerY = barShowerYCol.size();
		m_PreRec_NclusterX = barClusterX.size();
		m_PreRec_NclusterY = barClusterY.size();



		//Step2: Reconstruction in different cases
		std::vector<edm4hep::ConstCalorimeterHit> m_caloHits; m_caloHits.clear();

		//Case1(1*N or N*1): Use cross-locating directly. 
		if( barShowerXCol.size()<=1 || barShowerYCol.size()<=1 ){
			if(_Debug>=2) std::cout<<"PreRecAlg: Case1. Shower number X/Y: "<<barShowerXCol.size()<<'\t'<<barShowerYCol.size()<<std::endl;
			//m_caloHits = DigiHitsWithPos(m_block);
			std::vector<CRDEcalDigiEDM::DigiBar> m_ShowerBlock; m_ShowerBlock.clear();
			for(int i=0;i<barShowerXCol.size();i++)  m_ShowerBlock.insert(m_ShowerBlock.end(), barShowerXCol[i].Bars.begin(), barShowerXCol[i].Bars.end());
			for(int i=0;i<barShowerYCol.size();i++)  m_ShowerBlock.insert(m_ShowerBlock.end(), barShowerYCol[i].Bars.begin(), barShowerYCol[i].Bars.end());
			m_caloHits = DigiHitsWithPos(m_ShowerBlock);

		}

		//Case2(N*N): match bars with shower energy
		else if(barShowerXCol.size()==barShowerYCol.size()){
			if(_Debug>=2) std::cout<<"PreRecAlg: Case2. Shower number X/Y: "<<barShowerXCol.size()<<'\t'<<barShowerYCol.size()<<std::endl;
			//m_caloHits = DigiHitsWithEnergy(m_block, barShowerXCol, barShowerYCol);
			m_caloHits = DigiHitsWithMatching(barShowerXCol, barShowerYCol);
		}

		//Case3(M*N)
		else{
			if(_Debug>=2) std::cout<<"PreRecAlg: Case3. Shower number X/Y: "<<barShowerXCol.size()<<'\t'<<barShowerYCol.size()<<std::endl;
			m_caloHits = DigiHitsWithMatchingL2(barShowerXCol, barShowerYCol);
		}

		t_PreRec->Fill();


		//Finish caloHits reconstruction. Stored in m_caloHits. 
		if(_Debug>=2) std::cout<<"After PreRec. CaloHit number: "<<m_caloHits.size()<<std::endl;
		for(int ihit=0;ihit<m_caloHits.size();ihit++){
			totE_Digi+=m_caloHits[ihit].getEnergy();
			m_Rec_x.push_back(m_caloHits[ihit].getPosition()[0]);
			m_Rec_y.push_back(m_caloHits[ihit].getPosition()[1]);
			m_Rec_z.push_back(m_caloHits[ihit].getPosition()[2]);
			m_Rec_E.push_back(m_caloHits[ihit].getEnergy());
			caloVec->push_back(m_caloHits[ihit]);
		}



	}
	t_Rec->Fill();
	if(_Debug>=1) std::cout<<"RecoDigiHit Number: "<<caloVec->size()<<std::endl;


	if(_Debug>=1) std::cout<<"TruthSimHit Number: "<<SimHitCellCol->size()<<std::endl;
	for(int iter=0; iter<SimHitCellCol->size();iter++){
		edm4hep::SimCalorimeterHit m_simhit = SimHitCellCol->at(iter);
		if(!m_simhit.isAvailable()) continue;
		m_simTruth_x.push_back(m_simhit.getPosition()[0]);
		m_simTruth_y.push_back(m_simhit.getPosition()[1]);
		m_simTruth_z.push_back(m_simhit.getPosition()[2]);
		m_simTruth_E.push_back(m_simhit.getEnergy());
	}
	t_SimTruth->Fill();


	std::cout<<"Total Bar Energy: "<<totE_bar<<std::endl;
	std::cout<<"Total DigiHit Energy: "<<totE_Digi<<std::endl;


  _nEvt ++ ;
  return StatusCode::SUCCESS;
}

StatusCode CRDEcalDigiAlg::finalize()
{
	m_wfile->cd();
	t_SimCont->Write();
	t_SimBar->Write();
	t_SimTruth->Write();
	t_Rec->Write();
	t_PreRec->Write();
	m_wfile->Close();

  info() << "Processed " << _nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}

std::vector<edm4hep::SimCalorimeterHit> CRDEcalDigiAlg::MergeHits(const edm4hep::SimCalorimeterHitCollection* m_col){
	std::vector<edm4hep::SimCalorimeterHit> m_mergedhit;
	m_mergedhit.clear();

	for(int iter=0;iter<m_col->size();iter++){
		edm4hep::SimCalorimeterHit m_step = m_col->at(iter);
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



double CRDEcalDigiAlg::GetBarLength(CRDEcalDigiEDM::DigiBar bar){
	//NOTE: Didn't finish reading bar length from geosvc. 
	if(bar.slayer==1) return 460.;
	else return 470.-bar.dlayer*10.;
}

dd4hep::Position CRDEcalDigiAlg::GetCellPos(dd4hep::Position pos, CRDEcalDigiEDM::DigiBar bar){
	dd4hep::Position rpos = pos-bar.position;
	TVector3 vec(0,0,0); 
	if(bar.slayer==1) vec.SetXYZ(0, 0, floor(rpos.z()/10)*10+5 );
	else if(bar.slayer==0){
		if((bar.module==0||bar.module==4) && bar.dlayer%2==1) vec.SetXYZ(floor(rpos.x()/10)*10+5,0,0);
		if((bar.module==0||bar.module==4) && bar.dlayer%2==0) vec.SetXYZ(floor((rpos.x()-5)/10)*10+10,0,0);
		if((bar.module==2||bar.module==6) && bar.dlayer%2==1) vec.SetXYZ(0, floor(rpos.y()/10)*10+5,0);
		if((bar.module==2||bar.module==6) && bar.dlayer%2==0) vec.SetXYZ(0, floor((rpos.y()-5)/10)*10+10,0);
		if(bar.module==1 || bar.module==5){
			TVector3 unitv(1./sqrt(2), -1./sqrt(2), 0);
			if(bar.dlayer%2==1) vec = (floor(rpos.Dot(unitv)/10)*10+5)*unitv;
			if(bar.dlayer%2==0) vec = (floor((rpos.Dot(unitv)-5)/10)*10+10)*unitv;
		}
		if(bar.module==3 || bar.module==7){
			TVector3 unitv(1./sqrt(2), 1./sqrt(2), 0);
			if(bar.dlayer%2==1) vec = (floor(rpos.Dot(unitv)/10)*10+5)*unitv;
			if(bar.dlayer%2==0) vec = (floor((rpos.Dot(unitv)-5)/10)*10+10)*unitv;
		}
	}
	dd4hep::Position relv(vec.x(), vec.y(), vec.z());
	return relv+bar.position;
}

edm4hep::SimCalorimeterHit CRDEcalDigiAlg::find(edm4hep::SimCalorimeterHitCollection* m_col, dd4hep::Position pos){
   for(int i=0;i<m_col->size();i++){
      edm4hep::SimCalorimeterHit hit = m_col->at(i);
		dd4hep::Position ipos(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z);
		if(ipos==pos) return hit;
	}
   edm4hep::SimCalorimeterHit hit;
   hit.setCellID(0);
   return hit;
}

edm4hep::SimCalorimeterHit CRDEcalDigiAlg::find(std::vector<edm4hep::SimCalorimeterHit> m_col, unsigned long long cellid){
   for(int i=0;i<m_col.size();i++){
		edm4hep::SimCalorimeterHit hit=m_col.at(i);
		if(hit.getCellID() == cellid) return hit;
	}
	edm4hep::SimCalorimeterHit hit;
	hit.setCellID(0);
	return hit;
}

unsigned long int CRDEcalDigiAlg::coder(CRDEcalDigiEDM::DigiBar bar){
	return 10099*bar.module + 6911*bar.dlayer + 727*bar.part + 31*bar.block;
}


void CRDEcalDigiAlg::Clear(){
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
	m_simBar_dlayer.clear();
	m_simBar_part.clear();
	m_simBar_block.clear();
	m_simBar_slayer.clear();
	m_simTruth_x.clear();
	m_simTruth_y.clear();
	m_simTruth_z.clear();
	m_simTruth_E.clear();
	m_simTruth_dlayer.clear();
	m_simTruth_slayer.clear();
	m_Rec_x.clear();
	m_Rec_y.clear();
	m_Rec_z.clear();
	m_Rec_E.clear();
}

void CRDEcalDigiAlg::ClearPreRec(){
	m_PreRec_Bar0x.clear();
	m_PreRec_Bar0y.clear();
	m_PreRec_Bar0z.clear();
	m_PreRec_Bar0E.clear();
	m_PreRec_Bar1x.clear();
	m_PreRec_Bar1y.clear();
	m_PreRec_Bar1z.clear();
	m_PreRec_Bar1E.clear();
	m_PreRec_shower0X.clear();
	m_PreRec_shower0Y.clear();
	m_PreRec_shower0Z.clear();
	m_PreRec_shower0E.clear();
	m_PreRec_shower0T1.clear();
	m_PreRec_shower0T2.clear();
	m_PreRec_shower1X.clear();
	m_PreRec_shower1Y.clear();
	m_PreRec_shower1Z.clear();
	m_PreRec_shower1E.clear();
	m_PreRec_shower1T1.clear();
	m_PreRec_shower1T2.clear();
	m_PreRec_NshowerX=-999;
	m_PreRec_NshowerY=-999;
	m_PreRec_NclusterX=-999;
	m_PreRec_NclusterY=-999;
	m_ClusX_ScndM.clear();
	m_ClusY_ScndM.clear();
	m_chi2.clear();
	m_chi2E.clear();
	m_chi2Tx.clear();
	m_chi2Ty.clear();
	m_chi2comb.clear();

	m_showerX_Nbars.clear();
	m_showerX_barx.clear();
	m_showerX_bary.clear();
	m_showerX_barz.clear();
	m_showerX_barE.clear();
	m_showerY_Nbars.clear();
	m_showerY_barx.clear();
	m_showerY_bary.clear();
	m_showerY_barz.clear();
	m_showerY_barE.clear();
}
