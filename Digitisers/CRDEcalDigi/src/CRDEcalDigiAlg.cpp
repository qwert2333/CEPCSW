/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
// Unit in code: mm, ns. 
// NOTE: This digitialization highly matches detector geometry CRDEcalBarrel_v01. 
#include "CRDEcalDigiAlg.h"

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
	t_SimBar->Branch("simBar_stave", &m_simBar_stave);
	t_SimBar->Branch("simBar_dlayer", &m_simBar_dlayer);
	t_SimBar->Branch("simBar_part", &m_simBar_part);
	t_SimBar->Branch("simBar_slayer", &m_simBar_slayer);
	t_SimTruth->Branch("simTruth_x", &m_simTruth_x);
	t_SimTruth->Branch("simTruth_y", &m_simTruth_y);
	t_SimTruth->Branch("simTruth_z", &m_simTruth_z);
	t_SimTruth->Branch("simTruth_E", &m_simTruth_E);
	t_SimTruth->Branch("simTruth_dlayer", &m_simTruth_dlayer);
	t_SimTruth->Branch("simTruth_slayer", &m_simTruth_slayer);

	std::cout<<"CRDEcalDigiAlg::m_scale="<<m_scale<<std::endl;
	m_geosvc = service<IGeomSvc>("GeoSvc");
	if ( !m_geosvc )  throw "CRDEcalDigiAlg :Failed to find GeoSvc ...";
	dd4hep::Detector* m_dd4hep = m_geosvc->lcdd();
	if ( !m_dd4hep )  throw "CRDEcalDigiAlg :Failed to get dd4hep::Detector ...";
	m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep);
   m_decoder = m_geosvc->getDecoder(_readout);
	if (!m_decoder) {
		error() << "Failed to get the decoder. " << endmsg;
		return StatusCode::FAILURE;
	}
	
	m_edmsvc = service<ICRDEcalEDMSvc>("CRDEcalEDMSvc");
	if ( !m_edmsvc )  throw "CRDEcalDigiAlg :Failed to find CRDEcalEDMSvc ...";

	rndm.SetSeed(_seed);
	std::cout<<"CRDEcalDigiAlg::initialize"<<std::endl;
	return GaudiAlgorithm::initialize();
}

StatusCode CRDEcalDigiAlg::execute()
{
	if(_nEvt==0) std::cout<<"CRDEcalDigiAlg::execute Start"<<std::endl;
	std::cout<<"Processing event: "<<_nEvt<<std::endl;
   if(_nEvt<_Nskip){ _nEvt++; return StatusCode::SUCCESS; }

  m_edmsvc->ClearSystem(); 
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
	m_simhitCol = MergeHits(*SimHitCol);
	if(_Debug>=1) std::cout<<"Finish Hit Merge, with Nhit: "<<m_simhitCol.size()<<std::endl;

	std::map< unsigned long int, CRDEcalEDM::DigiBlock > DigiBlocks; DigiBlocks.clear();

	//Loop in SimHit, digitalize SimHit to DigiBar
	for(int i=0;i<m_simhitCol.size();i++){

		edm4hep::SimCalorimeterHit SimHit = m_simhitCol.at(i);
		if(SimHit.getEnergy()<_Eth) continue;


		unsigned long long id = SimHit.getCellID();
		CRDEcalEDM::CRDCaloBar hitbar;
		hitbar.setcellID( id);
		hitbar.setcellID(	m_decoder->get(id, "system"), 
											m_decoder->get(id, "module"), 
											m_decoder->get(id, "stave"), 
											m_decoder->get(id, "dlayer"), 
											m_decoder->get(id, "part"), 
											m_decoder->get(id, "slayer"),
											m_decoder->get(id, "bar"));

		double Lbar = GetBarLength(hitbar);  //NOTE: Is fixed with geometry CRDEcalBarrel_v01. 
		dd4hep::Position hitpos = m_cellIDConverter->position(id);
		dd4hep::Position barpos(10*hitpos.x(), 10*hitpos.y(), 10*hitpos.z());	//cm to mm.
		hitbar.setPosition(barpos);
		//hitbar.T1 = 99999; hitbar.T2 = 99999;
		//if(_Debug>=2) std::cout<<"SimHit contribution size: "<<SimHit.contributions_size()<<std::endl;


		std::vector<CRDEcalEDM::CRDCaloStep> DigiLvec; DigiLvec.clear();
		std::vector<CRDEcalEDM::CRDCaloStep> DigiRvec; DigiRvec.clear();
		double totQ1 = 0;
		double totQ2 = 0;

		//Loop in all SimHitContribution(G4Step). 
		for(int iCont=0; iCont < SimHit.contributions_size(); ++iCont){
			edm4hep::ConstCaloHitContribution conb = SimHit.getContributions(iCont);
			if( !conb.isAvailable() ) { std::cout<<"CRDEcalDigiAlg  Can not get SimHitContribution: "<<iCont<<std::endl; continue;}

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
				cout<<"Cell: "<<hitbar.getModule()<<"  "<<hitbar.getDlayer()<<"  "<<hitbar.getSlayer()<<endl;
			}

			//Get digitalized signal(Q1, Q2, T1, T2) from step
			//Define: 1 is left, 2 is right, clockwise direction in phi. 

			int sign=-999;
			if(hitbar.getSlayer()==1) sign = rpos.z()==0 ? 1 : rpos.z()/fabs(rpos.z());
			else{
				if(hitbar.getModule()==0 || hitbar.getModule()==1 || hitbar.getModule()==7) sign = rpos.x()==0 ?  1: rpos.x()/fabs(rpos.x());
				if(hitbar.getModule()==3 || hitbar.getModule()==4 || hitbar.getModule()==5) sign = rpos.x()==0 ? -1:-rpos.x()/fabs(rpos.x());
				else if(hitbar.getModule()==2) sign = rpos.y()==0 ?  1: rpos.y()/fabs(rpos.y());
				else if(hitbar.getModule()==6) sign = rpos.y()==0 ? -1:-rpos.y()/fabs(rpos.y());
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
			totQ1 += Qi_left;
			totQ2 += Qi_right;
	
			CRDEcalEDM::CRDCaloStep stepoutL, stepoutR;
			stepoutL.setQ(Qi_left); stepoutL.setT(Ti_left);
			stepoutR.setQ(Qi_right); stepoutR.setT(Ti_right);
			DigiLvec.push_back(stepoutL);
			DigiRvec.push_back(stepoutR);

			//Create truth SimHit(1cm*1cm*1cm cellsize)
			//NOTE: NO cellID for this truth SimHit. 
			dd4hep::Position truthpos = GetCellPos(steppos, hitbar); 
			edm4hep::SimCalorimeterHit simhitTruth = find(*SimHitCellCol, truthpos); 
			if(simhitTruth.getCellID()==0){ 
				simhitTruth = SimHitCellCol->create();
				edm4hep::Vector3f m_vec(truthpos.x(), truthpos.y(), truthpos.z());
				simhitTruth.setPosition(m_vec);
				simhitTruth.setCellID(1);
				m_simTruth_slayer.push_back(hitbar.getSlayer());
				m_simTruth_dlayer.push_back(hitbar.getDlayer());
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
		hitbar.setQ(totQ1, totQ2);
		hitbar.setT(thT1, thT2);


		totE_bar+=(hitbar.getQ1()+hitbar.getQ2())/2;
		unsigned long int blockID = coder(hitbar);
		DigiBlocks[blockID].push_back(hitbar);
		
		m_simBar_x.push_back(hitbar.getPosition().x());
		m_simBar_y.push_back(hitbar.getPosition().y());
		m_simBar_z.push_back(hitbar.getPosition().z());
		m_simBar_Q1.push_back(hitbar.getQ1());
		m_simBar_Q2.push_back(hitbar.getQ2());
		m_simBar_T1.push_back(hitbar.getT1());
		m_simBar_T2.push_back(hitbar.getT2());
		m_simBar_module.push_back(hitbar.getModule());
		m_simBar_stave.push_back(hitbar.getStave());
		m_simBar_dlayer.push_back(hitbar.getDlayer());
		m_simBar_part.push_back(hitbar.getPart());
		m_simBar_slayer.push_back(hitbar.getSlayer());

	}
	t_SimCont->Fill();
	t_SimBar->Fill();
	if(_Debug>=1) std::cout<<"End Loop: Bar Digitalization!"<<std::endl;


	if(_Debug>=1) std::cout<<"TruthSimHit Number: "<<SimHitCellCol->size()<<std::endl;
	for(int iter=0; iter<SimHitCellCol->size();iter++){
		edm4hep::SimCalorimeterHit m_simhit = SimHitCellCol->at(iter);
		if(!m_simhit.isAvailable()) continue;
		m_simTruth_x.push_back(m_simhit.getPosition()[0]);
		m_simTruth_y.push_back(m_simhit.getPosition()[1]);
		m_simTruth_z.push_back(m_simhit.getPosition()[2]);
		m_simTruth_E.push_back(m_simhit.getEnergy());

      edm4hep::CalorimeterHit hit;
      hit.setCellID(0);
      hit.setPosition(m_simhit.getPosition());
      hit.setEnergy(m_simhit.getEnergy());
      caloVec->push_back(hit);
	}
	t_SimTruth->Fill();
	std::cout<<"Total Bar Energy: "<<totE_bar<<std::endl;


	std::vector<CRDEcalEDM::CRDCaloBlock> blockVec; blockVec.clear();
	for(auto iter=DigiBlocks.begin(); iter!=DigiBlocks.end();iter++){ 
      CRDEcalEDM::DigiBlock m_block = iter->second;
      if(m_block.size()==0) continue; 
      CRDEcalEDM::CRDCaloBlock m_calobl; m_calobl.Clear();
      m_calobl.setForced(false);
      std::vector<CRDEcalEDM::CRDCaloBar> barXCol, barYCol; barXCol.clear(); barYCol.clear(); 
      for(int i=0;i<m_block.size();i++){
         if(m_block[i].getSlayer()==0) barXCol.push_back(m_block[i]);
         else barYCol.push_back(m_block[i]);
      }
      m_calobl.setBarCol( barXCol, barYCol );
      m_calobl.setIDInfo( m_block[0].getModule(), m_block[0].getStave(), m_block[0].getDlayer(), m_block[0].getPart());
      blockVec.push_back(m_calobl);
  }
	m_edmsvc->setDigiSystem( blockVec ); 


  _nEvt ++ ;
  return StatusCode::SUCCESS;
}

StatusCode CRDEcalDigiAlg::finalize()
{
	m_wfile->cd();
	t_SimCont->Write();
	t_SimBar->Write();
	t_SimTruth->Write();
	m_wfile->Close();

  info() << "Processed " << _nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}

std::vector<edm4hep::SimCalorimeterHit> CRDEcalDigiAlg::MergeHits(const edm4hep::SimCalorimeterHitCollection& m_col){
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



double CRDEcalDigiAlg::GetBarLength(CRDEcalEDM::CRDCaloBar& bar){
	//TODO: reading bar length from geosvc. 
	if(bar.getSlayer()==1) return 418.;
	else return 470.-bar.getDlayer()*10.;
}

dd4hep::Position CRDEcalDigiAlg::GetCellPos(dd4hep::Position& pos, CRDEcalEDM::CRDCaloBar& bar){
	dd4hep::Position rpos = pos-bar.getPosition();
	TVector3 vec(0,0,0); 
	if(bar.getSlayer()==1) vec.SetXYZ(0, 0, floor(rpos.z()/10)*10+5 );
	else if(bar.getSlayer()==0){
		if((bar.getModule()==0||bar.getModule()==4) && bar.getDlayer()%2==1) vec.SetXYZ(floor(rpos.x()/10)*10+5,0,0);
		if((bar.getModule()==0||bar.getModule()==4) && bar.getDlayer()%2==0) vec.SetXYZ(floor((rpos.x()-5)/10)*10+10,0,0);
		if((bar.getModule()==2||bar.getModule()==6) && bar.getDlayer()%2==1) vec.SetXYZ(0, floor(rpos.y()/10)*10+5,0);
		if((bar.getModule()==2||bar.getModule()==6) && bar.getDlayer()%2==0) vec.SetXYZ(0, floor((rpos.y()-5)/10)*10+10,0);
		if(bar.getModule()==1 || bar.getModule()==5){
			TVector3 unitv(1./sqrt(2), -1./sqrt(2), 0);
			if(bar.getDlayer()%2==1) vec = (floor(rpos.Dot(unitv)/10)*10+5)*unitv;
			if(bar.getDlayer()%2==0) vec = (floor((rpos.Dot(unitv)-5)/10)*10+10)*unitv;
		}
		if(bar.getModule()==3 || bar.getModule()==7){
			TVector3 unitv(1./sqrt(2), 1./sqrt(2), 0);
			if(bar.getDlayer()%2==1) vec = (floor(rpos.Dot(unitv)/10)*10+5)*unitv;
			if(bar.getDlayer()%2==0) vec = (floor((rpos.Dot(unitv)-5)/10)*10+10)*unitv;
		}
	}
	dd4hep::Position relv(vec.x(), vec.y(), vec.z());
	return relv+bar.getPosition();
}

edm4hep::SimCalorimeterHit CRDEcalDigiAlg::find(edm4hep::SimCalorimeterHitCollection& m_col, dd4hep::Position& pos){
   for(int i=0;i<m_col.size();i++){
      edm4hep::SimCalorimeterHit hit = m_col[i];
		dd4hep::Position ipos(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z);
		if(ipos==pos) return hit;
	}
   edm4hep::SimCalorimeterHit hit;
   hit.setCellID(0);
   return hit;
}

edm4hep::SimCalorimeterHit CRDEcalDigiAlg::find(std::vector<edm4hep::SimCalorimeterHit>& m_col, unsigned long long& cellid){
   for(int i=0;i<m_col.size();i++){
		edm4hep::SimCalorimeterHit hit=m_col.at(i);
		if(hit.getCellID() == cellid) return hit;
	}
	edm4hep::SimCalorimeterHit hit;
	hit.setCellID(0);
	return hit;
}

unsigned long int CRDEcalDigiAlg::coder(CRDEcalEDM::CRDCaloBar& bar){
   return 10099*bar.getModule() + 6911*bar.getDlayer() + 727*bar.getPart() + 31*bar.getStave();
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
	m_simBar_stave.clear();
	m_simBar_dlayer.clear();
	m_simBar_part.clear();
	m_simBar_slayer.clear();
	m_simTruth_x.clear();
	m_simTruth_y.clear();
	m_simTruth_z.clear();
	m_simTruth_E.clear();
	m_simTruth_dlayer.clear();
	m_simTruth_slayer.clear();
}

