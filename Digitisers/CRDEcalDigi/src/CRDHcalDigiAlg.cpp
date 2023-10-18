// /* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
// // Unit in code: mm, ns. 
// // NOTE: This digitialization highly matches detector geometry CRDEcalBarrel_v01.
// // TODO: read geometry info automatically.  

#include "CRDHcalDigiAlg.h" 

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

DECLARE_COMPONENT( CRDHcalDigiAlg )

CRDHcalDigiAlg::CRDHcalDigiAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{
  
	// Input collections
	declareProperty("SimCaloHitCollection", r_SimCaloCol, "Handle of the Input SimCaloHit collection");
  
	// Output collections
	declareProperty("CaloHitCollection", w_DigiCaloCol, "Handle of Digi CaloHit collection");
	declareProperty("CaloAssociationCollection", w_CaloAssociationCol, "Handle of CaloAssociation collection");
   
}

StatusCode CRDHcalDigiAlg::initialize()
{

	std::string s_outfile = _filename;
	m_wfile = new TFile(s_outfile.c_str(), "recreate");
	t_SimCont = new TTree("SimStep", "SimStep");
	t_simHit = new TTree("simHit", "simHit");
	t_SimCont->Branch("step_x", &m_step_x);
	t_SimCont->Branch("step_y", &m_step_y);
	t_SimCont->Branch("step_z", &m_step_z);
	t_SimCont->Branch("step_E", &m_step_E);

	t_simHit->Branch("simHit_x", &m_simHit_x);
	t_simHit->Branch("simHit_y", &m_simHit_y);
	t_simHit->Branch("simHit_z", &m_simHit_z);
	t_simHit->Branch("simHit_E", &m_simHit_E);
	t_simHit->Branch("simHit_steps", &m_simHit_steps);

	t_simHit->Branch("simHit_module", &m_simHit_module);
	t_simHit->Branch("simHit_stave", &m_simHit_stave);
	t_simHit->Branch("simHit_layer", &m_simHit_layer);
	t_simHit->Branch("simHit_tower", &m_simHit_tower);
	t_simHit->Branch("simHit_slice", &m_simHit_slice);
	t_simHit->Branch("simHit_cellID", &m_simHit_cellID);

	std::cout<<"CRDHcalDigiAlg::m_scale="<<m_scale<<std::endl;
	m_geosvc = service<IGeomSvc>("GeomSvc");
	if ( !m_geosvc )  throw "CRDHcalDigiAlg :Failed to find GeomSvc ...";
	dd4hep::Detector* m_dd4hep = m_geosvc->lcdd();
	if ( !m_dd4hep )  throw "CRDHcalDigiAlg :Failed to get dd4hep::Detector ...";
	m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep);
   	m_decoder = m_geosvc->getDecoder(_readoutName);
	if (!m_decoder) {
		error() << "Failed to get the decoder. " << endmsg;
		return StatusCode::FAILURE;
	}

	//m_edmsvc = service<ICRDEcalSvc>("CRDEcalSvc");
	//if ( !m_edmsvc )  throw "CRDHcalDigiAlg :Failed to find CRDEcalSvc ...";

	rndm.SetSeed(_seed);
	std::cout<<"CRDHcalDigiAlg::initialize"<<std::endl;
	return GaudiAlgorithm::initialize();
}

StatusCode CRDHcalDigiAlg::execute()
{
	if(_nEvt==0) std::cout<<"CRDHcalDigiAlg::execute Start"<<std::endl;
	std::cout<<"Processing event: "<<_nEvt<<std::endl;
   	if(_nEvt<_Nskip){ _nEvt++; return StatusCode::SUCCESS; }

	Clear();

 	const edm4hep::SimCalorimeterHitCollection* SimHitCol =  r_SimCaloCol.get();

	edm4hep::CalorimeterHitCollection* caloVec = w_DigiCaloCol.createAndPut();
	edm4hep::MCRecoCaloAssociationCollection* caloAssoVec = w_CaloAssociationCol.createAndPut();
 	std::vector<edm4hep::SimCalorimeterHit> m_simhitCol; m_simhitCol.clear();
  	// std::vector<CaloHit> m_barCol; m_barCol.clear(); 

	if(SimHitCol == 0) 
	{
		std::cout<<"not found SimCalorimeterHitCollection"<< std::endl;
		return StatusCode::SUCCESS;
	}
  	if(_Debug>=1) std::cout<<"digi, input sim hit size="<< SimHitCol->size() <<std::endl;

	double totE_hit=0;

	//Merge input simhit(steps) to real simhit(bar).
	MergeHits(*SimHitCol, m_simhitCol);
	if(_Debug>=1) std::cout<<"Finish Hit Merge, with Nhit: "<<m_simhitCol.size()<<std::endl;

	//Loop in SimHit, digitalize SimHit to DigiBar
	for(int i=0;i<m_simhitCol.size();i++)
	{

		auto SimHit = m_simhitCol.at(i);
		if(SimHit.getEnergy()<_Eth) continue;

		unsigned long long id = SimHit.getCellID();
		CaloHit hcalhit;
		hcalhit.setcellID( id);
		hcalhit.setcellID(	m_decoder->get(id, "system"), 
											m_decoder->get(id, "module"), 
											m_decoder->get(id, "layer"), 
											m_decoder->get(id, "stave"), 
											m_decoder->get(id, "tower"), 
											m_decoder->get(id, "slice"));

		dd4hep::Position hitpos = m_cellIDConverter->position(id);
    	TVector3 calohitpos(10*hitpos.x(), 10*hitpos.y(), 10*hitpos.z()); //cm to mm.
		hcalhit.setPosition(calohitpos);
		//hcalhit.T1 = 99999; hcalhit.T2 = 99999;
		//if(_Debug>=2) std::cout<<"SimHit contribution size: "<<SimHit.contributions_size()<<std::endl;

		double totQ = 0;

		//Loop in all SimHitContribution(G4Step). 
		for(int iCont=0; iCont < SimHit.contributions_size(); ++iCont){
			auto conb = SimHit.getContributions(iCont);
			if( !conb.isAvailable() ) { std::cout<<"CRDHcalDigiAlg  Can not get SimHitContribution: "<<iCont<<std::endl; continue;}

			double en = conb.getEnergy();
			if(en == 0) continue;

			TVector3 steppos(conb.getStepPosition().x, conb.getStepPosition().y, conb.getStepPosition().z);
			TVector3 rpos = steppos-hcalhit.getPosition();

			m_step_x.push_back(steppos.x());
			m_step_y.push_back(steppos.y());
			m_step_z.push_back(steppos.z());
			m_step_E.push_back(en);

			// if(_Debug>=3){
			// 	cout<<"Cell Pos: "<<hcalhit.getPosition().x()<<'\t'<<hcalhit.getPosition().y()<<'\t'<<hcalhit.getPosition().z()<<endl;
			// 	cout<<"step pos: "<<steppos.x()<<'\t'<<steppos.y()<<'\t'<<steppos.z()<<endl;
			// 	cout<<"Relative pos: "<<rpos.x()<<'\t'<<rpos.y()<<'\t'<<rpos.z()<<endl;
			// 	cout<<"Cell: "<<hcalhit.getModule()<<"  "<<hcalhit.getDlayer()<<"  "<<hcalhit.getSlayer()<<endl;
			// }

			totQ += en;
		}

		hcalhit.setE(totQ);
		hcalhit.setStep_numbers(SimHit.contributions_size());

    	//End bar digitization. 

		//2 hits with double-readout time. 
		edm4hep::Vector3f m_pos(hcalhit.getPosition().X(), hcalhit.getPosition().Y(), hcalhit.getPosition().Z());
		auto digiHit1 = caloVec->create();
		digiHit1.setCellID(hcalhit.getcellID());
		digiHit1.setEnergy(hcalhit.getEnergy());
		digiHit1.setPosition(m_pos);

		auto rel = caloAssoVec->create();
		rel.setRec(digiHit1);
		rel.setSim(SimHit);
		rel.setWeight(1.);
			
		//Temp: write into trees. 
		m_simHit_x.push_back(hcalhit.getPosition().x());
		m_simHit_y.push_back(hcalhit.getPosition().y());
		m_simHit_z.push_back(hcalhit.getPosition().z());
		m_simHit_E.push_back(hcalhit.getEnergy());
		m_simHit_steps.push_back(hcalhit.getStep_numbers());
		m_simHit_module.push_back(hcalhit.getModule());
		m_simHit_stave.push_back(hcalhit.getStave());
		m_simHit_layer.push_back(hcalhit.getLayer());
		m_simHit_slice.push_back(hcalhit.getSlice());
		m_simHit_tower.push_back(hcalhit.getTower());
		m_simHit_cellID.push_back(hcalhit.getcellID());
	}

	t_SimCont->Fill();
	t_simHit->Fill();
	// if(_Debug>=1) std::cout<<"End Loop: Bar Digitalization!"<<std::endl;

	_nEvt ++ ;
	//delete SimHitCol, caloVec, caloAssoVec; 
	m_simhitCol.clear();
	return StatusCode::SUCCESS;
}

StatusCode CRDHcalDigiAlg::finalize()
{
	m_wfile->cd();
	t_SimCont->Write();
	t_simHit->Write();
	m_wfile->Close();

	info() << "Processed " << _nEvt << " events " << endmsg;
	delete m_wfile, t_SimCont, t_simHit; 
	delete m_cellIDConverter, m_decoder, m_geosvc;
	return GaudiAlgorithm::finalize();
}

StatusCode CRDHcalDigiAlg::MergeHits( const edm4hep::SimCalorimeterHitCollection& m_col, std::vector<edm4hep::SimCalorimeterHit>& m_hits ){

  	m_hits.clear(); 
	std::vector<edm4hep::MutableSimCalorimeterHit> m_mergedhit;
	m_mergedhit.clear();

	for(int iter=0; iter<m_col.size(); iter++){
		edm4hep::SimCalorimeterHit m_step = m_col[iter];
		if(!m_step.isAvailable()){ cout<<"ERROR HIT!"<<endl; continue;}
		if(m_step.getEnergy()==0) continue;
		unsigned long long cellid = m_step.getCellID();
		dd4hep::Position hitpos = m_cellIDConverter->position(cellid);
		edm4hep::Vector3f pos(hitpos.x()*10, hitpos.y()*10, hitpos.z()*10);

		edm4hep::MutableCaloHitContribution conb;
		conb.setEnergy(m_step.getEnergy());
		conb.setStepPosition(m_step.getPosition());
		//---oooOOO000OOOooo---
		// if(m_step.contributions_size()==1)
		// 	conb.setTime(m_step.getContributions(0).getTime());
		// else{
		// 	cout<<"yyy: MergeHits(), m_step.contributions_size()!=1"<<endl;
		// }
		// //---oooOOO000OOOooo---
		edm4hep::MutableSimCalorimeterHit m_hit = find(m_mergedhit, cellid);
		if(m_hit.getCellID()==0){
			//m_hit = new edm4hep::SimCalorimeterHit();
			m_hit.setCellID(cellid);
			m_hit.setPosition(pos);
			m_mergedhit.push_back(m_hit);
		}
		m_hit.addToContributions(conb);
		m_hit.setEnergy(m_hit.getEnergy()+m_step.getEnergy());
	}

	for(auto iter = m_mergedhit.begin(); iter!=m_mergedhit.end(); iter++){
		edm4hep::SimCalorimeterHit constsimhit = *iter; 
		m_hits.push_back( constsimhit );  
	}
	return StatusCode::SUCCESS; 
}

edm4hep::MutableSimCalorimeterHit CRDHcalDigiAlg::find(const std::vector<edm4hep::MutableSimCalorimeterHit>& m_col, unsigned long long& cellid) const{
   for(int i=0;i<m_col.size();i++){
		edm4hep::MutableSimCalorimeterHit hit=m_col.at(i);
		if(hit.getCellID() == cellid) return hit;
	}
	edm4hep::MutableSimCalorimeterHit hit ;
	hit.setCellID(0);
	return hit;
}

void CRDHcalDigiAlg::Clear(){
	m_step_x.clear();
	m_step_y.clear();
	m_step_z.clear();
	m_step_E.clear();
	m_simHit_x.clear();
	m_simHit_y.clear();
	m_simHit_z.clear();
	m_simHit_E.clear();
	m_simHit_steps.clear();
	m_simHit_module.clear();
	m_simHit_stave.clear();
	m_simHit_layer.clear();
	m_simHit_slice.clear();
	m_simHit_tower.clear();
  	m_simHit_cellID.clear();
}

