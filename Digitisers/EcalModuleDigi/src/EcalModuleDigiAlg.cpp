/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
// Unit in code: mm, ns. 
// NOTE: This digitialization highly matches detector geometry CRDEcalBarrel_v01. 
#include "EcalModuleDigiAlg.h"

#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Vector3f.h"

#include "DD4hep/Detector.h"
#include <DD4hep/Objects.h>
#include <DDRec/CellIDPositionConverter.h>

#include "TRandom3.h"
#include "TVector3.h"
#include <math.h>
#include <cmath>
#include <algorithm>
#include <map>

#define C 299.79  // unit: mm/ns
#define PI 3.141592653
using namespace std;
using namespace dd4hep;
DECLARE_COMPONENT( EcalModuleDigiAlg )

EcalModuleDigiAlg::EcalModuleDigiAlg(const std::string& name, ISvcLocator* svcLoc)
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

StatusCode EcalModuleDigiAlg::initialize()
{

	std::string s_outfile = _filename;
	m_wfile = new TFile( s_outfile.c_str(), "recreate");
	t_SimBar = new TTree("SimBarHit", "SimBarHit");
	t_SimTruth = new TTree("TruthSimHit", "TruthSimHit");
	t_Rec = new TTree("RecHit", "RecHit");
	t_SimBar->Branch("simBar_x", &m_simBar_x);
	t_SimBar->Branch("simBar_y", &m_simBar_y);
	t_SimBar->Branch("simBar_z", &m_simBar_z);
	t_SimBar->Branch("simBar_T1", &m_simBar_T1);
	t_SimBar->Branch("simBar_T2", &m_simBar_T2);
	t_SimBar->Branch("simBar_Q1", &m_simBar_Q1);
	t_SimBar->Branch("simBar_Q2", &m_simBar_Q2);
	t_SimBar->Branch("simBar_group", &m_simBar_group);
	t_SimBar->Branch("simBar_block", &m_simBar_block);
	t_SimBar->Branch("simBar_slayer", &m_simBar_slayer);
	t_SimTruth->Branch("simTruth_x", &m_simTruth_x);
	t_SimTruth->Branch("simTruth_y", &m_simTruth_y);
	t_SimTruth->Branch("simTruth_z", &m_simTruth_z);
	t_SimTruth->Branch("simTruth_E", &m_simTruth_E);
	t_SimTruth->Branch("simTruth_group", &m_simTruth_group);
	t_SimTruth->Branch("simTruth_block", &m_simTruth_block);
	t_SimTruth->Branch("simTruth_slayer", &m_simTruth_slayer);
	t_Rec->Branch("Rec_x", &m_Rec_x);
	t_Rec->Branch("Rec_y", &m_Rec_y);
	t_Rec->Branch("Rec_z", &m_Rec_z);
	t_Rec->Branch("Rec_E", &m_Rec_E);

	std::cout<<"EcalModuleDigiAlg::m_scale="<<m_scale<<std::endl;
	m_geosvc = service<IGeomSvc>("GeoSvc");
	if ( !m_geosvc )  throw "EcalModuleDigiAlg :Failed to find GeoSvc ...";
	dd4hep::Detector* m_dd4hep = m_geosvc->lcdd();
	if ( !m_dd4hep )  throw "EcalModuleDigiAlg :Failed to get dd4hep::Detector ...";
	m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep);
	const std::string name_readout = "CaloHitsCollection";
	m_decoder = m_geosvc->getDecoder(name_readout);
	if (!m_decoder) {
		error() << "Failed to get the decoder. " << endmsg;
		return StatusCode::FAILURE;
	}

	std::cout<<"EcalModuleDigiAlg::initialize"<<std::endl;
	return GaudiAlgorithm::initialize();
}

StatusCode EcalModuleDigiAlg::execute()
{
	if(_nEvt==0) std::cout<<"EcalModuleDigiAlg::execute Start"<<std::endl;
	std::cout<<"Processing event: "<<_nEvt<<std::endl;
	Clear();
	edm4hep::CalorimeterHitCollection* caloVec = w_DigiCaloCol.createAndPut();
	edm4hep::MCRecoCaloAssociationCollection* caloAssoVec = w_CaloAssociationCol.createAndPut();
	edm4hep::SimCalorimeterHitCollection* SimHitCellCol = w_SimCaloTruth.createAndPut();
 	const edm4hep::SimCalorimeterHitCollection* SimHitCol =  r_SimCaloCol.get();
 	std::vector<edm4hep::SimCalorimeterHit> m_simhitCol;
	m_simhitCol.clear();
  if(SimHitCol == 0) 
  {
     std::cout<<"not found SimCalorimeterHitCollection"<< std::endl;
     return StatusCode::SUCCESS;
  }
  if(_Debug>=1) std::cout<<"digi, input sim hit size="<< SimHitCol->size() <<std::endl;

	double totE_bar=0;
	double totE_Sim=0;
	double totE_Digi=0;

	//Merge input simhit(steps) to real simhit(bar).
	m_simhitCol = MergeHits(SimHitCol);
	if(_Debug>=1) std::cout<<"Finish Hit Merge, with Nhit: "<<m_simhitCol.size()<<std::endl;

	std::map< unsigned long int, std::vector<DigiBar> > DigiBlocks;
	std::map< int, std::vector<DigiBar> > DigiLayer;
	DigiBlocks.clear();
	DigiLayer.clear();

	//Loop in SimHit, digitalize SimHit to DigiBar
	for(int i=0;i<m_simhitCol.size();i++){

		edm4hep::SimCalorimeterHit SimHit = m_simhitCol.at(i);
		if(SimHit.getEnergy()<_Eth) continue;
		unsigned long long id = SimHit.getCellID();

		DigiBar hitbar;
		hitbar.cellID = id;
		hitbar.system = m_decoder->get(id, "system");
		hitbar.group   = m_decoder->get(id, "group");
		hitbar.block  = m_decoder->get(id, "block");
		hitbar.slayer = m_decoder->get(id, "slayer");
		hitbar.bar    = m_decoder->get(id, "bar");

		double Lbar = GetBarLength(hitbar);  //NOTE: Is fixed with geometry CRDEcalBarrel_v01. 
		dd4hep::Position hitpos = m_cellIDConverter->position(id);
		dd4hep::Position barpos(10*hitpos.x(), 10*hitpos.y(), 10*hitpos.z());	//cm to mm.
		hitbar.position = barpos;
		hitbar.T1 = 99999; hitbar.T2 = 99999;
		if(_Debug>=2) std::cout<<"SimHit contribution size: "<<SimHit.contributions_size()<<std::endl;

		for(int iCont=0; iCont < SimHit.contributions_size(); ++iCont){
			edm4hep::ConstCaloHitContribution conb = SimHit.getContributions(iCont);
			if( !conb.isAvailable() ) { std::cout<<"EcalModuleDigiAlg  Can not get SimHitContribution: "<<iCont<<std::endl; continue;}
			double en = conb.getEnergy();
			if(en == 0) continue;
			dd4hep::Position steppos(conb.getStepPosition().x, conb.getStepPosition().y, conb.getStepPosition().z);
			dd4hep::Position rpos = steppos-hitbar.position;
			if(_Debug>=3){
				cout<<"Cell Pos: "<<hitbar.position.x()<<'\t'<<hitbar.position.y()<<'\t'<<hitbar.position.z()<<endl;
				cout<<"step pos: "<<steppos.x()<<'\t'<<steppos.y()<<'\t'<<steppos.z()<<endl;
				cout<<"Relative pos: "<<rpos.x()<<'\t'<<rpos.y()<<'\t'<<rpos.z()<<endl;
			}

			//Get digitalized signal(Q1, Q2, T1, T2) from step
			//Define: 1 is left, 2 is right, clockwise direction in phi. 
			TRandom3 rndm(i*_seed+iCont);

			int sign=-999;
			if(hitbar.slayer==1) sign = rpos.z()==0 ? 1 : rpos.z()/fabs(rpos.z());
			if(hitbar.slayer==0) sign = rpos.z()==0 ? 1 : rpos.y()/fabs(rpos.y());

			if(!fabs(sign)) {std::cout<<"ERROR: Wrong bar direction/position!"<<std::endl; continue;}
			double Qi_left = en*exp(-(Lbar/2 + sign*sqrt(rpos.Mag2()))/Latt);	//FIXME: Need to use z rather than magnitude.
			double Qi_right = en*exp(-(Lbar/2 - sign*sqrt(rpos.Mag2()))/Latt);
			if(_Debug>=3){
				cout<<"Bar length: "<<Lbar<<'\t'<<sign*sqrt(rpos.Mag2())<<endl;
				cout<<"Readout Charge: "<<Qi_left<<'\t'<<Qi_right<<endl;
			}
			double Ti_left = -1; int looptime=0;
			while(Ti_left<0){ 
				Ti_left = rndm.Gaus(nMat*(Lbar/2 + sign*sqrt(rpos.Mag2()))/C, Tres);
				looptime++;
				if(looptime>500){ std::cout<<"ERROR: Step "<<iCont<<" can not get a positive left-side time!"<<std::endl; break;}
			}
			if(looptime>500) continue;		
			double Ti_right = -1; looptime=0;
			while(Ti_right<0){ 
				Ti_right = rndm.Gaus(nMat*(Lbar/2 - sign*sqrt(rpos.Mag2()))/C, Tres);
				looptime++;
            if(looptime>500){ std::cout<<"ERROR: Step "<<iCont<<" can not get a positive right-side time!"<<std::endl; break;}
			}
			if(looptime>500) continue;		

			hitbar.Q1 += Qi_left;
			hitbar.Q2 += Qi_right;
			hitbar.T1 = hitbar.T1>Ti_left ? Ti_left : hitbar.T1;
			hitbar.T2 = hitbar.T2>Ti_right ? Ti_right : hitbar.T2;

			//Create truth SimHit(1cm*1cm*1cm cellsize)
			//NOTE: NO cellID for this truth SimHit. 
			dd4hep::Position truthpos = GetCellPos(steppos, hitbar); 
			edm4hep::SimCalorimeterHit simhitTruth = find(SimHitCellCol, truthpos); 
			if(simhitTruth.getCellID()==0){ 
				simhitTruth = SimHitCellCol->create();
				edm4hep::Vector3f m_vec(truthpos.x(), truthpos.y(), truthpos.z());
				simhitTruth.setPosition(m_vec);
				simhitTruth.setCellID(1);
				//SimHitCellCol->push_back(simhitTruth);
				m_simTruth_slayer.push_back(hitbar.slayer);
				m_simTruth_group.push_back(hitbar.group);
				m_simTruth_block.push_back(hitbar.block);
			}
			simhitTruth.addToContributions(conb);
			simhitTruth.setEnergy(simhitTruth.getEnergy()+en );

		}
		totE_bar+=(hitbar.Q1+hitbar.Q2)/2;
		unsigned long int blockID = coder(hitbar);
		int ilayer = 4*(hitbar.group-1) + hitbar.block==1?0:2 + hitbar.slayer+1;  //from 1 to 28. 
		DigiBlocks[blockID].push_back(hitbar);
		DigiLayer[ilayer].push_back(hitbar);

		m_simBar_x.push_back(hitbar.position.x());
		m_simBar_y.push_back(hitbar.position.y());
		m_simBar_z.push_back(hitbar.position.z());
		m_simBar_Q1.push_back(hitbar.Q1);
		m_simBar_Q2.push_back(hitbar.Q2);
		m_simBar_T1.push_back(hitbar.T1);
		m_simBar_T2.push_back(hitbar.T2);
		m_simBar_group.push_back(hitbar.group);
		m_simBar_slayer.push_back(hitbar.slayer);
		m_simBar_block.push_back(hitbar.block);

	}
	t_SimBar->Fill();
	if(_Debug>=1) std::cout<<"End Loop: Merged SimHit!"<<std::endl;
	
	if(_Debug>=1) std::cout<<"TruthSimHit Number: "<<SimHitCellCol->size()<<std::endl;
	for(int iter=0; iter<SimHitCellCol->size();iter++){
		edm4hep::SimCalorimeterHit m_simhit = SimHitCellCol->at(iter);
		if(!m_simhit.isAvailable()) continue;
		m_simTruth_x.push_back(m_simhit.getPosition()[0]);
		m_simTruth_y.push_back(m_simhit.getPosition()[1]);
		m_simTruth_z.push_back(m_simhit.getPosition()[2]);
		m_simTruth_E.push_back(m_simhit.getEnergy());
		totE_Sim+=m_simhit.getEnergy();
	}
	t_SimTruth->Fill();

	if(_Debug>=1) std::cout<<"Block number: "<<DigiBlocks.size()<<std::endl;
	for(auto iter=DigiBlocks.begin(); iter!=DigiBlocks.end();iter++){
		std::vector<edm4hep::CalorimeterHit> m_caloHits = CreateDigiHits(iter->second);
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
	std::cout<<"Total Bar Energy: "<<totE_bar<<std::endl;
	std::cout<<"Total SimHit Energy: "<<totE_Sim<<std::endl;
	std::cout<<"Total DigiHit Energy: "<<totE_Digi<<std::endl;


  _nEvt ++ ;
  return StatusCode::SUCCESS;
}

StatusCode EcalModuleDigiAlg::finalize()
{
	m_wfile->cd();
	t_SimBar->Write();
	t_SimTruth->Write();
	t_Rec->Write();
	m_wfile->Close();

  info() << "Processed " << _nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}

std::vector<edm4hep::SimCalorimeterHit> EcalModuleDigiAlg::MergeHits(const edm4hep::SimCalorimeterHitCollection* m_col){
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



double EcalModuleDigiAlg::GetBarLength(EcalModuleDigiAlg::DigiBar bar){
	//NOTE: Didn't finish reading bar length from geosvc. 
	return 400;
}

dd4hep::Position EcalModuleDigiAlg::GetCellPos(dd4hep::Position pos, EcalModuleDigiAlg::DigiBar bar){
	dd4hep::Position rpos = pos-bar.position;
	TVector3 vec(0,0,0); 
	if(bar.slayer==0) vec.SetXYZ(0, floor(rpos.y()/10)*10+5, 0 );
	if(bar.slayer==1) vec.SetXYZ(0, 0, floor(rpos.z()/10)*10+5 );
	dd4hep::Position relv(vec.x(), vec.y(), vec.z());
	return relv+bar.position;
}

edm4hep::SimCalorimeterHit EcalModuleDigiAlg::find(edm4hep::SimCalorimeterHitCollection* m_col, dd4hep::Position pos){
   for(int i=0;i<m_col->size();i++){
    edm4hep::SimCalorimeterHit hit = m_col->at(i);
		dd4hep::Position ipos(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z);
		if(ipos==pos) return hit;
	}
  edm4hep::SimCalorimeterHit hit;
  hit.setCellID(0);
  return hit;
}

edm4hep::SimCalorimeterHit EcalModuleDigiAlg::find(std::vector<edm4hep::SimCalorimeterHit> m_col, unsigned long long cellid){
   for(int i=0;i<m_col.size();i++){
		edm4hep::SimCalorimeterHit hit=m_col.at(i);
		if(hit.getCellID() == cellid) return hit;
	}
	edm4hep::SimCalorimeterHit hit;
	hit.setCellID(0);
	return hit;
}

unsigned long int EcalModuleDigiAlg::coder(EcalModuleDigiAlg::DigiBar bar){
	return 10099*bar.group + 31*bar.block;
}


std::vector<edm4hep::CalorimeterHit> EcalModuleDigiAlg::CreateDigiHits(std::vector<EcalModuleDigiAlg::DigiBar> m_block){
   int Nbars = m_block.size();
	if(_Debug>=2) std::cout<<"Bar number in this block: "<<Nbars<<std::endl;
   if(Nbars==0){std::cout<<"No bars in this block!"<<std::endl; exit(0);}

   std::vector<edm4hep::CalorimeterHit> m_digiCol;
	m_digiCol.clear();
	float rotAngle = 2*PI/4.;
   float Edes[2][50]={0};           //Deposited energy in slayer 0/1. 
   float totE[2]={0};
   for(int ibar=0;ibar<Nbars;ibar++){
      float m_En = r_cali*(m_block[ibar].Q1+m_block[ibar].Q2)/2;
      Edes[m_block[ibar].slayer][m_block[ibar].bar] = m_En;
      totE[m_block[ibar].slayer] += m_En;
   }

	if(_Debug>=2) std::cout<<"totE in sub-layer 0/1: "<<totE[0]<<'\t'<<totE[1]<<std::endl;

   TVector3 m_vec(0,0,0);
   dd4hep::Position m_pos(0,0,0);
   for(int ibar=0;ibar<Nbars;ibar++){
      if(m_block[ibar].slayer==1) continue;     //Outer loop in slayer==0.
      DigiBar bar_s0 = m_block[ibar];

      m_vec.SetXYZ(bar_s0.position.x(), bar_s0.position.y(), bar_s0.position.z());
		m_vec.RotateZ(rotAngle);
      bar_s0.position.SetXYZ(m_vec.x(), m_vec.y(), m_vec.z());

      for(int jbar=0;jbar<Nbars;jbar++){
         if(m_block[jbar].slayer==0) continue;  //Inner loop in slayer==1.
         DigiBar bar_s1 = m_block[jbar];
			double Lbar = GetBarLength(bar_s0);
			double dis_s0 = C*(bar_s0.T1-bar_s0.T2)/(2*nMat);
			double dis_s1 = C*(bar_s1.T1-bar_s1.T2)/(2*nMat);
         m_vec.SetXYZ(bar_s1.position.x(), bar_s1.position.y(), bar_s1.position.z());
			m_vec.RotateZ(rotAngle);
         bar_s1.position.SetXYZ(m_vec.x(), m_vec.y(), m_vec.z());

			if(fabs( dis_s0+bar_s0.position.x() - bar_s1.position.x())>_DeltaZth*sqrt(2)*C*Tres/(2*nMat) ) continue;	//If cross hit is out of _Zth sigma. 
			if(fabs( dis_s1+bar_s1.position.z() - bar_s0.position.z())>_DeltaZth*sqrt(2)*C*Tres/(2*nMat) ) continue;	//If cross hit is out of _Zth sigma. 

         TVector3 p_hit(bar_s1.position.x(), (bar_s0.position.y()+bar_s1.position.y())/2., bar_s0.position.z() );
			p_hit.RotateZ(-rotAngle);
         edm4hep::Vector3f m_vec3f(p_hit.x(), p_hit.y(), p_hit.z());
         float m_En = Edes[bar_s0.slayer][bar_s0.bar]*Edes[bar_s1.slayer][bar_s1.bar]/totE[bar_s1.slayer]
                           + Edes[bar_s1.slayer][bar_s1.bar]*Edes[bar_s0.slayer][bar_s0.bar]/totE[bar_s0.slayer];
         edm4hep::CalorimeterHit hit;
         hit.setCellID(0);
         hit.setPosition(m_vec3f);
         hit.setEnergy(m_En);
         m_digiCol.push_back(hit);
      }
   }
   return m_digiCol;
}

void EcalModuleDigiAlg::Clear(){
	m_simBar_x.clear();
	m_simBar_y.clear();
	m_simBar_z.clear();
	m_simBar_T1.clear();
	m_simBar_T2.clear();
	m_simBar_Q1.clear();
	m_simBar_Q2.clear();
	m_simBar_slayer.clear();
	m_simBar_group.clear();
	m_simBar_block.clear();
	m_simTruth_x.clear();
	m_simTruth_y.clear();
	m_simTruth_z.clear();
	m_simTruth_E.clear();
	m_simTruth_block.clear();
	m_simTruth_group.clear();
	m_simTruth_slayer.clear();
	m_Rec_x.clear();
	m_Rec_y.clear();
	m_Rec_z.clear();
	m_Rec_E.clear();
}
