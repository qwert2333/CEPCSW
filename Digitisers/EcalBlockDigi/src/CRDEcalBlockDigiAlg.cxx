/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
// Unit in code: mm, ns. 
// NOTE: This digitialization highly matches detector geometry CRDEcalBarrel_v01. 
#include "CRDEcalDigiEDM.h"
#include "CRDEcalBlockDigiAlg.h"
#include "CRDEcalPreRecAlg.h"

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
DECLARE_COMPONENT( CRDEcalBlockDigiAlg )

CRDEcalBlockDigiAlg::CRDEcalBlockDigiAlg(const std::string& name, ISvcLocator* svcLoc)
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

StatusCode CRDEcalBlockDigiAlg::initialize()
{

	std::string s_outfile = _filename;
	m_wfile = new TFile(s_outfile.c_str(), "recreate");
	t_SimCont = new TTree("SimStep", "SimStep");
	t_PreRec = new TTree("SimBarHit", "SimBarHit");
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

	t_PreRec->Branch("simBar_x", &m_simBar_x);
	t_PreRec->Branch("simBar_y", &m_simBar_y);
	t_PreRec->Branch("simBar_z", &m_simBar_z);
	t_PreRec->Branch("simBar_T1", &m_simBar_T1);
	t_PreRec->Branch("simBar_T2", &m_simBar_T2);
	t_PreRec->Branch("simBar_Q1", &m_simBar_Q1);
	t_PreRec->Branch("simBar_Q2", &m_simBar_Q2);
	t_PreRec->Branch("simBar_Nstep", &m_simBar_Nstep);
   t_PreRec->Branch("PreRec_showerX", &m_PreRec_showerX);
   t_PreRec->Branch("PreRec_showerY", &m_PreRec_showerY);
   t_PreRec->Branch("PreRec_showerZ", &m_PreRec_showerZ);
   t_PreRec->Branch("PreRec_showerE", &m_PreRec_showerE);
   t_PreRec->Branch("PreRec_showerT1", &m_PreRec_showerT1);
   t_PreRec->Branch("PreRec_showerT2", &m_PreRec_showerT2);
   t_PreRec->Branch("PreRec_Nshower", &m_PreRec_Nshower);
   t_PreRec->Branch("PreRec_Ncluster", &m_PreRec_Ncluster);
   t_PreRec->Branch("PreRec_Clus_ScndM", &m_Clus_ScndM);
   t_PreRec->Branch("PreRec_shower_Nbars", &m_shower_Nbars);
   t_PreRec->Branch("PreRec_shower_barx", &m_shower_barx);
   t_PreRec->Branch("PreRec_shower_bary", &m_shower_bary);
   t_PreRec->Branch("PreRec_shower_barz", &m_shower_barz);
   t_PreRec->Branch("PreRec_shower_barE", &m_shower_barE);


	std::cout<<"CRDEcalBlockDigiAlg::m_scale="<<m_scale<<std::endl;
	m_geosvc = service<IGeomSvc>("GeoSvc");
	if ( !m_geosvc )  throw "CRDEcalBlockDigiAlg :Failed to find GeoSvc ...";
	dd4hep::Detector* m_dd4hep = m_geosvc->lcdd();
	if ( !m_dd4hep )  throw "CRDEcalBlockDigiAlg :Failed to get dd4hep::Detector ...";
	m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep);
	const std::string name_readout = "CaloHitsCollection";
	m_decoder = m_geosvc->getDecoder(name_readout);
	if (!m_decoder) {
		error() << "Failed to get the decoder. " << endmsg;
		return StatusCode::FAILURE;
	}

	rndm.SetSeed(_seed);
	std::cout<<"CRDEcalBlockDigiAlg::initialize"<<std::endl;
	return GaudiAlgorithm::initialize();
}

StatusCode CRDEcalBlockDigiAlg::execute()
{
	if(_nEvt==0) std::cout<<"CRDEcalBlockDigiAlg::execute Start"<<std::endl;
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
	std::vector<CRDEcalDigiEDM::DigiBar> barCol; barCol.clear();


	//Loop in SimHit, digitalize SimHit to DigiBar
	for(int i=0;i<m_simhitCol.size();i++){


		edm4hep::SimCalorimeterHit SimHit = m_simhitCol.at(i);
		if(SimHit.getEnergy()<_Eth) continue;


		CRDEcalDigiEDM::DigiBar hitbar;
		unsigned long long id = SimHit.getCellID();
		hitbar.cellID = id;
		hitbar.bar    = m_decoder->get(id, "slice");

		double Lbar = GetBarLength(hitbar);  //NOTE: Is fixed with geometry CRDEcalBarrel_v01. 
		dd4hep::Position hitpos = m_cellIDConverter->position(id);
		dd4hep::Position barpos(10*hitpos.x(), 10*hitpos.y(), 10*hitpos.z());	//cm to mm.
		hitbar.position = barpos;
		hitbar.T1 = 99999; hitbar.T2 = 99999;
		if(_Debug>=2) std::cout<<"SimHit contribution size: "<<SimHit.contributions_size()<<std::endl;


		std::vector<CRDEcalDigiEDM::StepDigiOut> DigiLvec; DigiLvec.clear();
		std::vector<CRDEcalDigiEDM::StepDigiOut> DigiRvec; DigiRvec.clear();
		double totQ1 = 0;
		double totQ2 = 0;

		//Loop in all SimHitContribution(G4Step). 
		for(int iCont=0; iCont < SimHit.contributions_size(); ++iCont){
			edm4hep::ConstCaloHitContribution conb = SimHit.getContributions(iCont);
			if( !conb.isAvailable() ) { std::cout<<"CRDEcalBlockDigiAlg  Can not get SimHitContribution: "<<iCont<<std::endl; continue;}

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

			if(_Debug>=2){
				cout<<"Cell Pos: "<<hitbar.position.x()<<'\t'<<hitbar.position.y()<<'\t'<<hitbar.position.z()<<endl;
				cout<<"step pos: "<<steppos.x()<<'\t'<<steppos.y()<<'\t'<<steppos.z()<<endl;
				cout<<"Relative pos: "<<rpos.x()<<'\t'<<rpos.y()<<'\t'<<rpos.z()<<endl;
			}

			//Get digitalized signal(Q1, Q2, T1, T2) from step
			//Define: 1 is left, 2 is right, clockwise direction in phi. 


			double Qi_left = en*exp(-(Lbar/2 + steppos.y() )/Latt);	//FIXME: Need to use z rather than magnitude.
			double Qi_right = en*exp(-(Lbar/2 - steppos.y() )/Latt);

			if(_Debug>=2){
				cout<<Qi_left<<'\t'<<Qi_right<<endl;
				cout<<Lbar<<'\t'<<sqrt(rpos.Mag2())<<endl;
			}


			double Ti_left = -1; int looptime=0;
			while(Ti_left<0){ 
				Ti_left = Tinit + rndm.Gaus(nMat*(Lbar/2 + steppos.y() )/C, Tres); //FIXME: same as Qi.
				looptime++;
				if(looptime>500){ std::cout<<"ERROR: Step "<<iCont<<" can not get a positive left-side time!"<<std::endl; break;}
			}
			if(looptime>500) continue;		
			double Ti_right = -1; looptime=0;
			while(Ti_right<0){ 
				Ti_right = Tinit + rndm.Gaus(nMat*(Lbar/2 - steppos.y() )/C, Tres); //FIXME: same as Qi.
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
		
		m_simBar_x.push_back(hitbar.position.x());
		m_simBar_y.push_back(hitbar.position.y());
		m_simBar_z.push_back(hitbar.position.z());
		m_simBar_Q1.push_back(hitbar.Q1);
		m_simBar_Q2.push_back(hitbar.Q2);
		m_simBar_T1.push_back(hitbar.T1);
		m_simBar_T2.push_back(hitbar.T2);
		m_simBar_Nstep.push_back(SimHit.contributions_size());

		barCol.push_back(hitbar);
	}
	t_SimCont->Fill();
	if(_Debug>=1) std::cout<<"End Loop: Bar Digitalization!"<<std::endl;
	if(_Debug>=1) std::cout<<"Fired bar number: "<<barCol.size()<<std::endl;

	
	//1D clustering 
   CRDEcalPreRecAlg preRec;
   preRec.Eth_SeedWithNeigh  = _Eth_SeedWithNeigh;
   preRec.Eth_ShowerWithTot  = _Eth_ShowerWithTot;
   preRec.Eth_ClusterWithTot = _Eth_ClusterWithTot;
   preRec.Sth_split          = _Sth_split;
   preRec.Eth_SeedAbs        =  _Eth_SeedAbs;
   preRec.Eth_ShowerAbs      = _Eth_ShowerAbs;
   preRec.Eth_ClusAbs        = _Eth_ClusAbs;

	std::vector<CRDEcalDigiEDM::BarCollection> barShowerCol; barShowerCol.clear();
	std::vector<CRDEcalDigiEDM::BarCluster> barCluster = preRec.Clustering(barCol);
   for(int i=0;i<barCluster.size();i++){
      double _Etot=0;
      for(int i=0;i<barCol.size();i++) _Etot+=barCol[i].getEnergy();
      if(barCluster[i].getE()/_Etot < _Eth_ClusterWithTot || barCluster[i].getE()<_Eth_ClusAbs ) continue;
      std::vector<CRDEcalDigiEDM::BarCollection> showers = preRec.ClusterSplitting(barCluster[i]);
      if(showers.size()==0) continue;
      barShowerCol.insert(barShowerCol.end(), showers.begin(), showers.end());
      m_Clus_ScndM.push_back(barCluster[i].getScndMoment());
   }

   for(int i=0;i<barShowerCol.size();i++){
      m_PreRec_showerX.push_back(barShowerCol[i].getPos().x());
      m_PreRec_showerY.push_back(barShowerCol[i].getPos().y());
      m_PreRec_showerZ.push_back(barShowerCol[i].getPos().z());
      m_PreRec_showerE.push_back(barShowerCol[i].getE());
      m_PreRec_showerT2.push_back(barShowerCol[i].getT2());
      m_PreRec_showerT1.push_back(barShowerCol[i].getT1());

      m_shower_Nbars.push_back(barShowerCol[i].Bars.size());
      for(int ii=0;ii<barShowerCol[i].Bars.size();ii++){
         m_shower_barx.push_back(barShowerCol[i].Bars[ii].position.x());
         m_shower_bary.push_back(barShowerCol[i].Bars[ii].position.y());
         m_shower_barz.push_back(barShowerCol[i].Bars[ii].position.z());
         m_shower_barE.push_back(barShowerCol[i].Bars[ii].getEnergy());
      }
   }

	m_PreRec_Nshower = barShowerCol.size();
	m_PreRec_Ncluster = barCluster.size();

	t_PreRec->Fill();

	std::cout<<"Total Bar Energy: "<<totE_bar<<std::endl;
	std::cout<<"Total DigiHit Energy: "<<totE_Digi<<std::endl;


  _nEvt ++ ;
  return StatusCode::SUCCESS;
}

StatusCode CRDEcalBlockDigiAlg::finalize()
{
	m_wfile->cd();
	t_SimCont->Write();
	t_PreRec->Write();
	m_wfile->Close();

  info() << "Processed " << _nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}

std::vector<edm4hep::SimCalorimeterHit> CRDEcalBlockDigiAlg::MergeHits(const edm4hep::SimCalorimeterHitCollection* m_col){
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



double CRDEcalBlockDigiAlg::GetBarLength(CRDEcalDigiEDM::DigiBar bar){
	//NOTE: Didn't finish reading bar length from geosvc. 
	return 400.;
}

dd4hep::Position CRDEcalBlockDigiAlg::GetCellPos(dd4hep::Position pos, CRDEcalDigiEDM::DigiBar bar){
	dd4hep::Position rpos = pos-bar.position;
	TVector3 vec(0,0,0); 
	vec.SetXYZ(0, 0, floor(rpos.z()/10)*10+5 );
	dd4hep::Position relv(vec.x(), vec.y(), vec.z());
	return relv+bar.position;
}

edm4hep::SimCalorimeterHit CRDEcalBlockDigiAlg::find(edm4hep::SimCalorimeterHitCollection* m_col, dd4hep::Position pos){
   for(int i=0;i<m_col->size();i++){
      edm4hep::SimCalorimeterHit hit = m_col->at(i);
		dd4hep::Position ipos(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z);
		if(ipos==pos) return hit;
	}
   edm4hep::SimCalorimeterHit hit;
   hit.setCellID(0);
   return hit;
}

edm4hep::SimCalorimeterHit CRDEcalBlockDigiAlg::find(std::vector<edm4hep::SimCalorimeterHit> m_col, unsigned long long cellid){
   for(int i=0;i<m_col.size();i++){
		edm4hep::SimCalorimeterHit hit=m_col.at(i);
		if(hit.getCellID() == cellid) return hit;
	}
	edm4hep::SimCalorimeterHit hit;
	hit.setCellID(0);
	return hit;
}


void CRDEcalBlockDigiAlg::Clear(){
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
	m_simBar_Nstep.clear();

   m_PreRec_showerX.clear();
   m_PreRec_showerY.clear();
   m_PreRec_showerZ.clear();
   m_PreRec_showerE.clear();
   m_PreRec_showerT1.clear();
   m_PreRec_showerT2.clear();
   m_PreRec_Nshower=-999;
   m_PreRec_Ncluster=-999;
   m_Clus_ScndM.clear();
   m_shower_Nbars.clear();
   m_shower_barx.clear();
   m_shower_bary.clear();
   m_shower_barz.clear();
   m_shower_barE.clear();

}
