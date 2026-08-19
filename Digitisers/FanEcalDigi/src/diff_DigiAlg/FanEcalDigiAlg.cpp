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

#include "CLUEAlgo.cc"
#include "CLUEAlgo.h"

#define C 299.79  // unit: mm/ns
#define PI 3.141592653

using namespace std;
using namespace dd4hep;

DECLARE_COMPONENT( FanEcalDigiAlg )

FanEcalDigiAlg::FanEcalDigiAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0),
    recEvt(0)
{
	// Input collections
	declareProperty("SimCaloHitCollection", r_SimCaloCol, "Handle of the Input SimCaloHit collection");
	declareProperty("CaloHitContCol", r_CaloContCol, "Handle of the Input SimCaloHit contribution collection");
//	declareProperty("MCParticleCol", m_mcParCol, "MCParticle collection (input)");
	
	// Output collections
	declareProperty("CaloHitCollection", w_DigiCaloCol, "Handle of Digi CaloHit collection");
	declareProperty("CaloAssociationCollection", w_CaloAssociationCol, "Handle of CaloAssociation collection");
        declareProperty("CaloMCPAssociationCollection", w_MCPCaloAssociationCol, "Handle of CaloAssociation collection");
}

StatusCode FanEcalDigiAlg::initialize()
{
  if(_writeNtuple){
   	std::string s_outfile = _filename;
	  m_wfile = new TFile(s_outfile.c_str(), "recreate");
	  t_SimStep = new TTree("SimStep", "SimStep");
	  t_SimBar = new TTree("SimBar", "SimBar");
	  t_ClusBar = new TTree("ClusBar", "ClusBar");
	  t_MCPdata = new TTree("MCPdata", "MCPdata");
	  t_MCHitInfo = new TTree("MCHitInfo", "MCHitInfo");
          
	  t_MCHitInfo->Branch("MChit_shower_pdg", &m_MChit_shower_pdg);
 	  t_MCHitInfo->Branch("MChit_crystal", &m_MChit_crystal);
	  t_MCHitInfo->Branch("MChit_module", &m_MChit_module);
	  t_MCHitInfo->Branch("MChit_E", &m_MChit_E);
	  t_MCHitInfo->Branch("MChit_x", &m_MChit_x);
	  t_MCHitInfo->Branch("MChit_y", &m_MChit_y);
	  t_MCHitInfo->Branch("MChit_z", &m_MChit_z);

	  t_MCHitInfo->Branch("MChit_primary_pdg", &m_MChit_primary_pdg);
	  t_MCHitInfo->Branch("MChit_primary_E", &m_MChit_primary_E);
	  t_MCHitInfo->Branch("MChit_primary_x", &m_MChit_primary_x);
	  t_MCHitInfo->Branch("MChit_primary_y", &m_MChit_primary_y);
	  t_MCHitInfo->Branch("MChit_primary_z", &m_MChit_primary_z);

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

	t_ClusBar->Branch("clusBar_x", &m_clusBar_x);
        t_ClusBar->Branch("clusBar_y", &m_clusBar_y);
        t_ClusBar->Branch("clusBar_z", &m_clusBar_z);
        t_ClusBar->Branch("clusBar_T1", &m_clusBar_T1);
        t_ClusBar->Branch("clusBar_T2", &m_clusBar_T2);
        t_ClusBar->Branch("clusBar_Q1", &m_clusBar_Q1);
        t_ClusBar->Branch("clusBar_Q2", &m_clusBar_Q2);
        t_ClusBar->Branch("clusBar_module", &m_clusBar_module);
        t_ClusBar->Branch("clusBar_crystal", &m_clusBar_crystal);
 //       t_ClusBar->Branch("Photon_E_MC", &m_Photon_E_MC);

	t_MCPdata->Branch("MCP_pdg", &m_MCP_pdg);
	t_MCPdata->Branch("MCP_charge", &m_MCP_charge);
	t_MCPdata->Branch("MCP_crystal", &m_MCP_crystal);
	t_MCPdata->Branch("MCP_module", &m_MCP_module);
	t_MCPdata->Branch("MCP_endpoint_x", &m_MCP_endpoint_x);
	t_MCPdata->Branch("MCP_endpoint_y", &m_MCP_endpoint_y);
	t_MCPdata->Branch("MCP_endpoint_z", &m_MCP_endpoint_z);
  }

	std::cout<<"FanEcalDigiAlg::m_scale="<<m_scale<<std::endl;
	m_geosvc = service<IGeomSvc>("GeomSvc");
	if ( !m_geosvc )  throw "FanEcalDigiAlg :Failed to find GeoSvc ...";
	dd4hep::Detector* m_dd4hep = m_geosvc->lcdd();
	if ( !m_dd4hep )  throw "FanEcalDigiAlg :Failed to get dd4hep::Detector ...";
	m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep);
        m_decoder = m_geosvc->getDecoder(_readoutName);
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
	edm4hep::CalorimeterHitCollection* caloVec = w_DigiCaloCol.createAndPut();
	edm4hep::MCRecoCaloAssociationCollection* caloAssoVec = w_CaloAssociationCol.createAndPut();
 	edm4hep::MCRecoCaloParticleAssociationCollection* caloMCPAssoVec = w_MCPCaloAssociationCol.createAndPut(); 

	const edm4hep::SimCalorimeterHitCollection* SimHitCol =  r_SimCaloCol.get();
        std::vector<edm4hep::SimCalorimeterHit> m_simhitCol; m_simhitCol.clear();

    //===========MCdata=============================
/*	 auto mcCol = m_mcParCol.get();
   	int m_nParticles = 0;
   	std::cout<<mcCol->size()<<std::endl;
   	for ( auto p : *mcCol ) {
	    m_Photon_E_MC.push_back(p.getEnergy());
	}
*/
    //================Digitization=================
    if(SimHitCol == 0) 
    {
       std::cout<<"not found SimCalorimeterHitCollection"<< std::endl;
       return StatusCode::SUCCESS;
    }
    if(_Debug>=1) std::cout<<"digi, input sim hit size="<< SimHitCol->size() <<std::endl;
  
    double totE_bar=0;

  /*  //================get each hit info===================
    for(int i_hitcon=0; i_hitcon<SimHitCol->size(); i_hitcon++)
   {  cout<<"This is "<<i_hitcon<<" hit in total "<<SimHitCol->size()<<" sim hits."<<endl;
      edm4hep::SimCalorimeterHit SimCalHit = SimHitCol->at(i_hitcon);
      if(SimCalHit.getEnergy()<_Eth) continue;
      for(int con_size=0; con_size<SimCalHit.contributions_size(); con_size++)
      {  
	edm4hep::CaloHitContribution calo_contribution = SimCalHit.getContributions(con_size);
	if(calo_contribution.getEnergy()>0){
	  unsigned long long hit_id = SimCalHit.getCellID();
          m_MChit_crystal.push_back(m_decoder->get(hit_id,"crystal"));
          m_MChit_module.push_back(m_decoder->get(hit_id,"module"));
          m_MChit_E.push_back(calo_contribution.getEnergy());
          m_MChit_x.push_back(calo_contribution.getStepPosition().x);
          m_MChit_y.push_back(calo_contribution.getStepPosition().y);
          m_MChit_z.push_back(calo_contribution.getStepPosition().z);
	  m_MChit_shower_pdg.push_back(calo_contribution.getPDG());// PDG code of the shower particle that caused this contribution.
 
          m_MChit_primary_pdg.push_back(calo_contribution.getParticle().getPDG());//primary MCParticle that caused the shower responsible for this contribution to the hit
	  const edm4hep::Vector3f& momentumVector = calo_contribution.getParticle().getMomentum();
	  double primary_mass = calo_contribution.getParticle().getMass();
	  float primary_E = sqrt(momentumVector[0]*momentumVector[0] + momentumVector[1]*momentumVector[1] + momentumVector[2]*momentumVector[2] + primary_mass*primary_mass);
          m_MChit_primary_E.push_back(primary_E);
	  m_MChit_primary_x.push_back(calo_contribution.getParticle().getEndpoint()[0]);
	  m_MChit_primary_y.push_back(calo_contribution.getParticle().getEndpoint()[1]);
	  m_MChit_primary_z.push_back(calo_contribution.getParticle().getEndpoint()[2]);   

          cout<<"Info from SimHitContribution: "<<endl;
          cout<<"CellID: "<<hit_id<<endl;
          cout<<"Energy: "<<calo_contribution.getEnergy()<<endl;
          cout<<"Position: "<<calo_contribution.getStepPosition().x<<'\t'<<calo_contribution.getStepPosition().y<<'\t'<<calo_contribution.getStepPosition().z<<endl;
          cout<<"Particle PDG: "<<calo_contribution.getParticle().getPDG()<<endl;
          cout<<"Particle Energy: "<<primary_E<<endl;
          cout<<"PDG: "<<calo_contribution.getPDG()<<endl;
	}
      }
   }
    t_MCHitInfo->Fill();
*/
   std::map<std::tuple<int, double, double, double>, std::vector<edm4hep::SimCalorimeterHit>>  m_HitCollectionMap;
   std::vector<std::pair<std::tuple<int, double, double, double>, edm4hep::SimCalorimeterHit>> m_MapHits; m_MapHits.clear();
   HitClassification(*SimHitCol, m_HitCollectionMap);
   SplitCaloHits(m_HitCollectionMap, m_MapHits);
   
   for(const auto&entry : m_HitCollectionMap) {
	int incoming_pdg = std::get<0>(entry.first);
	double MC_px = std::get<1>(entry.first);
 	double MC_py = std::get<2>(entry.first);
	double MC_pz = std::get<3>(entry.first);
	cout << "Key: PDG=" << std::get<0>(entry.first) << ", Position=("
                  << std::get<1>(entry.first) << ", " << std::get<2>(entry.first) << ", " << std::get<3>(entry.first) << ")\n";
	const std::vector<edm4hep::SimCalorimeterHit>& hitCollection = entry.second;
	cout<<"this group pdg is: "<<incoming_pdg<<endl;
	for (const auto& Hits : hitCollection) {
	    double energy = Hits.getEnergy();
	    std::cout << "Incoming PDG: " << incoming_pdg << std::endl;
	    std::cout << "Energy: " << energy << std::endl;
	}
	
    }
    //======================= digitalize SimHit to DigiBar from same incoming particle===============
    std::vector<CaloBar> MCHitVec; MCHitVec.clear();
    for(int imch=0; imch<m_MapHits.size(); imch++){
	auto MC_hit = m_MapHits.at(imch);
	if(MC_hit.getEnergy()<_Eth) continue;
	
    }

    //==================================================== 
    MergeHits(*SimHitCol, m_simhitCol);
    std::vector<CaloBar> HitBarVec; HitBarVec.clear(); 
    //Loop in SimHit, digitalize SimHit to DigiBar
    for(int i=0;i<m_simhitCol.size();i++){

	auto SimHit = m_simhitCol.at(i);
	if(SimHit.getEnergy()<_Eth) continue;

	unsigned long long id = SimHit.getCellID();
	CaloBar hitbar;
	hitbar.setcellID(id);
	hitbar.setcellID(m_decoder->get(id, "system"), 
			 m_decoder->get(id, "module"), 
			 m_decoder->get(id, "crystal"));
	dd4hep::Position hitpos = m_cellIDConverter->position(id);
	dd4hep::Position barpos(10*hitpos.x(), 10*hitpos.y(), 10*hitpos.z());	//cm to mm.
	hitbar.setPosition(barpos);

	MCParticleToEnergyWeightMap MCPEnMap; MCPEnMap.clear();

	std::vector<CaloStep> DigiLvec; DigiLvec.clear();
	std::vector<CaloStep> DigiRvec; DigiRvec.clear();
	double totQ1 = 0;
	double totQ2 = 0;

	//Loop in all SimHitContribution(G4Step). 
	for(int iCont=0; iCont < SimHit.contributions_size(); ++iCont){
		auto conb = SimHit.getContributions(iCont);
	        if( !conb.isAvailable() ) { std::cout<<"FanEcalDigiAlg  Can not get SimHitContribution: "<<iCont<<std::endl; continue;}

		double en = conb.getEnergy();
		if(en == 0) continue;

		auto mcp = conb.getParticle();
		MCPEnMap[mcp] += en;

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
//==================create calo association==================================================================    
        edm4hep::Vector3f m_pos(hitbar.getPosition().X(), hitbar.getPosition().Y(), hitbar.getPosition().Z());
        auto digiHit1 = caloVec->create();
        digiHit1.setCellID(hitbar.getcellID());
    	digiHit1.setEnergy(hitbar.getQ1());
    	digiHit1.setTime(hitbar.getT1());
    	digiHit1.setPosition(m_pos);
    	auto digiHit2 = caloVec->create();
    	digiHit2.setCellID(hitbar.getcellID());
    	digiHit2.setEnergy(hitbar.getQ2());
    	digiHit2.setTime(hitbar.getT2());
    	digiHit2.setPosition(m_pos);

    	//SimHit - CaloHit association
    	auto rel1 = caloAssoVec->create();
    	rel1.setRec(digiHit1);
    	rel1.setSim(SimHit);
    	rel1.setWeight( hitbar.getQ1()/(hitbar.getQ1()+hitbar.getQ2()) );
    	auto rel2 = caloAssoVec->create();
    	rel2.setRec(digiHit2);
    	rel2.setSim(SimHit);
    	rel2.setWeight( hitbar.getQ2()/(hitbar.getQ1()+hitbar.getQ2()) );

    	//MCParticle - CaloHit association
    	for(auto iter : MCPEnMap){
    		auto rel_MCP1 = caloMCPAssoVec->create();
      	 	rel_MCP1.setRec(digiHit1);
      		rel_MCP1.setSim(iter.first);
      		rel_MCP1.setWeight(iter.second/digiHit1.getEnergy());
      		auto rel_MCP2 = caloMCPAssoVec->create();
      		rel_MCP2.setRec(digiHit2);
      		rel_MCP2.setSim(iter.first);
      		rel_MCP2.setWeight(iter.second/digiHit2.getEnergy());      
    	}
//=================================================================================================
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
         cout<<"caloMCPAssoVec->size : "<<caloMCPAssoVec->size()<<endl;
     for(int con_caloMCPAssoVec = 0; con_caloMCPAssoVec < caloMCPAssoVec->size(); con_caloMCPAssoVec++){
	 auto calo_hit = caloMCPAssoVec->at(con_caloMCPAssoVec).getRec();
         unsigned long long hit_id = calo_hit.getCellID();
         cout<<"cell ID: "<<hit_id<<endl;
	 auto calo_mcp = caloMCPAssoVec->at(con_caloMCPAssoVec).getSim();
	 cout<<"PDG: "<<calo_mcp.getPDG()<<endl;
	 cout<<"endpoint X: "<<calo_mcp.getEndpoint().x<<"	endpoint Y: "<<calo_mcp.getEndpoint().y<<"      endpoint Z: "<<calo_mcp.getEndpoint().z<<endl;  
	 m_MCP_pdg.push_back(calo_mcp.getPDG());
	 m_MCP_charge.push_back(calo_mcp.getCharge());
	 m_MCP_crystal.push_back(m_decoder->get(hit_id, "crystal"));
	 m_MCP_module.push_back(m_decoder->get(hit_id, "module"));
	 m_MCP_endpoint_x.push_back(calo_mcp.getEndpoint().x);
	 m_MCP_endpoint_y.push_back(calo_mcp.getEndpoint().y);
	 m_MCP_endpoint_z.push_back(calo_mcp.getEndpoint().z);
      }
      t_MCPdata->Fill();
  //====================Reconstruction=================================
  std::vector<CaloCluster> ClusVec; ClusVec.clear();
  NeighborClustering(HitBarVec,ClusVec);

  double tot_clus_E = 0;
  m_Ncluster = ClusVec.size();
  for(int ic=0; ic<m_Ncluster; ic++){

  /*    m_clus_phi_start.push_back(ClusVec[ic].getPhiStart());
      m_clus_phi.push_back(ClusVec[ic].getPhiStartFit());
      m_clus_Z.push_back(ClusVec[ic].getZ());
      m_clus_aphi_E.push_back(ClusVec[ic].getAvePhi(true));
      m_clus_aphi.push_back(ClusVec[ic].getAvePhi(false));
      m_clus_E.push_back(ClusVec[ic].getEnergy());
      m_clus_chi2.push_back(ClusVec[ic].getFitChi2());
      m_clus_alpha.push_back(ClusVec[ic].getFitAlpha());
      m_clus_beta.push_back(ClusVec[ic].getFitBeta());
      m_clus_Nbars.push_back(ClusVec[ic].getBars().size());
*/
      vector<CaloBar> n_hitBars;  n_hitBars.clear();
      n_hitBars = ClusVec[ic].getBars();
      for(int ibar=0;ibar<n_hitBars.size();ibar++){
	m_clusBar_x.push_back(n_hitBars[ibar].getPosition().x());
        m_clusBar_y.push_back(n_hitBars[ibar].getPosition().y());
        m_clusBar_z.push_back(n_hitBars[ibar].getPosition().z());
        m_clusBar_Q1.push_back(n_hitBars[ibar].getQ1());
        m_clusBar_Q2.push_back(n_hitBars[ibar].getQ2());
        m_clusBar_T1.push_back(n_hitBars[ibar].getT1());
        m_clusBar_T2.push_back(n_hitBars[ibar].getT2());
    	m_clusBar_module.push_back(n_hitBars[ibar].getModule());
        m_clusBar_crystal.push_back(n_hitBars[ibar].getCrystal());
      } 

  }
//  t_RecClus->Fill();
  t_ClusBar->Fill();

  _nEvt ++ ;
  m_simhitCol.clear();
  return StatusCode::SUCCESS;
}

StatusCode FanEcalDigiAlg::finalize()
{
	m_wfile->cd();
	t_SimStep->Write();
	t_SimBar->Write();
	t_ClusBar->Write();
	t_MCPdata->Write();
	t_MCHitInfo->Write();
	m_wfile->Close();
	delete m_wfile, t_SimStep, t_SimBar;

    info() << "Processed " << _nEvt << " events " << endmsg;
    cout<<"there are "<<double(recEvt)/_nEvt<<" evts reconstructed."<<endl;
    return GaudiAlgorithm::finalize();
}

StatusCode FanEcalDigiAlg::HitClassification(const edm4hep::SimCalorimeterHitCollection& m_col, std::map<std::tuple<int, double, double, double>, std::vector<edm4hep::SimCalorimeterHit>>& HitCollectionMap ){
	
   HitCollectionMap.clear();

   for(int i_hitcon=0; i_hitcon<m_col.size(); i_hitcon++){
      edm4hep::SimCalorimeterHit SimCalHit = m_col[i_hitcon];
      if(SimCalHit.getEnergy()<_Eth) continue;

      for(int con_size=0; con_size<SimCalHit.contributions_size(); con_size++){
        edm4hep::CaloHitContribution calo_contribution = SimCalHit.getContributions(con_size);
	int incoming_pdg = calo_contribution.getParticle().getPDG();
	double positionX = calo_contribution.getParticle().getEndpoint()[0];
        double positionY = calo_contribution.getParticle().getEndpoint()[1];
        double positionZ = calo_contribution.getParticle().getEndpoint()[2];

	auto key = std::make_tuple(incoming_pdg, positionX, positionY, positionZ);
        auto it = HitCollectionMap.find(key);
	if (it != HitCollectionMap.end()) {
	       it->second.push_back(SimCalHit);
            } else {
                // If not found, create a new vector and move it into the map
               HitCollectionMap[key] = {SimCalHit};
            }
        }
    }

    return StatusCode::SUCCESS;
}

StatusCode FanEcalDigiAlg::MergeHits( const edm4hep::SimCalorimeterHitCollection& m_col, std::vector<edm4hep::SimCalorimeterHit>& m_hits ){

        m_hits.clear(); 
	std::vector<edm4hep::MutableSimCalorimeterHit> m_mergedhit;
	m_mergedhit.clear();

	for(int iter=0; iter<m_col.size(); iter++){
		edm4hep::SimCalorimeterHit m_step = m_col[iter];
		if(!m_step.isAvailable()){ cout<<"ERROR HIT!"<<endl; continue;}
		if(m_step.getEnergy()==0 || m_step.contributions_size()<1) continue;

		edm4hep::CaloHitContribution m_calohitcon = m_step.getContributions(0);
		unsigned long long cellid = m_step.getCellID();
		dd4hep::Position hitpos = m_cellIDConverter->position(cellid);
		edm4hep::Vector3f pos(hitpos.x()*10, hitpos.y()*10, hitpos.z()*10);

		edm4hep::MutableCaloHitContribution conb;
		conb.setEnergy(m_step.getEnergy());
		conb.setStepPosition(m_step.getPosition());
    		conb.setParticle( m_step.getContributions(0).getParticle() );
		conb.setTime(m_step.getContributions(0).getTime());

		edm4hep::MutableSimCalorimeterHit m_hit = find(m_mergedhit, cellid);
		if(m_hit.getCellID()==0){
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

StatusCode FanEcalDigiAlg::SplitCaloHits(const std::map<std::tuple<int, double, double, double>, std::vector<edm4hep::SimCalorimeterHit>>& m_HitColMap, std::vector<std::pair<std::tuple<int, double, double, double>, edm4hep::SimCalorimeterHit>>& m_hits ){
        m_hits.clear();
        std::vector<edm4hep::MutableSimCalorimeterHit> m_mergedhit;
        m_mergedhit.clear();

  for (const auto& entry : m_HitColMap) {
        const std::vector<edm4hep::SimCalorimeterHit>& hitVector = entry.second;
        if (hitVector.empty()) continue;

	edm4hep::MutableSimCalorimeterHit m_hit;
        m_hit.setCellID(hitVector[0].getCellID());

        for (const auto& m_step : hitVector) {
            if (!m_step.isAvailable()) {
                cout << "ERROR HIT!" << endl;
                continue;
            }

            if (m_step.getEnergy() == 0 || m_step.contributions_size() < 1)
                continue;

            edm4hep::CaloHitContribution m_calohitcon = m_step.getContributions(0);
            unsigned long long cellid = m_step.getCellID();
            dd4hep::Position hitpos = m_cellIDConverter->position(cellid);
            edm4hep::Vector3f pos(hitpos.x() * 10, hitpos.y() * 10, hitpos.z() * 10);

            edm4hep::MutableCaloHitContribution conb;
            conb.setEnergy(m_step.getEnergy());
            conb.setStepPosition(m_step.getPosition());
            conb.setParticle(m_step.getContributions(0).getParticle());
            conb.setTime(m_step.getContributions(0).getTime());

            if (m_hit.getCellID() == 0) {
                m_hit.setCellID(cellid);
                m_hit.setPosition(pos);
                m_mergedhit.push_back(m_hit);
            }
            m_hit.addToContributions(conb);
            m_hit.setEnergy(m_hit.getEnergy() + m_step.getEnergy());
        }
    }

  for (const auto& entry : HitCollectionMap) {
        for (const auto& mergedHit : m_mergedhit) {
            if (mergedHit.getCellID() == std::get<0>(entry.first)) {
                edm4hep::SimCalorimeterHit constsimhit = mergedHit;
                m_hits.push_back({entry.first, constsimhit});
                break;  // Move to the next key in HitCollectionMap
             }
        }
    }
  return StatusCode::SUCCESS;
}


edm4hep::MutableSimCalorimeterHit FanEcalDigiAlg::find(std::vector<edm4hep::MutableSimCalorimeterHit>& m_col, unsigned long long& cellid){
   for(int i=0;i<m_col.size();i++){
                edm4hep::MutableSimCalorimeterHit hit=m_col.at(i);
                if(hit.getCellID() == cellid) return hit;
        }
        edm4hep::MutableSimCalorimeterHit hit;
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

StatusCode FanEcalDigiAlg::CLUEClustering(std::vector<CaloBar>& m_barVec, std::vector<CaloCluster>& m_clus, string input, double dc, double deltao, double deltac, double rho){
   
   if(m_barVec.size()==0) return StatusCode::SUCCESS;

   std::vector<float> x;        x.clear();
   std::vector<float> y;        y.clear();
   std::vector<int> layer;      layer.clear();
   std::vector<float> weight;   weight.clear();

   for(int ib=0; ib<m_barVec.size(); ib++ ){
      if(input.compare("even") == 0){
	if((int)m_barVec[ib].getModule()%2 == 0){
	   x.push_back(m_barVec[ib].getModule());
	   y.push_back(m_barVec[ib].getCrystal());
	   layer.push_back(0);
	   weight.push_back(m_barVec[ib].getEnergy());
	}
      }
      if(input.compare("odd") == 0){
	if((int)m_barVec[ib].getModule()%2 != 0){
           x.push_back(m_barVec[ib].getModule());
           y.push_back(m_barVec[ib].getCrystal());
           layer.push_back(0);
           weight.push_back(m_barVec[ib].getEnergy());
        }
      }
   }
   //dc: critical distance used to compute the local density; 
   //deltao: maximum distance for a point to be linked to a nearest higher point
   //deltac: minimum distance for a local high density point to be promoted as a Seed
   //rho: minimum local density for a point to be promoted as a Seed;
   CLUEAlgo clueAlgo(dc,deltao,deltac,rho,-1);  //dc deltao deltac rho 1,5,2,0.2
//   for (int r = 0; r<4; r++){
      clueAlgo.setPoints(x.size(), &x[0],&y[0],&layer[0],&weight[0]);
      clueAlgo.makeClusters();
//   }
   Points resultPoints; resultPoints.clear();
   string outputFileType = "csvFile";
   string outputFileName = "result";
   clueAlgo.verboseResults(resultPoints, outputFileType, outputFileName, -1);
   clueAlgo.verboseToCluster(resultPoints, m_clus); 
   
   return StatusCode::SUCCESS;
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
  m_clus_phi.clear();
  m_clus_phi_start.clear();
  m_clus_Z.clear();
  m_clus_aphi_E.clear();
  m_clus_aphi.clear();
  m_clus_E.clear();
  m_clus_chi2.clear();
  m_clus_alpha.clear();
  m_clus_beta.clear();
  m_clus_Nbars.clear();
 
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
	m_Photon_E_MC.clear();
	
	m_MCP_pdg.clear();
	m_MCP_charge.clear();
	m_MCP_module.clear();
	m_MCP_crystal.clear();
	m_MCP_endpoint_x.clear();
	m_MCP_endpoint_y.clear();
	m_MCP_endpoint_z.clear();

   m_MChit_crystal.clear();
   m_MChit_module.clear();
   m_MChit_E.clear();
   m_MChit_x.clear();
   m_MChit_y.clear();
   m_MChit_z.clear();
   m_MChit_primary_pdg.clear();
   m_MChit_primary_E.clear();
   m_MChit_primary_x.clear();
   m_MChit_primary_y.clear();
   m_MChit_primary_z.clear();
   m_MChit_shower_pdg.clear();

}

#endif
