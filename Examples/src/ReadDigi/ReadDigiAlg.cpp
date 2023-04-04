#ifndef READ_DIGI_C
#define READ_DIGI_C

#include "ReadDigiAlg.h"

DECLARE_COMPONENT(ReadDigiAlg)

ReadDigiAlg::ReadDigiAlg(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc), 
     _nEvt(0)
{
  //declareProperty("MCParticle",     m_MCParticleCol, "MCParticle collection (input)");
  //declareProperty("EcalBarrelCollectionMerged", m_ECalBarrelHitCol, "ECal Barrel");
  //declareProperty("HcalBarrelCollection", m_HCalBarrelHitCol, "HCal Barrel");
  //declareProperty("HcalEndcapsCollection", m_HCalEndcapHitCol, "HCal Endcap");
  //declareProperty("MCRecoCaloAssociationCollection", m_MCRecoCaloAssociationCol, "MC-Reco association");

}

StatusCode ReadDigiAlg::initialize()
{
  debug() << "begin initialize ReadDigiAlg" << endmsg;
  std::string s_outfile = _filename;
  m_wfile = new TFile(s_outfile.c_str(), "recreate");
  m_wtree = new TTree("DigiHits", "DigiHits");

  m_wtree->Branch("N_MCP", &N_MCP);
  m_wtree->Branch("MCP_px", &MCP_px);
  m_wtree->Branch("MCP_py", &MCP_py);
  m_wtree->Branch("MCP_pz", &MCP_pz);
  m_wtree->Branch("MCP_E", &MCP_E);
  m_wtree->Branch("MCP_pdgid", &MCP_pdgid);
  m_wtree->Branch("MCP_gStatus", &MCP_gStatus);

  m_wtree->Branch("Nhit_EcalB", &Nhit_EcalB);
  m_wtree->Branch("EcalBHit_MCtag", &EcalBHit_MCtag);
  m_wtree->Branch("EcalBHit_MCpid", &EcalBHit_MCpid);
  m_wtree->Branch("EcalBHit_x", &EcalBHit_x);
  m_wtree->Branch("EcalBHit_y", &EcalBHit_y);
  m_wtree->Branch("EcalBHit_z", &EcalBHit_z);
  m_wtree->Branch("EcalBHit_E", &EcalBHit_E);
  m_wtree->Branch("EcalBHit_T", &EcalBHit_T);

  m_wtree->Branch("Nhit_EcalE", &Nhit_EcalE);
  m_wtree->Branch("EcalEHit_MCtag", &EcalEHit_MCtag);
  m_wtree->Branch("EcalEHit_MCpid", &EcalEHit_MCpid);
  m_wtree->Branch("EcalEHit_x", &EcalEHit_x);
  m_wtree->Branch("EcalEHit_y", &EcalEHit_y);
  m_wtree->Branch("EcalEHit_z", &EcalEHit_z);
  m_wtree->Branch("EcalEHit_E", &EcalEHit_E);
  m_wtree->Branch("EcalEHit_T", &EcalEHit_T);

  m_wtree->Branch("Nhit_HcalB", &Nhit_HcalB);
  m_wtree->Branch("HcalBHit_MCtag", &HcalBHit_MCtag);
  m_wtree->Branch("HcalBHit_MCpid", &HcalBHit_MCpid);
  m_wtree->Branch("HcalBHit_x", &HcalBHit_x);
  m_wtree->Branch("HcalBHit_y", &HcalBHit_y);
  m_wtree->Branch("HcalBHit_z", &HcalBHit_z);
  m_wtree->Branch("HcalBHit_E", &HcalBHit_E);
  m_wtree->Branch("HcalBHit_T", &HcalBHit_T);

  m_wtree->Branch("Nhit_HcalE", &Nhit_HcalE);
  m_wtree->Branch("HcalEHit_MCtag", &HcalEHit_MCtag);
  m_wtree->Branch("HcalEHit_MCpid", &HcalEHit_MCpid);
  m_wtree->Branch("HcalEHit_x", &HcalEHit_x);
  m_wtree->Branch("HcalEHit_y", &HcalEHit_y);
  m_wtree->Branch("HcalEHit_z", &HcalEHit_z);
  m_wtree->Branch("HcalEHit_E", &HcalEHit_E);
  m_wtree->Branch("HcalEHit_T", &HcalEHit_T);


  return GaudiAlgorithm::initialize();
}

StatusCode ReadDigiAlg::execute()
{
  if(_nEvt==0) std::cout<<"ReadDigi::execute Start"<<std::endl;
  debug() << "Processing event: "<<_nEvt<< endmsg; 
  Clear(); 

  const edm4hep::MCRecoCaloAssociationCollection* const_CaloSimDigiRelCol = m_MCRecoCaloAssociationCol.get();
  const edm4hep::MCParticleCollection* const_MCPCol = m_MCParticleCol.get();
  N_MCP = const_MCPCol->size(); 
  for(int i=0; i<N_MCP; i++){
    edm4hep::MCParticle m_MCp = const_MCPCol->at(i);
    MCP_px.push_back(m_MCp.getMomentum().x);
    MCP_py.push_back(m_MCp.getMomentum().y);
    MCP_pz.push_back(m_MCp.getMomentum().z);
    MCP_E.push_back(m_MCp.getEnergy());
    MCP_pdgid.push_back(m_MCp.getPDG());
    MCP_gStatus.push_back(m_MCp.getGeneratorStatus());
  }


  const edm4hep::CalorimeterHitCollection* ECalBHitCol = m_ECalBarrelHitCol.get();
  Nhit_EcalB = ECalBHitCol->size();
  for(int i=0; i<ECalBHitCol->size(); i++){
    edm4hep::CalorimeterHit CaloHit = ECalBHitCol->at(i);
    EcalBHit_x.push_back(CaloHit.getPosition().x);
    EcalBHit_y.push_back(CaloHit.getPosition().y);
    EcalBHit_z.push_back(CaloHit.getPosition().z);
    EcalBHit_E.push_back(CaloHit.getEnergy());

    double aveT = 0;
    int totConb = 0; 
    MCParticleToEnergyWeightMap MCPEnMap; MCPEnMap.clear();
    for(int irel=0; irel<const_CaloSimDigiRelCol->size(); irel++){
      if(const_CaloSimDigiRelCol->at(irel).getRec().id()!=CaloHit.id()) continue; 
      auto pSimHit = const_CaloSimDigiRelCol->at(irel).getSim(); 
      
      for(int iCont=0; iCont < pSimHit.contributions_size(); ++iCont){
        auto conb = pSimHit.getContributions(iCont);
        if( !conb.isAvailable() ) continue;
        if(conb.getEnergy() == 0) continue;
   
        auto mcp = conb.getParticle();
        MCPEnMap[mcp] += conb.getEnergy();
        float step_time = conb.getTime();
        aveT += step_time;
        totConb ++;
      }
    }

    float EnMax = -99;
    int mcPid = 0; 
    int mcTag = -1; 
    for(auto iter : MCPEnMap){
      if(iter.second>EnMax){ 
        mcPid = iter.first.getPDG();
        mcTag = iter.first.id();
      }
    }

    EcalBHit_MCtag.push_back(mcTag);
    EcalBHit_MCpid.push_back(mcPid);
    EcalBHit_T.push_back(aveT/totConb);
  }


  const edm4hep::CalorimeterHitCollection* HCalBHitCol = m_HCalBarrelHitCol.get();
  Nhit_HcalB = HCalBHitCol->size();
  for(int i=0; i<HCalBHitCol->size(); i++){
    edm4hep::CalorimeterHit CaloHit = HCalBHitCol->at(i);
    HcalBHit_x.push_back(CaloHit.getPosition().x);
    HcalBHit_y.push_back(CaloHit.getPosition().y);
    HcalBHit_z.push_back(CaloHit.getPosition().z);
    HcalBHit_E.push_back(CaloHit.getEnergy());

    double aveT = 0;
    int totConb = 0;
    MCParticleToEnergyWeightMap MCPEnMap; MCPEnMap.clear();
    for(int irel=0; irel<const_CaloSimDigiRelCol->size(); irel++){
      if(const_CaloSimDigiRelCol->at(irel).getRec().id()!=CaloHit.id()) continue;
      auto pSimHit = const_CaloSimDigiRelCol->at(irel).getSim();

      for(int iCont=0; iCont < pSimHit.contributions_size(); ++iCont){
        auto conb = pSimHit.getContributions(iCont);
        if( !conb.isAvailable() ) continue;
        if(conb.getEnergy() == 0) continue;

        auto mcp = conb.getParticle();
        MCPEnMap[mcp] += conb.getEnergy();
        float step_time = conb.getTime();
        aveT += step_time;
        totConb ++;
      }
    }

    float EnMax = -99;
    int mcPid = 0;
    int mcTag = -1;
    for(auto iter : MCPEnMap){
      if(iter.second>EnMax){
        mcPid = iter.first.getPDG();
        mcTag = iter.first.id();
      }
    }

    HcalBHit_MCtag.push_back(mcTag);
    HcalBHit_MCpid.push_back(mcPid);
    HcalBHit_T.push_back(aveT/totConb);

  }

  const edm4hep::CalorimeterHitCollection* HCalEHitCol = m_HCalEndcapHitCol.get();
  Nhit_HcalE = HCalEHitCol->size();
  for(int i=0; i<HCalEHitCol->size(); i++){
    edm4hep::CalorimeterHit CaloHit = HCalEHitCol->at(i);
    HcalEHit_x.push_back(CaloHit.getPosition().x);
    HcalEHit_y.push_back(CaloHit.getPosition().y);
    HcalEHit_z.push_back(CaloHit.getPosition().z);
    HcalEHit_E.push_back(CaloHit.getEnergy());

    double aveT = 0;
    int totConb = 0;
    MCParticleToEnergyWeightMap MCPEnMap; MCPEnMap.clear();
    for(int irel=0; irel<const_CaloSimDigiRelCol->size(); irel++){
      if(const_CaloSimDigiRelCol->at(irel).getRec().id()!=CaloHit.id()) continue;
      auto pSimHit = const_CaloSimDigiRelCol->at(irel).getSim();

      for(int iCont=0; iCont < pSimHit.contributions_size(); ++iCont){
        auto conb = pSimHit.getContributions(iCont);
        if( !conb.isAvailable() ) continue;
        if(conb.getEnergy() == 0) continue;

        auto mcp = conb.getParticle();
        MCPEnMap[mcp] += conb.getEnergy();
        float step_time = conb.getTime();
        aveT += step_time;
        totConb ++;
      }
    }

    float EnMax = -99;
    int mcPid = 0;
    int mcTag = -1;
    for(auto iter : MCPEnMap){
      if(iter.second>EnMax){
        mcPid = iter.first.getPDG();
        mcTag = iter.first.id();
      }
    }

    HcalEHit_MCtag.push_back(mcTag);
    HcalEHit_MCpid.push_back(mcPid);
    HcalEHit_T.push_back(aveT/totConb);
  }  


  m_wtree->Fill(); 
  _nEvt++; 
  return StatusCode::SUCCESS;
}

StatusCode ReadDigiAlg::finalize()
{
  debug() << "begin finalize ReadDigiAlg" << endmsg;
  m_wfile->cd();
  m_wtree->Write();
  m_wfile->Close(); 

  return GaudiAlgorithm::finalize();
}

StatusCode ReadDigiAlg::Clear()
{
  N_MCP = -99;
  MCP_px.clear();
  MCP_py.clear();
  MCP_pz.clear();
  MCP_E.clear();
  MCP_pdgid.clear();
  MCP_gStatus.clear(); 

  Nhit_EcalB = -99;
  EcalBHit_MCtag.clear();
  EcalBHit_MCpid.clear();
  EcalBHit_x.clear();
  EcalBHit_y.clear();
  EcalBHit_z.clear();
  EcalBHit_E.clear();
  EcalBHit_T.clear();

  Nhit_EcalE = -99;
  EcalEHit_MCtag.clear();
  EcalEHit_MCpid.clear();
  EcalEHit_x.clear();
  EcalEHit_y.clear();
  EcalEHit_z.clear();
  EcalEHit_E.clear();
  EcalEHit_T.clear();

  Nhit_HcalB = -99;
  HcalBHit_MCtag.clear();
  HcalBHit_MCpid.clear();
  HcalBHit_x.clear();
  HcalBHit_y.clear();
  HcalBHit_z.clear();
  HcalBHit_E.clear();
  HcalBHit_T.clear();

  Nhit_HcalE = -99;
  HcalEHit_MCtag.clear();
  HcalEHit_MCpid.clear();
  HcalEHit_x.clear();
  HcalEHit_y.clear();
  HcalEHit_z.clear();
  HcalEHit_E.clear();
  HcalEHit_T.clear();

  return StatusCode::SUCCESS;
}

#endif
