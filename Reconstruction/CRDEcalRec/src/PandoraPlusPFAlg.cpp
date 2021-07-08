#ifndef PANDORAPLUS_ALG_C
#define PANDORAPLUS_ALG_C

#include "PandoraPlusPFAlg.h"


#define C 299.79  // unit: mm/ns
#define PI 3.141592653
using namespace std;
using namespace dd4hep;


DECLARE_COMPONENT( PandoraPlusPFAlg )

PandoraPlusPFAlg::PandoraPlusPFAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{
   declareProperty("MCParticle",      r_MCParticleCol, "Handle of the Input MCParticle collection");
   declareProperty("MarlinTrkTracks", r_MarlinTrkCol,  "Handle of the Input Track collection");
 
   
}

StatusCode PandoraPlusPFAlg::initialize()
{
  //Initialize Settings
  m_pMCParticleCreatorSettings = new MCParticleCreator::Settings;
  m_pTrackCreatorSettings      = new TrackCreator::Settings;
  m_pVertexCreatorSettings     = new VertexCreator::Settings;
  m_EcalHitsCreatorSettings    = new EcalHitsCreator::Settings;
  m_pHcalHitsCreatorSettings   = new HcalHitsCreator::Settings;
  m_pEcalClusterRecSettings    = new EcalClusterReconstruction::Settings;
  m_pPfoCreatorSettings        = new PFOCreator::Settings; 

  //Initialize Creators
  m_pMCParticleCreator = new MCParticleCreator( *m_pMCParticleCreatorSettings );
  m_pTrackCreator      = new TrackCreator( *m_pTrackCreatorSettings );
  m_pVertexCreator     = new VertexCreator( *m_pVertexCreatorSettings );
  m_pEcalHitsCreator   = new EcalHitsCreator( *m_EcalHitsCreatorSettings );
  m_pHcalHitsCreator   = new HcalHitsCreator( *m_pHcalHitsCreatorSettings );
  m_pPfoCreator        = new PFOCreator( *m_pPfoCreatorSettings );

  //Initialize Algorithm
  m_pEcalClusterRec = new EcalClusterReconstruction();
  m_pEcalClusterRec->Initialize();


  //Initialize services
  m_edmsvc = service<ICRDEcalEDMSvc>("CRDEcalEDMSvc");
  if ( !m_edmsvc )  throw "CRDEcalDigiAlg :Failed to find CRDEcalEDMSvc ...";
	
  rndm.SetSeed(_seed);
  std::cout<<"PandoraPlusPFAlg::initialize"<<std::endl;

  //Output file for code testing. 
  std::string s_outfile = _filename;
  m_wfile = new TFile(s_outfile.c_str(), "recreate");
  t_SimBar = new TTree("SimBarHit", "SimBarHit");
  t_ArborTree0 = new TTree("ArborTreeX", "ArborTreeX");
  t_ArborTree1 = new TTree("ArborTreeY", "ArborTreeY");
  t_recoPFO = new TTree("RecPFO",  "RecPFO");

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
  t_SimBar->Branch("simBar_stave", &m_simBar_stave);
  t_SimBar->Branch("simBar_slayer", &m_simBar_slayer);

  t_ArborTree0->Branch("tree0_centx", &m_tree0_centx); 
  t_ArborTree0->Branch("tree0_centy", &m_tree0_centy); 
  t_ArborTree0->Branch("tree0_centz", &m_tree0_centz); 
  t_ArborTree0->Branch("tree0_minLayer", &m_tree0_minLayer); 
  t_ArborTree0->Branch("tree0_maxLayer", &m_tree0_maxLayer); 
  t_ArborTree0->Branch("Nnode0", &m_Nnode0);
  t_ArborTree0->Branch("node0x", &m_node0x);
  t_ArborTree0->Branch("node0y", &m_node0y);
  t_ArborTree0->Branch("node0z", &m_node0z);
  t_ArborTree0->Branch("node0E", &m_node0E);
  t_ArborTree0->Branch("node0Type", &m_node0Type);

  t_ArborTree1->Branch("tree1_centx", &m_tree1_centx);
  t_ArborTree1->Branch("tree1_centy", &m_tree1_centy);
  t_ArborTree1->Branch("tree1_centz", &m_tree1_centz);
  t_ArborTree1->Branch("tree1_minLayer", &m_tree1_minLayer);
  t_ArborTree1->Branch("tree1_maxLayer", &m_tree1_maxLayer);
  t_ArborTree1->Branch("Nnode1", &m_Nnode1);
  t_ArborTree1->Branch("node1x", &m_node1x);
  t_ArborTree1->Branch("node1y", &m_node1y);
  t_ArborTree1->Branch("node1z", &m_node1z);
  t_ArborTree1->Branch("node1E", &m_node1E);
  t_ArborTree1->Branch("node1Type", &m_node1Type);

  t_recoPFO->Branch("Npfo", &m_Npfo);
  t_recoPFO->Branch("recPFO_px", &m_recPFO_px);
  t_recoPFO->Branch("recPFO_py", &m_recPFO_py);
  t_recoPFO->Branch("recPFO_pz", &m_recPFO_pz);
  t_recoPFO->Branch("recPFO_En", &m_recPFO_En);
  t_recoPFO->Branch("recPFO_pid", &m_recPFO_pid);
  t_recoPFO->Branch("N3dclus", &m_N3dclus);
  t_recoPFO->Branch("N2dshInClus", &m_N2dshInClus);
  t_recoPFO->Branch("Clus_x", &m_Clus_x);
  t_recoPFO->Branch("Clus_y", &m_Clus_y);
  t_recoPFO->Branch("Clus_z", &m_Clus_z);
  t_recoPFO->Branch("Clus_E", &m_Clus_E);
  t_recoPFO->Branch("Clus_type", &m_Clus_type);
  t_recoPFO->Branch("mcPdgid",     &m_mcPdgid);
  t_recoPFO->Branch("mcStatus",    &m_mcStatus);
  t_recoPFO->Branch("mcNdaughter", &m_mcNdaughter);
  t_recoPFO->Branch("mcNparent",   &m_mcNparent);
  t_recoPFO->Branch("mcPx", &m_mcPx);
  t_recoPFO->Branch("mcPy", &m_mcPy);
  t_recoPFO->Branch("mcPz", &m_mcPz);
  t_recoPFO->Branch("mcEn", &m_mcEn);
  t_recoPFO->Branch("Nmc",  &m_Nmc);

  return GaudiAlgorithm::initialize();
}

StatusCode PandoraPlusPFAlg::execute()
{
  if(_nEvt==0) std::cout<<"PandoraPlusPFAlg::execute Start"<<std::endl;
  std::cout<<"Processing event: "<<_nEvt<<std::endl;
  if(_nEvt<_Nskip){ _nEvt++;  return GaudiAlgorithm::initialize(); }


  //InitializeForNewEvent(); 
  m_DataCol.Clear();
  //m_edmsvc->ClearSystem();

  //Get dataCol from service
  const edm4hep::MCParticleCollection* const_MCPCol = r_MCParticleCol.get();
  //const edm4hep::TrackCollection*      const_TrkCol = r_MarlinTrkCol.get(); 
  if( const_MCPCol!=NULL ) m_pMCParticleCreator->GetMCParticle( m_DataCol, *const_MCPCol);
  //if( const_TrkCol!=NULL ) m_pTrackCreator->     GetTracks( m_DataCol, *const_TrkCol );
  m_pTrackCreator->     GetTracks( m_DataCol );
  m_pVertexCreator->    GetVertex( m_DataCol );
  m_pEcalHitsCreator->  GetEcalBars( m_DataCol, *m_edmsvc); 
  m_pHcalHitsCreator->  GetHcalHits( m_DataCol );


  //Link MCParticle-Track and MCParticle-Hit. 
  m_pMCParticleCreator->CreateTrackMCParticleRelation();
  m_pMCParticleCreator->CreateEcalBarMCParticleRelation();
  m_pMCParticleCreator->CreateHcalHitsMCParticleRelation();

  //Match track with Ecal super-cell.
  m_pTrackCreator->MatchTrkEcalRelation( m_DataCol);


  //Perform PFA algorithm
  m_pEcalClusterRec->Initialize();
  m_pEcalClusterRec->RunAlgorithm( *m_pEcalClusterRecSettings, m_DataCol );
  m_pEcalClusterRec->ClearAlgorithm();


  //std::cout<<"Reconstruction is done. Create PFO"<<std::endl;
  //m_pPfoCreator->CreatePFO( m_DataCol ); 


  //PFA algorithm is end!
  //-----------------------------------------------------
  //Followings are for code check. 


  //Print out information for algorithm check
//  m_DataCol.PrintLayer();
//  m_DataCol.PrintShower();
//  m_DataCol.Print3DClus();


  //Save Raw bars information
  ClearBar();
  std::vector<CRDEcalEDM::CRDCaloBlock>  tmp_stavevec = m_DataCol.BlockVec;
  for(int ibl=0;ibl<tmp_stavevec.size();ibl++){
    CRDEcalEDM::CRDCaloBlock tmp_barcol = tmp_stavevec[ibl];
    for(int ibar=0;ibar<tmp_barcol.getBarXCol().size();ibar++){
      CRDEcalEDM::CRDCaloBar hitbar = tmp_barcol.getBarXCol()[ibar];
      m_simBar_x.push_back(hitbar.getPosition().x());
      m_simBar_y.push_back(hitbar.getPosition().y());
      m_simBar_z.push_back(hitbar.getPosition().z());
      m_simBar_Q1.push_back(hitbar.getQ1());
      m_simBar_Q2.push_back(hitbar.getQ2());
      m_simBar_T1.push_back(hitbar.getT1());
      m_simBar_T2.push_back(hitbar.getT2());
      m_simBar_module.push_back(hitbar.getModule());
      m_simBar_dlayer.push_back(hitbar.getDlayer());
      m_simBar_part.push_back(hitbar.getPart());
      m_simBar_stave.push_back(hitbar.getStave());
      m_simBar_slayer.push_back(hitbar.getSlayer());
    }
    for(int ibar=0;ibar<tmp_barcol.getBarYCol().size();ibar++){
      CRDEcalEDM::CRDCaloBar hitbar = tmp_barcol.getBarYCol()[ibar];
      m_simBar_x.push_back(hitbar.getPosition().x());
      m_simBar_y.push_back(hitbar.getPosition().y());
      m_simBar_z.push_back(hitbar.getPosition().z());
      m_simBar_Q1.push_back(hitbar.getQ1());
      m_simBar_Q2.push_back(hitbar.getQ2());
      m_simBar_T1.push_back(hitbar.getT1());
      m_simBar_T2.push_back(hitbar.getT2());
      m_simBar_module.push_back(hitbar.getModule());
      m_simBar_dlayer.push_back(hitbar.getDlayer());
      m_simBar_part.push_back(hitbar.getPart());
      m_simBar_stave.push_back(hitbar.getStave());
      m_simBar_slayer.push_back(hitbar.getSlayer());
    }
  }
  t_SimBar->Fill();

  //Save Arbor Tree Info
  std::vector<CRDEcalEDM::CRDArborTree> m_treeColX = m_DataCol.ArborTreeColX;
  for(int it=0; it<m_treeColX.size(); it++){
    ClearTree0();

    m_tree0_centx = m_treeColX[it].GetBarycenter().x(); 
    m_tree0_centy = m_treeColX[it].GetBarycenter().y(); 
    m_tree0_centz = m_treeColX[it].GetBarycenter().z();
    m_tree0_minLayer = (float)m_treeColX[it].GetMinDlayer();  
    m_tree0_maxLayer = (float)m_treeColX[it].GetMaxDlayer();  
    m_Nnode0 = m_treeColX[it].GetNodes().size(); 
    for(int in=0; in<m_Nnode0; in++){
      m_node0x.push_back( m_treeColX[it].GetNodes()[in]->GetPosition().x() );
      m_node0y.push_back( m_treeColX[it].GetNodes()[in]->GetPosition().y() );
      m_node0z.push_back( m_treeColX[it].GetNodes()[in]->GetPosition().z() );
      m_node0E.push_back( m_treeColX[it].GetNodes()[in]->GetEnergy() );
      m_node0Type.push_back( m_treeColX[it].GetNodes()[in]->GetType() );
    }    
    t_ArborTree0->Fill();
  }

  std::vector<CRDEcalEDM::CRDArborTree> m_treeColY = m_DataCol.ArborTreeColY;
  for(int it=0; it<m_treeColY.size(); it++){
    ClearTree1();

    m_tree1_centx = m_treeColY[it].GetBarycenter().x();
    m_tree1_centy = m_treeColY[it].GetBarycenter().y();
    m_tree1_centz = m_treeColY[it].GetBarycenter().z();
    m_tree1_minLayer = (float)m_treeColY[it].GetMinDlayer();
    m_tree1_maxLayer = (float)m_treeColY[it].GetMaxDlayer();
    m_Nnode1 = m_treeColY[it].GetNodes().size();
    for(int in=0; in<m_Nnode1; in++){
      m_node1x.push_back( m_treeColY[it].GetNodes()[in]->GetPosition().x() );
      m_node1y.push_back( m_treeColY[it].GetNodes()[in]->GetPosition().y() );
      m_node1z.push_back( m_treeColY[it].GetNodes()[in]->GetPosition().z() );
      m_node1E.push_back( m_treeColY[it].GetNodes()[in]->GetEnergy() );
      m_node1Type.push_back( m_treeColY[it].GetNodes()[in]->GetType() );
    }
    t_ArborTree1->Fill();
  }

  
  //Save PFO information
  ClearRecPFO();
  std::vector< CRDEcalEDM::PFObject > m_pfoCol = m_DataCol.PFOCol;
  m_Npfo = m_pfoCol.size(); 
  for(int ipfo=0;ipfo<m_Npfo;ipfo++){
    m_recPFO_px.push_back(m_pfoCol[ipfo].getP4().Px());
    m_recPFO_py.push_back(m_pfoCol[ipfo].getP4().Py());
    m_recPFO_pz.push_back(m_pfoCol[ipfo].getP4().Pz());
    m_recPFO_En.push_back(m_pfoCol[ipfo].getEnergy()); 
    m_recPFO_pid.push_back(m_pfoCol[ipfo].getPdgID());
  }

  std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_clus = m_DataCol.Clus3DCol;
  m_N3dclus = m_clus.size();
  for(int i=0;i<m_N3dclus;i++){ 
    m_N2dshInClus.push_back(m_clus[i].get2DShowers().size());
    m_Clus_x.push_back(m_clus[i].getShowerCenter().x());
    m_Clus_y.push_back(m_clus[i].getShowerCenter().y());
    m_Clus_z.push_back(m_clus[i].getShowerCenter().z());
    m_Clus_E.push_back(m_clus[i].getShowerE());
    m_Clus_type.push_back(m_clus[i].getType());
  }


  //Save MC information
  std::vector<edm4hep::MCParticle> m_MCPCol = m_DataCol.MCParticleCol; 
  m_Nmc = m_MCPCol.size(); 
  for(int imc=0; imc<m_MCPCol.size(); imc++){
    m_mcPdgid.push_back( m_MCPCol[imc].getPDG() );
    m_mcNdaughter.push_back( m_MCPCol[imc].daughters_size() );
    m_mcNparent.push_back( m_MCPCol[imc].parents_size() );
    m_mcStatus.push_back( m_MCPCol[imc].getGeneratorStatus() );
    m_mcPx.push_back( m_MCPCol[imc].getMomentum()[0] );
    m_mcPy.push_back( m_MCPCol[imc].getMomentum()[1] );
    m_mcPz.push_back( m_MCPCol[imc].getMomentum()[2] );
    m_mcEn.push_back( m_MCPCol[imc].getEnergy() );
  }
  t_recoPFO->Fill();

  //Reset
  m_edmsvc->ClearSystem();

  std::cout<<"Event: "<<_nEvt<<" is done"<<std::endl;
  _nEvt ++ ;
  return StatusCode::SUCCESS;
}

StatusCode PandoraPlusPFAlg::finalize()
{

  m_wfile->cd();
  t_SimBar->Write();
  t_ArborTree0->Write();
  t_ArborTree1->Write();
  t_recoPFO->Write();
  m_wfile->Close();

  delete m_pMCParticleCreator;
  delete m_pTrackCreator;
  delete m_pVertexCreator;
  delete m_pEcalHitsCreator;
  delete m_pHcalHitsCreator;
  delete m_pEcalClusterRec;
  delete m_pPfoCreator;

  delete m_pMCParticleCreatorSettings;
  delete m_pTrackCreatorSettings;
  delete m_pVertexCreatorSettings;
  delete m_EcalHitsCreatorSettings;
  delete m_pHcalHitsCreatorSettings;
  delete m_pPfoCreatorSettings;
  delete m_pEcalClusterRecSettings;

  delete m_wfile, t_SimBar, t_recoPFO, t_ArborTree0, t_ArborTree1; 

  info() << "Processed " << _nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}


void PandoraPlusPFAlg::ClearBar(){
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
  m_simBar_stave.clear();
  m_simBar_slayer.clear();
}


void PandoraPlusPFAlg::ClearTree0(){
  m_tree0_centx = -99.;
  m_tree0_centy = -99.;
  m_tree0_centz = -99.;
  m_tree0_minLayer = -99.;
  m_tree0_maxLayer = -99.;
  m_Nnode0 = -99; 
  m_node0x.clear(); 
  m_node0y.clear(); 
  m_node0z.clear(); 
  m_node0E.clear(); 
  m_node0Type.clear(); 
}

void PandoraPlusPFAlg::ClearTree1(){
  m_tree1_centx = -99.;
  m_tree1_centy = -99.;
  m_tree1_centz = -99.;
  m_tree1_minLayer = -99.;
  m_tree1_maxLayer = -99.;
  m_Nnode1 = -99;
  m_node1x.clear();
  m_node1y.clear();
  m_node1z.clear();
  m_node1E.clear();
  m_node1Type.clear();
}

void PandoraPlusPFAlg::ClearRecPFO(){
  m_recPFO_px.clear(); 
  m_recPFO_py.clear(); 
  m_recPFO_pz.clear(); 
  m_recPFO_En.clear(); 
  m_recPFO_pid.clear();
  m_Clus_x.clear();
  m_Clus_y.clear();
  m_Clus_z.clear();
  m_Clus_E.clear();
  m_Clus_type.clear();
  m_mcPdgid.clear(); 
  m_mcStatus.clear();
  m_mcNdaughter.clear(); 
  m_mcNparent.clear();
  m_N2dshInClus.clear();
  m_mcPx.clear();
  m_mcPy.clear();
  m_mcPz.clear();
  m_mcEn.clear();
  m_Npfo=-99;
  m_Nmc=-99;
  m_N3dclus=-99;
}

#endif
