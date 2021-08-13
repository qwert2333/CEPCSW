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
  m_edmsvc = service<ICRDEcalSvc>("CRDEcalSvc");
  if ( !m_edmsvc )  throw "CRDEcalDigiAlg :Failed to find CRDEcalSvc ...";
	
  rndm.SetSeed(_seed);
  std::cout<<"PandoraPlusPFAlg::initialize"<<std::endl;

  //Output file for code testing. 
  std::string s_outfile = _filename;
  m_wfile = new TFile(s_outfile.c_str(), "recreate");
  t_SimBar = new TTree("SimBarHit", "SimBarHit");
  t_ArborTree = new TTree("ArborTree", "ArborTree");
  t_dataColIter0 = new TTree("dataColIter0", "dataColIter0");
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

  t_ArborTree->Branch("tree_centx", &m_tree_centx); 
  t_ArborTree->Branch("tree_centy", &m_tree_centy); 
  t_ArborTree->Branch("tree_centz", &m_tree_centz); 
  t_ArborTree->Branch("tree_minLayer", &m_tree_minLayer); 
  t_ArborTree->Branch("tree_maxLayer", &m_tree_maxLayer); 
  t_ArborTree->Branch("Nnode", &m_Nnode);
  t_ArborTree->Branch("nodex", &m_nodex);
  t_ArborTree->Branch("nodey", &m_nodey);
  t_ArborTree->Branch("nodez", &m_nodez);
  t_ArborTree->Branch("nodeE", &m_nodeE);
  t_ArborTree->Branch("nodeType", &m_nodeType);

  t_dataColIter0->Branch("Ngoodclus", &m_Iter0_Ngoodclus );
  t_dataColIter0->Branch("Nbadclus", &m_Iter0_Nbadclus );
  t_dataColIter0->Branch("clus_Nly", &m_Iter0_clus_Nly );
  t_dataColIter0->Branch("clus_x", &m_Iter0_clus_x );
  t_dataColIter0->Branch("clus_y", &m_Iter0_clus_y );
  t_dataColIter0->Branch("clus_z", &m_Iter0_clus_z );
  t_dataColIter0->Branch("clus_E", &m_Iter0_clus_E );
  t_dataColIter0->Branch("clus_px", &m_Iter0_clus_px );
  t_dataColIter0->Branch("clus_py", &m_Iter0_clus_py );
  t_dataColIter0->Branch("clus_pz", &m_Iter0_clus_pz );
  t_dataColIter0->Branch("clus_chi2", &m_Iter0_clus_chi2 );
  t_dataColIter0->Branch("clus_aveE", &m_Iter0_clus_aveE );
  t_dataColIter0->Branch("clus_stdDevE", &m_Iter0_clus_stdDevE );
  t_dataColIter0->Branch("clus_start", &m_Iter0_clus_start);
  t_dataColIter0->Branch("clus_end", &m_Iter0_clus_end);
  t_dataColIter0->Branch("clus_alpha", &m_Iter0_clus_alpha);
  t_dataColIter0->Branch("clus_beta", &m_Iter0_clus_beta);
  t_dataColIter0->Branch("clus_showermax", &m_Iter0_clus_showermax);
  t_dataColIter0->Branch("gclus_2dshx", &m_Iter0_gclus_2dshx );
  t_dataColIter0->Branch("gclus_2dshy", &m_Iter0_gclus_2dshy );
  t_dataColIter0->Branch("gclus_2dshz", &m_Iter0_gclus_2dshz );
  t_dataColIter0->Branch("gclus_2dshE", &m_Iter0_gclus_2dshE );

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
  t_recoPFO->Branch("scndM", &m_scndM);

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
  m_pPfoCreator->CreatePFO( m_DataCol ); 


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

  //Save Cluster info
  m_Iter0_Ngoodclus = m_DataCol.Clus3DCol.size();
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_3dclusCol_iter0 = m_DataCol.Clus3DCol;
  for(int icl=0; icl<m_3dclusCol_iter0.size(); icl++){
  ClearIter();
    CRDEcalEDM::CRDCaloHit3DCluster m_clus = m_3dclusCol_iter0[icl];
    m_clus.FitProfile(); 
    m_Iter0_clus_Nly.push_back(m_clus.get2DShowers().size());
    m_Iter0_clus_x.push_back(m_clus.getShowerCenter().x());
    m_Iter0_clus_y.push_back(m_clus.getShowerCenter().y());
    m_Iter0_clus_z.push_back(m_clus.getShowerCenter().z());
    m_Iter0_clus_E.push_back(m_clus.getShowerE());
    m_Iter0_clus_px.push_back(m_clus.getAxis().x());
    m_Iter0_clus_py.push_back(m_clus.getAxis().y());
    m_Iter0_clus_pz.push_back(m_clus.getAxis().z());
    m_Iter0_clus_chi2.push_back(m_clus.getChi2()); 
    m_Iter0_clus_aveE.push_back(m_clus.getAveE());
    m_Iter0_clus_stdDevE.push_back(m_clus.getStdDevE());
    m_Iter0_clus_alpha.push_back(m_clus.getFitAlpha());
    m_Iter0_clus_beta.push_back(m_clus.getFitBeta());

    //std::vector<double> m_widthVec = m_clus.getClusterWidth();
    m_Iter0_clus_showermax.push_back(m_clus.getMaxWidth());
    std::vector<double> m_EnVec = m_clus.getEnInLayer(); 
    int startLayer = 0; 
    for(int i=0; i<m_EnVec.size(); i++){
      if(m_EnVec[i]<0.1) continue; 
      startLayer=i; 
      break; 
    }
    m_Iter0_clus_start.push_back(startLayer);
    m_Iter0_clus_end.push_back(m_clus.getEndDlayer());

    for(int ig=0; ig<m_clus.get2DShowers().size(); ig++){
      m_Iter0_gclus_2dshx.push_back(m_clus.get2DShowers()[ig].getPos().x() );
      m_Iter0_gclus_2dshy.push_back(m_clus.get2DShowers()[ig].getPos().y() );
      m_Iter0_gclus_2dshz.push_back(m_clus.get2DShowers()[ig].getPos().z() );
      m_Iter0_gclus_2dshE.push_back(m_clus.get2DShowers()[ig].getShowerE() );
    }
  t_dataColIter0->Fill();
  }

  //Save second moment in each layer
  ClearRecPFO();
  m_scndM.resize(14);
  std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_2DshowerCol = m_DataCol.Shower2DCol;
  std::map<int, std::vector<CRDEcalEDM::CRDCaloHit2DShower> > m_orderedShower; m_orderedShower.clear(); 
  for(int is=0;is<m_2DshowerCol.size();is++){
    m_orderedShower[m_2DshowerCol[is].getDlayer()].push_back(m_2DshowerCol[is]);
  }
  for(int il=0; il<14; il++){
    std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_showers = m_orderedShower[il];
    if(m_showers.size()==0) { m_scndM[il]=0; continue; }
    if(m_showers.size()==1) { m_scndM[il]=m_showers[0].getHitsWidth(); continue; }

    std::vector<edm4hep::ConstCalorimeterHit> m_hits; m_hits.clear(); 
    double totE = 0;
    for(int is=0; is<m_showers.size(); is++){
      std::vector<edm4hep::ConstCalorimeterHit> m_calohits = m_showers[is].getCaloHits(); 
      m_hits.insert(m_hits.end(), m_calohits.begin(), m_calohits.end() );
      totE += m_showers[is].getHitsE();
    }

    TVector3 cent(0., 0., 0.);
    for(int i=0;i<m_hits.size();i++){
      TVector3 pos(m_hits[i].getPosition().x, m_hits[i].getPosition().y, m_hits[i].getPosition().z);
      cent += pos* (m_hits[i].getEnergy()/totE);
    }

    double width=0;
    for(int i=0;i<m_hits.size();i++){
      TVector3 pos(m_hits[i].getPosition().x, m_hits[i].getPosition().y, m_hits[i].getPosition().z);
      double r2 = (pos-cent).Mag2();
      width += r2*m_hits[i].getEnergy()/totE;
    }
    m_scndM[il] = width; 
  }  

  //Save PFO
  std::vector< CRDEcalEDM::PFObject > m_pfoCol = m_DataCol.PFOCol;
  m_Npfo = m_pfoCol.size(); 
  for(int ipfo=0;ipfo<m_Npfo;ipfo++){
    m_recPFO_px.push_back(m_pfoCol[ipfo].getP4().Px());
    m_recPFO_py.push_back(m_pfoCol[ipfo].getP4().Py());
    m_recPFO_pz.push_back(m_pfoCol[ipfo].getP4().Pz());
    m_recPFO_En.push_back(m_pfoCol[ipfo].getEnergy()); 
    m_recPFO_pid.push_back(m_pfoCol[ipfo].getPdgID());
  }

  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_clus = m_DataCol.Clus3DCol;
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
  t_ArborTree->Write();
  t_dataColIter0->Write();
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

  delete m_wfile, t_SimBar, t_recoPFO, t_ArborTree, t_dataColIter0;

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


void PandoraPlusPFAlg::ClearTree(){
  m_tree_centx = -99.;
  m_tree_centy = -99.;
  m_tree_centz = -99.;
  m_tree_minLayer = -99.;
  m_tree_maxLayer = -99.;
  m_Nnode = -99; 
  m_nodex.clear(); 
  m_nodey.clear(); 
  m_nodez.clear(); 
  m_nodeE.clear(); 
  m_nodeType.clear(); 
}

void PandoraPlusPFAlg::ClearIter(){
  m_Iter0_Ngoodclus=-99;
  m_Iter0_Nbadclus=-99;
  m_Iter0_clus_Nly.clear();
  m_Iter0_clus_x.clear();
  m_Iter0_clus_y.clear();
  m_Iter0_clus_z.clear();
  m_Iter0_clus_E.clear();
  m_Iter0_clus_px.clear();
  m_Iter0_clus_py.clear();
  m_Iter0_clus_pz.clear();
  m_Iter0_clus_chi2.clear(); 
  m_Iter0_clus_aveE.clear();
  m_Iter0_clus_stdDevE.clear();
  m_Iter0_clus_start.clear(); 
  m_Iter0_clus_end.clear(); 
  m_Iter0_clus_alpha.clear();
  m_Iter0_clus_beta.clear();
  m_Iter0_clus_showermax.clear();
  m_Iter0_gclus_2dshx.clear();
  m_Iter0_gclus_2dshy.clear();
  m_Iter0_gclus_2dshz.clear();
  m_Iter0_gclus_2dshE.clear();
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
  m_scndM.clear(); 
}

#endif
