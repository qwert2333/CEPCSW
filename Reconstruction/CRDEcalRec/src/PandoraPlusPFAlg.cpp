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
  t_Layers = new TTree("RecLayers","RecLayers");
  t_Showers = new TTree("RecShowers", "RecShowers");
  t_Clusters = new TTree("RecClusrers", "RecClusrers");
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

  t_Layers->Branch("NshowerX", &m_NshowerX);
  t_Layers->Branch("NshowerY", &m_NshowerY);
  t_Layers->Branch("dlayer", &m_dlayer);
  t_Layers->Branch("module", &m_module);
  t_Layers->Branch("stave", &m_stave);
  t_Layers->Branch("part", &m_part);
  t_Layers->Branch("barShowerX_x", &m_barShowerX_x);
  t_Layers->Branch("barShowerX_y", &m_barShowerX_y);
  t_Layers->Branch("barShowerX_z", &m_barShowerX_z);
  t_Layers->Branch("barShowerX_E", &m_barShowerX_E);
  t_Layers->Branch("barShowerY_x", &m_barShowerY_x);
  t_Layers->Branch("barShowerY_y", &m_barShowerY_y);
  t_Layers->Branch("barShowerY_z", &m_barShowerY_z);
  t_Layers->Branch("barShowerY_E", &m_barShowerY_E);

  t_Showers->Branch("Nshower2D", &m_Nshower2D);
  t_Showers->Branch("shower2D_dlayer", &m_shower2D_dlayer);
  t_Showers->Branch("shower2D_part", &m_shower2D_part);
  t_Showers->Branch("shower2D_stave", &m_shower2D_stave);
  t_Showers->Branch("shower2D_module", &m_shower2D_module);
  t_Showers->Branch("shower2D_type", &m_shower2D_type);
  t_Showers->Branch("shower2D_x", &m_shower2D_x);
  t_Showers->Branch("shower2D_y", &m_shower2D_y);
  t_Showers->Branch("shower2D_z", &m_shower2D_z);
  t_Showers->Branch("shower2D_E", &m_shower2D_E);

  t_Clusters->Branch("Ngoodclus", &m_Ngoodclus);
  t_Clusters->Branch("Nbadclus", &m_Nbadclus);
  t_Clusters->Branch("Nclus", &m_Nclus);
  t_Clusters->Branch("clus_type", &m_clus_type);
  t_Clusters->Branch("clus_Nlayer", &m_clus_Nlayer);
  t_Clusters->Branch("clus_Nshower", &m_clus_Nshower);
  t_Clusters->Branch("clus_x", &m_clus_x);
  t_Clusters->Branch("clus_y", &m_clus_y);
  t_Clusters->Branch("clus_z", &m_clus_z);
  t_Clusters->Branch("clus_E", &m_clus_E);
  t_Clusters->Branch("clus_Lstart", &m_clus_Lstart);
  t_Clusters->Branch("clus_Lend", &m_clus_Lend);
  t_Clusters->Branch("clus_stdDevE", &m_clus_stdDevE);

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

  //Save Layer info
  std::vector<CRDEcalEDM::CRDCaloLayer> m_layerCol = m_DataCol.LayerCol; 
  for(int il=0; il<m_layerCol.size(); il++){
    ClearLayerTree(); 
    m_NshowerX = m_layerCol[il].barShowerXCol.size(); 
    m_NshowerY = m_layerCol[il].barShowerYCol.size();
    m_dlayer = m_layerCol[il].getDlayer(); 
    m_module = m_layerCol[il].getModule(); 
    m_stave = m_layerCol[il].getStave();
    m_part = m_layerCol[il].getPart();
    for(int is=0; is<m_NshowerX; is++){
      m_barShowerX_x.push_back(m_layerCol[il].barShowerXCol[is].getPos().x());
      m_barShowerX_y.push_back(m_layerCol[il].barShowerXCol[is].getPos().y());
      m_barShowerX_z.push_back(m_layerCol[il].barShowerXCol[is].getPos().z());
      m_barShowerX_E.push_back(m_layerCol[il].barShowerXCol[is].getE());
    }
    for(int is=0; is<m_NshowerY; is++){
      m_barShowerY_x.push_back(m_layerCol[il].barShowerYCol[is].getPos().x());
      m_barShowerY_y.push_back(m_layerCol[il].barShowerYCol[is].getPos().y());
      m_barShowerY_z.push_back(m_layerCol[il].barShowerYCol[is].getPos().z());
      m_barShowerY_E.push_back(m_layerCol[il].barShowerYCol[is].getE());
    }
    t_Layers->Fill(); 
  }

  //Save Shower info
  ClearShowerTree();
  m_Nshower2D = m_DataCol.MIPShower2DCol.size() + m_DataCol.EMShower2DCol.size() + m_DataCol.Shower2DCol.size(); 
  for(int is=0; is<m_DataCol.MIPShower2DCol.size(); is++){
    m_shower2D_type.push_back(0); 
    m_shower2D_dlayer.push_back( m_DataCol.MIPShower2DCol[is].getDlayer() ); 
    m_shower2D_part.push_back( m_DataCol.MIPShower2DCol[is].getPart() ); 
    m_shower2D_stave.push_back( m_DataCol.MIPShower2DCol[is].getStave() ); 
    m_shower2D_module.push_back( m_DataCol.MIPShower2DCol[is].getModule() ); 
    m_shower2D_x.push_back( m_DataCol.MIPShower2DCol[is].getPos().x() );
    m_shower2D_y.push_back( m_DataCol.MIPShower2DCol[is].getPos().y() );
    m_shower2D_z.push_back( m_DataCol.MIPShower2DCol[is].getPos().z() );
    m_shower2D_E.push_back( m_DataCol.MIPShower2DCol[is].getShowerE() );
  }
  for(int is=0; is<m_DataCol.EMShower2DCol.size(); is++){
    m_shower2D_type.push_back(1);
    m_shower2D_dlayer.push_back( m_DataCol.EMShower2DCol[is].getDlayer() );
    m_shower2D_part.push_back( m_DataCol.EMShower2DCol[is].getPart() );
    m_shower2D_stave.push_back( m_DataCol.EMShower2DCol[is].getStave() );
    m_shower2D_module.push_back( m_DataCol.EMShower2DCol[is].getModule() );
    m_shower2D_x.push_back( m_DataCol.EMShower2DCol[is].getPos().x() );
    m_shower2D_y.push_back( m_DataCol.EMShower2DCol[is].getPos().y() );
    m_shower2D_z.push_back( m_DataCol.EMShower2DCol[is].getPos().z() );
    m_shower2D_E.push_back( m_DataCol.EMShower2DCol[is].getShowerE() );
  }
  for(int is=0; is<m_DataCol.Shower2DCol.size(); is++){
    m_shower2D_type.push_back(2);
    m_shower2D_dlayer.push_back( m_DataCol.Shower2DCol[is].getDlayer() );
    m_shower2D_part.push_back( m_DataCol.Shower2DCol[is].getPart() );
    m_shower2D_stave.push_back( m_DataCol.Shower2DCol[is].getStave() );
    m_shower2D_module.push_back( m_DataCol.Shower2DCol[is].getModule() );
    m_shower2D_x.push_back( m_DataCol.Shower2DCol[is].getPos().x() );
    m_shower2D_y.push_back( m_DataCol.Shower2DCol[is].getPos().y() );
    m_shower2D_z.push_back( m_DataCol.Shower2DCol[is].getPos().z() );
    m_shower2D_E.push_back( m_DataCol.Shower2DCol[is].getShowerE() );
  }
  t_Showers->Fill();

  //Save Cluster info
  ClearClusterTree(); 
  m_Ngoodclus = m_DataCol.GoodClus3DCol.size(); 
  m_Nbadclus  = m_DataCol.BadClus3DCol.size(); 
  m_Nclus     = m_DataCol.Clus3DCol.size(); 
  for(int ic=0; ic<m_Nclus; ic++){
    m_clus_type.push_back( m_DataCol.Clus3DCol[ic].getType() );
    m_clus_Nlayer.push_back( m_DataCol.Clus3DCol[ic].getEndDlayer()-m_DataCol.Clus3DCol[ic].getBeginningDlayer() );
    m_clus_Nshower.push_back( m_DataCol.Clus3DCol[ic].get2DShowers().size() );
    m_clus_x.push_back( m_DataCol.Clus3DCol[ic].getShowerCenter().x() );
    m_clus_y.push_back( m_DataCol.Clus3DCol[ic].getShowerCenter().y() );
    m_clus_z.push_back( m_DataCol.Clus3DCol[ic].getShowerCenter().z() );
    m_clus_E.push_back( m_DataCol.Clus3DCol[ic].getShowerE() );
    m_clus_Lend.push_back(m_DataCol.Clus3DCol[ic].getEndDlayer());
    m_clus_stdDevE.push_back(m_DataCol.Clus3DCol[ic].getStdDevE());

    float Lstart=0;
    bool f_found = false;
    std::vector<double> m_EnVec = m_DataCol.Clus3DCol[ic].getEnInLayer();
    for(int i=0; i<m_EnVec.size(); i++)
      if(!f_found && m_EnVec[i]>0.1) { Lstart=i; f_found==true; }

    m_clus_Lstart.push_back(Lstart);
  }
  t_Clusters->Fill();

  //Save PFO
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
  t_Layers->Write();
  t_Showers->Write();
  t_Clusters->Write();
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

  delete m_wfile, t_SimBar, t_Layers, t_Showers, t_Clusters, t_recoPFO; 

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

void PandoraPlusPFAlg::ClearLayerTree(){
  m_NshowerX = -99;
  m_NshowerY = -99;
  m_dlayer = -99;
  m_module = -99;
  m_stave = -99;
  m_part = -99;
  m_barShowerX_x.clear(); 
  m_barShowerX_y.clear();
  m_barShowerX_z.clear(); 
  m_barShowerX_E.clear(); 
  m_barShowerY_x.clear();
  m_barShowerY_y.clear();
  m_barShowerY_z.clear();
  m_barShowerY_E.clear();
}

void PandoraPlusPFAlg::ClearShowerTree(){
  m_Nshower2D = -99;
  m_shower2D_dlayer.clear(); 
  m_shower2D_part.clear(); 
  m_shower2D_stave.clear(); 
  m_shower2D_module.clear(); 
  m_shower2D_type.clear(); 
  m_shower2D_x.clear(); 
  m_shower2D_y.clear(); 
  m_shower2D_z.clear(); 
  m_shower2D_E.clear(); 
}

void PandoraPlusPFAlg::ClearClusterTree(){
  m_Ngoodclus = -99;
  m_Nbadclus = -99;
  m_Nclus = -99;
  m_clus_type.clear(); 
  m_clus_Nlayer.clear(); 
  m_clus_Nshower.clear(); 
  m_clus_x.clear(); 
  m_clus_y.clear(); 
  m_clus_z.clear(); 
  m_clus_E.clear(); 
  m_clus_Lstart.clear();
  m_clus_Lend.clear();
  m_clus_stdDevE.clear(); 
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
