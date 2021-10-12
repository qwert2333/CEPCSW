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
  t_dataColIter0 = new TTree("dataColIter0", "dataColIter0");
  t_dataColIter1 = new TTree("dataColIter1", "dataColIter1");
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

  t_dataColIter0->Branch("barshower0X", &m_barshower0X );
  t_dataColIter0->Branch("barshower0Y", &m_barshower0Y );
  t_dataColIter0->Branch("barshower0Z", &m_barshower0Z );
  t_dataColIter0->Branch("barshower0E", &m_barshower0E );
  t_dataColIter0->Branch("barshower0_layer", &m_barshower0_layer );
  t_dataColIter0->Branch("barshower1X", &m_barshower1X );
  t_dataColIter0->Branch("barshower1Y", &m_barshower1Y );
  t_dataColIter0->Branch("barshower1Z", &m_barshower1Z );
  t_dataColIter0->Branch("barshower1E", &m_barshower1E );
  t_dataColIter0->Branch("barshower1_layer", &m_barshower1_layer );
/*  t_dataColIter0->Branch("trkX_x", &m_trkX_x);
  t_dataColIter0->Branch("trkX_y", &m_trkX_y);
  t_dataColIter0->Branch("trkX_z", &m_trkX_z);
  t_dataColIter0->Branch("trkX_px", &m_trkX_px);
  t_dataColIter0->Branch("trkX_py", &m_trkX_py);
  t_dataColIter0->Branch("trkX_pz", &m_trkX_pz);
  t_dataColIter0->Branch("trkX_Nsh", &m_trkX_Nsh);
  t_dataColIter0->Branch("trkY_x", &m_trkY_x);
  t_dataColIter0->Branch("trkY_y", &m_trkY_y);
  t_dataColIter0->Branch("trkY_z", &m_trkY_z);
  t_dataColIter0->Branch("trkY_px", &m_trkY_px);
  t_dataColIter0->Branch("trkY_py", &m_trkY_py);
  t_dataColIter0->Branch("trkY_pz", &m_trkY_pz);
  t_dataColIter0->Branch("trkY_Nsh", &m_trkY_Nsh);
*/

  t_dataColIter1->Branch("barshower0X", &m_barshower0X_iter1);
  t_dataColIter1->Branch("barshower0Y", &m_barshower0Y_iter1);
  t_dataColIter1->Branch("barshower0Z", &m_barshower0Z_iter1);
  t_dataColIter1->Branch("barshower0E", &m_barshower0E_iter1);
  t_dataColIter1->Branch("barshower1X", &m_barshower1X_iter1);
  t_dataColIter1->Branch("barshower1Y", &m_barshower1Y_iter1);
  t_dataColIter1->Branch("barshower1Z", &m_barshower1Z_iter1);
  t_dataColIter1->Branch("barshower1E", &m_barshower1E_iter1);
  t_dataColIter1->Branch("barshower0_layer", &m_barshower0_layer_iter1 );
  t_dataColIter1->Branch("barshower1_layer", &m_barshower1_layer_iter1 );
  t_dataColIter1->Branch("gclus_2dshx", &m_Iter1_gclus_2dshx );
  t_dataColIter1->Branch("gclus_2dshy", &m_Iter1_gclus_2dshy );
  t_dataColIter1->Branch("gclus_2dshz", &m_Iter1_gclus_2dshz );
  t_dataColIter1->Branch("gclus_2dshE", &m_Iter1_gclus_2dshE );


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

  //Save Layer info
  ClearIter0();
  std::vector< CRDEcalEDM::CRDCaloLayer > m_layerCol = m_DataCol.LayerCol_tmp;
  for(int ily=0; ily<m_layerCol.size(); ily++){
    CRDEcalEDM::CRDCaloLayer m_layer = m_layerCol[ily];
    for(int ishx=0; ishx<m_layer.barShowerXCol.size(); ishx++){
      m_barshower0_layer.push_back(m_layer.getDlayer());
      m_barshower0X.push_back( m_layer.barShowerXCol[ishx].getPos().x() );
      m_barshower0Y.push_back( m_layer.barShowerXCol[ishx].getPos().y() );
      m_barshower0Z.push_back( m_layer.barShowerXCol[ishx].getPos().z() );
      m_barshower0E.push_back( m_layer.barShowerXCol[ishx].getE() );
    }
    for(int ishy=0; ishy<m_layer.barShowerYCol.size(); ishy++){
      m_barshower1_layer.push_back(m_layer.getDlayer());
      m_barshower1X.push_back( m_layer.barShowerYCol[ishy].getPos().x() );
      m_barshower1Y.push_back( m_layer.barShowerYCol[ishy].getPos().y() );
      m_barshower1Z.push_back( m_layer.barShowerYCol[ishy].getPos().z() );
      m_barshower1E.push_back( m_layer.barShowerYCol[ishy].getE() );
    }
  }
  t_dataColIter0->Fill();

/*  ClearIter1(); 
  m_layerCol = m_DataCol.LayerCol;
  for(int ily=0; ily<m_layerCol.size(); ily++){
    CRDEcalEDM::CRDCaloLayer m_layer = m_layerCol[ily];
    for(int ishx=0; ishx<m_layer.barShowerXCol.size(); ishx++){
      m_barshower0_layer_iter1.push_back(m_layer.getDlayer());
      m_barshower0X_iter1.push_back( m_layer.barShowerXCol[ishx].getPos().x() );
      m_barshower0Y_iter1.push_back( m_layer.barShowerXCol[ishx].getPos().y() );
      m_barshower0Z_iter1.push_back( m_layer.barShowerXCol[ishx].getPos().z() );
      m_barshower0E_iter1.push_back( m_layer.barShowerXCol[ishx].getE() );
    }
    for(int ishy=0; ishy<m_layer.barShowerYCol.size(); ishy++){
      m_barshower1_layer_iter1.push_back(m_layer.getDlayer());
      m_barshower1X_iter1.push_back( m_layer.barShowerYCol[ishy].getPos().x() );
      m_barshower1Y_iter1.push_back( m_layer.barShowerYCol[ishy].getPos().y() );
      m_barshower1Z_iter1.push_back( m_layer.barShowerYCol[ishy].getPos().z() );
      m_barshower1E_iter1.push_back( m_layer.barShowerYCol[ishy].getE() );
    }
  }
  t_dataColIter1->Fill();
*/

  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_3dclusCol_iter0 = m_DataCol.Clus3DCol;
  for(int icl=0; icl<m_3dclusCol_iter0.size(); icl++){
  ClearIter1();
    CRDEcalEDM::CRDCaloHit3DCluster m_clus = m_3dclusCol_iter0[icl];

    for(int ig=0; ig<m_clus.get2DShowers().size(); ig++){
      m_Iter1_gclus_2dshx.push_back(m_clus.get2DShowers()[ig].getPos().x() );
      m_Iter1_gclus_2dshy.push_back(m_clus.get2DShowers()[ig].getPos().y() );
      m_Iter1_gclus_2dshz.push_back(m_clus.get2DShowers()[ig].getPos().z() );
      m_Iter1_gclus_2dshE.push_back(m_clus.get2DShowers()[ig].getShowerE() );
    }
  t_dataColIter1->Fill();
  }

/*
  //Save ClusTrk info
  std::vector< CRDEcalEDM::CRDCaloHitLongiCluster > m_ClusTrk_X = m_DataCol.LongiClusXCol;
  for(int itrk=0; itrk<m_ClusTrk_X.size(); itrk++){
    //m_trkX_x.push_back(m_ClusTrk_X[itrk].getPos().x());
    //m_trkX_y.push_back(m_ClusTrk_X[itrk].getPos().y());
    //m_trkX_z.push_back(m_ClusTrk_X[itrk].getPos().z());
    //m_trkX_px.push_back(m_ClusTrk_X[itrk].getAxis().x());
    //m_trkX_py.push_back(m_ClusTrk_X[itrk].getAxis().y());
    //m_trkX_pz.push_back(m_ClusTrk_X[itrk].getAxis().z());
    //m_trkX_Nsh.push_back(m_ClusTrk_X[itrk].getBarShowers().size());
    ClearIter();

    for(int is=0; is<m_ClusTrk_X[itrk].getBarShowers().size(); is++){
      m_barshower0_layer.push_back( m_ClusTrk_X[itrk].getBarShowers()[is].getDlayer() );
      m_barshower0X.push_back( m_ClusTrk_X[itrk].getBarShowers()[is].getPos().x() );
      m_barshower0Y.push_back( m_ClusTrk_X[itrk].getBarShowers()[is].getPos().y() );
      m_barshower0Z.push_back( m_ClusTrk_X[itrk].getBarShowers()[is].getPos().z() );
      m_barshower0E.push_back( m_ClusTrk_X[itrk].getBarShowers()[is].getE() );    
    }
    t_dataColIter0->Fill();
  }
  std::vector< CRDEcalEDM::CRDCaloHitLongiCluster > m_ClusTrk_Y = m_DataCol.LongiClusYCol;
  for(int itrk=0; itrk<m_ClusTrk_Y.size(); itrk++){
    m_trkY_x.push_back(m_ClusTrk_Y[itrk].getPos().x());
    m_trkY_y.push_back(m_ClusTrk_Y[itrk].getPos().y());
    m_trkY_z.push_back(m_ClusTrk_Y[itrk].getPos().z());
    m_trkY_px.push_back(m_ClusTrk_Y[itrk].getAxis().x());
    m_trkY_py.push_back(m_ClusTrk_Y[itrk].getAxis().y());
    m_trkY_pz.push_back(m_ClusTrk_Y[itrk].getAxis().z());
    m_trkY_Nsh.push_back(m_ClusTrk_Y[itrk].getBarShowers().size());
  }
  t_dataColIter0->Fill();
*/

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
    int m_type = -1; 
    if(m_clus[i].getStdDevE()<0.05) m_type = 0; 
    else if(m_clus[i].getStdDevE()>0.5) m_type = 1; 
    else m_type = 2; 
    m_Clus_type.push_back(m_type);
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
  t_dataColIter0->Write();
  t_dataColIter1->Write();
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

  delete m_wfile, t_SimBar, t_recoPFO, t_dataColIter0, t_dataColIter1;

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


void PandoraPlusPFAlg::ClearIter0(){
  m_barshower0X.clear(); 
  m_barshower0Y.clear(); 
  m_barshower0Z.clear(); 
  m_barshower0E.clear(); 
  m_barshower0_layer.clear();
  m_barshower1X.clear();
  m_barshower1Y.clear();
  m_barshower1Z.clear();
  m_barshower1E.clear();
  m_barshower1_layer.clear(); 
  m_trkX_x.clear(); 
  m_trkX_y.clear();
  m_trkX_z.clear();
  m_trkX_px.clear();
  m_trkX_py.clear();
  m_trkX_pz.clear();
  m_trkX_Nsh.clear();
  m_trkY_x.clear();
  m_trkY_y.clear();
  m_trkY_z.clear();
  m_trkY_px.clear();
  m_trkY_py.clear();
  m_trkY_pz.clear();
  m_trkY_Nsh.clear();
}


void PandoraPlusPFAlg::ClearIter1(){
  m_barshower0X_iter1.clear();
  m_barshower0Y_iter1.clear();
  m_barshower0Z_iter1.clear();
  m_barshower0E_iter1.clear();
  m_barshower1X_iter1.clear();
  m_barshower1Y_iter1.clear();
  m_barshower1Z_iter1.clear();
  m_barshower1E_iter1.clear();
  m_barshower0_layer_iter1.clear();
  m_barshower1_layer_iter1.clear();

  m_Iter1_gclus_2dshx.clear();
  m_Iter1_gclus_2dshy.clear();
  m_Iter1_gclus_2dshz.clear();
  m_Iter1_gclus_2dshE.clear();
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
