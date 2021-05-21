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
  t_dataColIter0 = new TTree("dataColIter0", "dataColIter0");
  t_dataColIter1 = new TTree("dataColIter1", "dataColIter1");
  t_dataColIter2 = new TTree("dataColIter2", "dataColIter2");
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

  t_dataColIter0->Branch("barx", &m_Iter0_barx);
  t_dataColIter0->Branch("bary", &m_Iter0_bary);
  t_dataColIter0->Branch("barz", &m_Iter0_barz);
  t_dataColIter0->Branch("barE", &m_Iter0_barE);
  t_dataColIter0->Branch("barslayer", &m_Iter0_barslayer);
  t_dataColIter0->Branch("shower0E", &m_Iter0_shower0E);
  t_dataColIter0->Branch("shower0X", &m_Iter0_shower0X);
  t_dataColIter0->Branch("shower0Y", &m_Iter0_shower0Y);
  t_dataColIter0->Branch("shower0Z", &m_Iter0_shower0Z);
  t_dataColIter0->Branch("shower0T1", &m_Iter0_shower0T1);
  t_dataColIter0->Branch("shower0T2", &m_Iter0_shower0T2);
  t_dataColIter0->Branch("shower1E", &m_Iter0_shower1E);
  t_dataColIter0->Branch("shower1X", &m_Iter0_shower1X);
  t_dataColIter0->Branch("shower1Y", &m_Iter0_shower1Y);
  t_dataColIter0->Branch("shower1Z", &m_Iter0_shower1Z);
  t_dataColIter0->Branch("shower1T1", &m_Iter0_shower1T1);
  t_dataColIter0->Branch("shower1T2", &m_Iter0_shower1T2);
  t_dataColIter0->Branch("ly_dlayer", &m_Iter0_ly_dlayer );
  t_dataColIter0->Branch("ly_part", &m_Iter0_ly_part );
  t_dataColIter0->Branch("ly_stave", &m_Iter0_ly_stave );
  t_dataColIter0->Branch("ly_module", &m_Iter0_ly_module );
  t_dataColIter0->Branch("ly_x", &m_Iter0_ly_x );
  t_dataColIter0->Branch("ly_y", &m_Iter0_ly_y );
  t_dataColIter0->Branch("ly_z", &m_Iter0_ly_z );
  t_dataColIter0->Branch("ly_NshowerX", &m_Iter0_ly_NshowerX );
  t_dataColIter0->Branch("ly_NshowerY", &m_Iter0_ly_NshowerY );
  t_dataColIter0->Branch("ly_NclusterX", &m_Iter0_ly_NclusterX );
  t_dataColIter0->Branch("ly_NclusterY", &m_Iter0_ly_NclusterY );
  t_dataColIter0->Branch("Nexpsh", &m_Iter0_Nexpsh );
  t_dataColIter0->Branch("Expsh_x", &m_Iter0_Expsh_x );
  t_dataColIter0->Branch("Expsh_y", &m_Iter0_Expsh_y );
  t_dataColIter0->Branch("Expsh_z", &m_Iter0_Expsh_z );
  t_dataColIter0->Branch("Expsh_E", &m_Iter0_Expsh_E );
  t_dataColIter0->Branch("N2dshower", &m_Iter0_N2dshower);
  t_dataColIter0->Branch("2ds_dlayer", &m_Iter0_2ds_dlayer );
  t_dataColIter0->Branch("2ds_part", &m_Iter0_2ds_part );
  t_dataColIter0->Branch("2ds_stave", &m_Iter0_2ds_stave );
  t_dataColIter0->Branch("2ds_x", &m_Iter0_2ds_x );
  t_dataColIter0->Branch("2ds_y", &m_Iter0_2ds_y );
  t_dataColIter0->Branch("2ds_z", &m_Iter0_2ds_z );
  t_dataColIter0->Branch("2ds_E", &m_Iter0_2ds_E );
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

  t_dataColIter1->Branch("barx", &m_Iter1_barx);
  t_dataColIter1->Branch("bary", &m_Iter1_bary);
  t_dataColIter1->Branch("barz", &m_Iter1_barz);
  t_dataColIter1->Branch("barE", &m_Iter1_barE);
  t_dataColIter1->Branch("barslayer", &m_Iter1_barslayer);
  t_dataColIter1->Branch("shower0E", &m_Iter1_shower0E);
  t_dataColIter1->Branch("shower0X", &m_Iter1_shower0X);
  t_dataColIter1->Branch("shower0Y", &m_Iter1_shower0Y);
  t_dataColIter1->Branch("shower0Z", &m_Iter1_shower0Z);
  t_dataColIter1->Branch("shower0T1", &m_Iter1_shower0T1);
  t_dataColIter1->Branch("shower0T2", &m_Iter1_shower0T2);
  t_dataColIter1->Branch("shower1E", &m_Iter1_shower1E);
  t_dataColIter1->Branch("shower1X", &m_Iter1_shower1X);
  t_dataColIter1->Branch("shower1Y", &m_Iter1_shower1Y);
  t_dataColIter1->Branch("shower1Z", &m_Iter1_shower1Z);
  t_dataColIter1->Branch("shower1T1", &m_Iter1_shower1T1);
  t_dataColIter1->Branch("shower1T2", &m_Iter1_shower1T2);
  t_dataColIter1->Branch("ly_dlayer", &m_Iter1_ly_dlayer );
  t_dataColIter1->Branch("ly_part", &m_Iter1_ly_part );
  t_dataColIter1->Branch("ly_stave", &m_Iter1_ly_stave );
  t_dataColIter1->Branch("ly_module", &m_Iter1_ly_module );
  t_dataColIter1->Branch("ly_x", &m_Iter1_ly_x );
  t_dataColIter1->Branch("ly_y", &m_Iter1_ly_y );
  t_dataColIter1->Branch("ly_z", &m_Iter1_ly_z );
  t_dataColIter1->Branch("ly_NshowerX", &m_Iter1_ly_NshowerX );
  t_dataColIter1->Branch("ly_NshowerY", &m_Iter1_ly_NshowerY );
  t_dataColIter1->Branch("ly_NclusterX", &m_Iter1_ly_NclusterX );
  t_dataColIter1->Branch("ly_NclusterY", &m_Iter1_ly_NclusterY );
  t_dataColIter1->Branch("Nexpsh", &m_Iter1_Nexpsh );
  t_dataColIter1->Branch("Expsh_x", &m_Iter1_Expsh_x );
  t_dataColIter1->Branch("Expsh_y", &m_Iter1_Expsh_y );
  t_dataColIter1->Branch("Expsh_z", &m_Iter1_Expsh_z );
  t_dataColIter1->Branch("Expsh_E", &m_Iter1_Expsh_E );
  t_dataColIter1->Branch("N2dshower", &m_Iter1_N2dshower);
  t_dataColIter1->Branch("2ds_dlayer", &m_Iter1_2ds_dlayer );
  t_dataColIter1->Branch("2ds_part", &m_Iter1_2ds_part );
  t_dataColIter1->Branch("2ds_stave", &m_Iter1_2ds_stave );
  t_dataColIter1->Branch("2ds_x", &m_Iter1_2ds_x );
  t_dataColIter1->Branch("2ds_y", &m_Iter1_2ds_y );
  t_dataColIter1->Branch("2ds_z", &m_Iter1_2ds_z );
  t_dataColIter1->Branch("2ds_E", &m_Iter1_2ds_E );
  t_dataColIter1->Branch("Ngoodclus", &m_Iter1_Ngoodclus );
  t_dataColIter1->Branch("Nbadclus", &m_Iter1_Nbadclus );
  t_dataColIter1->Branch("clus_Nly", &m_Iter1_clus_Nly );
  t_dataColIter1->Branch("clus_x", &m_Iter1_clus_x );
  t_dataColIter1->Branch("clus_y", &m_Iter1_clus_y );
  t_dataColIter1->Branch("clus_z", &m_Iter1_clus_z );
  t_dataColIter1->Branch("clus_E", &m_Iter1_clus_E );
  t_dataColIter1->Branch("clus_px", &m_Iter1_clus_px );
  t_dataColIter1->Branch("clus_py", &m_Iter1_clus_py );
  t_dataColIter1->Branch("clus_pz", &m_Iter1_clus_pz );
  t_dataColIter1->Branch("gclus_2dshx", &m_Iter1_gclus_2dshx );
  t_dataColIter1->Branch("gclus_2dshy", &m_Iter1_gclus_2dshy );
  t_dataColIter1->Branch("gclus_2dshz", &m_Iter1_gclus_2dshz );
  t_dataColIter1->Branch("gclus_2dshE", &m_Iter1_gclus_2dshE );
  t_dataColIter1->Branch("bclus_2dshx", &m_Iter1_bclus_2dshx );
  t_dataColIter1->Branch("bclus_2dshy", &m_Iter1_bclus_2dshy );
  t_dataColIter1->Branch("bclus_2dshz", &m_Iter1_bclus_2dshz );
  t_dataColIter1->Branch("bclus_2dshE", &m_Iter1_bclus_2dshE );

  t_dataColIter2->Branch("barx", &m_Iter2_barx);
  t_dataColIter2->Branch("bary", &m_Iter2_bary);
  t_dataColIter2->Branch("barz", &m_Iter2_barz);
  t_dataColIter2->Branch("barE", &m_Iter2_barE);
  t_dataColIter2->Branch("barslayer", &m_Iter2_barslayer);
  t_dataColIter2->Branch("shower0E", &m_Iter2_shower0E);
  t_dataColIter2->Branch("shower0X", &m_Iter2_shower0X);
  t_dataColIter2->Branch("shower0Y", &m_Iter2_shower0Y);
  t_dataColIter2->Branch("shower0Z", &m_Iter2_shower0Z);
  t_dataColIter2->Branch("shower0T1", &m_Iter2_shower0T1);
  t_dataColIter2->Branch("shower0T2", &m_Iter2_shower0T2);
  t_dataColIter2->Branch("shower1E", &m_Iter2_shower1E);
  t_dataColIter2->Branch("shower1X", &m_Iter2_shower1X);
  t_dataColIter2->Branch("shower1Y", &m_Iter2_shower1Y);
  t_dataColIter2->Branch("shower1Z", &m_Iter2_shower1Z);
  t_dataColIter2->Branch("shower1T1", &m_Iter2_shower1T1);
  t_dataColIter2->Branch("shower1T2", &m_Iter2_shower1T2);
  t_dataColIter2->Branch("ly_dlayer", &m_Iter2_ly_dlayer );
  t_dataColIter2->Branch("ly_part", &m_Iter2_ly_part );
  t_dataColIter2->Branch("ly_stave", &m_Iter2_ly_stave );
  t_dataColIter2->Branch("ly_module", &m_Iter2_ly_module );
  t_dataColIter2->Branch("ly_x", &m_Iter2_ly_x );
  t_dataColIter2->Branch("ly_y", &m_Iter2_ly_y );
  t_dataColIter2->Branch("ly_z", &m_Iter2_ly_z );
  t_dataColIter2->Branch("ly_NshowerX", &m_Iter2_ly_NshowerX );
  t_dataColIter2->Branch("ly_NshowerY", &m_Iter2_ly_NshowerY );
  t_dataColIter2->Branch("ly_NclusterX", &m_Iter2_ly_NclusterX );
  t_dataColIter2->Branch("ly_NclusterY", &m_Iter2_ly_NclusterY );
  t_dataColIter2->Branch("Nexpsh", &m_Iter2_Nexpsh );
  t_dataColIter2->Branch("Expsh_x", &m_Iter2_Expsh_x );
  t_dataColIter2->Branch("Expsh_y", &m_Iter2_Expsh_y );
  t_dataColIter2->Branch("Expsh_z", &m_Iter2_Expsh_z );
  t_dataColIter2->Branch("Expsh_E", &m_Iter2_Expsh_E );
  t_dataColIter2->Branch("N2dshower", &m_Iter2_N2dshower);
  t_dataColIter2->Branch("2ds_dlayer", &m_Iter2_2ds_dlayer );
  t_dataColIter2->Branch("2ds_part", &m_Iter2_2ds_part );
  t_dataColIter2->Branch("2ds_stave", &m_Iter2_2ds_stave );
  t_dataColIter2->Branch("2ds_x", &m_Iter2_2ds_x );
  t_dataColIter2->Branch("2ds_y", &m_Iter2_2ds_y );
  t_dataColIter2->Branch("2ds_z", &m_Iter2_2ds_z );
  t_dataColIter2->Branch("2ds_E", &m_Iter2_2ds_E );
  t_dataColIter2->Branch("Ngoodclus", &m_Iter2_Ngoodclus );
  t_dataColIter2->Branch("Nbadclus", &m_Iter2_Nbadclus );
  t_dataColIter2->Branch("clus_Nly", &m_Iter2_clus_Nly );
  t_dataColIter2->Branch("clus_x", &m_Iter2_clus_x );
  t_dataColIter2->Branch("clus_y", &m_Iter2_clus_y );
  t_dataColIter2->Branch("clus_z", &m_Iter2_clus_z );
  t_dataColIter2->Branch("clus_E", &m_Iter2_clus_E );
  t_dataColIter2->Branch("clus_px", &m_Iter2_clus_px );
  t_dataColIter2->Branch("clus_py", &m_Iter2_clus_py );
  t_dataColIter2->Branch("clus_pz", &m_Iter2_clus_pz );
  t_dataColIter2->Branch("gclus_2dshx", &m_Iter2_gclus_2dshx );
  t_dataColIter2->Branch("gclus_2dshy", &m_Iter2_gclus_2dshy );
  t_dataColIter2->Branch("gclus_2dshz", &m_Iter2_gclus_2dshz );
  t_dataColIter2->Branch("gclus_2dshE", &m_Iter2_gclus_2dshE );
  t_dataColIter2->Branch("bclus_2dshx", &m_Iter2_bclus_2dshx );
  t_dataColIter2->Branch("bclus_2dshy", &m_Iter2_bclus_2dshy );
  t_dataColIter2->Branch("bclus_2dshz", &m_Iter2_bclus_2dshz );
  t_dataColIter2->Branch("bclus_2dshE", &m_Iter2_bclus_2dshE );

  t_recoPFO->Branch("Npfo", &m_Npfo);
  t_recoPFO->Branch("recPFO_px", &m_recPFO_px);
  t_recoPFO->Branch("recPFO_py", &m_recPFO_py);
  t_recoPFO->Branch("recPFO_pz", &m_recPFO_pz);
  t_recoPFO->Branch("recPFO_En", &m_recPFO_En);
  t_recoPFO->Branch("N3dclus", &m_N3dclus);
  t_recoPFO->Branch("N2dshInClus", &m_N2dshInClus);
  t_recoPFO->Branch("Shower2D_x", &m_2DShower_x);
  t_recoPFO->Branch("Shower2D_y", &m_2DShower_y);
  t_recoPFO->Branch("Shower2D_z", &m_2DShower_z);
  t_recoPFO->Branch("Shower2D_E", &m_2DShower_E);
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

/*  std::cout<<"Good cluster fit chi2: "<<std::endl;
  int NgoodCl = m_DataCol.GoodClus3DCol.size();
  for(int igc=0;igc<NgoodCl; igc++){
    m_DataCol.GoodClus3DCol[igc].FitProfile();
    std::cout<<"Good Cluster #"<<igc<<'\t'<<m_DataCol.GoodClus3DCol[igc].getChi2()<<std::endl;
  }
*/

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

  //Save intermediate info
  //Iteration0 info
/*  std::vector<CRDEcalEDM::CRDCaloBlock> m_staveCol_iter0 = m_DataCol.BlockVec_iter0;
  std::vector<CRDEcalEDM::CRDCaloLayer> m_layerCol_iter0 = m_DataCol.LayerCol_iter0;
  for(int il=0;il<m_layerCol_iter0.size();il++){
    ClearIter0();
    CRDEcalEDM::CRDCaloLayer m_layer = m_layerCol_iter0[il];

    CRDEcalEDM::CRDCaloBlock tmp_barcol; tmp_barcol.Clear();
    for(int ib=0;ib<m_staveCol_iter0.size(); ib++){
      if(m_layer.getDlayer()==m_staveCol_iter0[ib].getDlayer()){ tmp_barcol = m_staveCol_iter0[ib]; break; }
    }
    if(tmp_barcol.getBarXCol().size()==0 && tmp_barcol.getBarYCol().size()==0) continue;

    m_Iter0_ly_NshowerX.push_back( m_layer.barShowerXCol.size());
    m_Iter0_ly_NshowerY.push_back( m_layer.barShowerYCol.size());
    m_Iter0_ly_NclusterX.push_back( m_layer.barClusXCol.size());
    m_Iter0_ly_NclusterY.push_back( m_layer.barClusYCol.size());
    m_Iter0_Nexpsh.push_back( tmp_barcol.getEMCandidateCol().size());

    for(int ib=0; ib<m_layer.barXCol.size(); ib++){
      m_Iter0_barx.push_back(m_layer.barXCol[ib].getPosition().x());
      m_Iter0_bary.push_back(m_layer.barXCol[ib].getPosition().y());
      m_Iter0_barz.push_back(m_layer.barXCol[ib].getPosition().z());
      m_Iter0_barE.push_back( m_layer.barXCol[ib].getEnergy() );
      m_Iter0_barslayer.push_back(m_layer.barXCol[ib].getSlayer());
    }
    for(int ib=0; ib<m_layer.barYCol.size(); ib++){
      m_Iter0_barx.push_back(m_layer.barYCol[ib].getPosition().x());
      m_Iter0_bary.push_back(m_layer.barYCol[ib].getPosition().y());
      m_Iter0_barz.push_back(m_layer.barYCol[ib].getPosition().z());
      m_Iter0_barE.push_back(m_layer.barYCol[ib].getEnergy());
      m_Iter0_barslayer.push_back(m_layer.barYCol[ib].getSlayer());
    }

    for(int ie=0; ie<tmp_barcol.getEMCandidateCol().size(); ie++){
      m_Iter0_Expsh_x.push_back(tmp_barcol.getEMCandidateCol()[ie].ExpPos.x());
      m_Iter0_Expsh_y.push_back(tmp_barcol.getEMCandidateCol()[ie].ExpPos.y());
      m_Iter0_Expsh_z.push_back(tmp_barcol.getEMCandidateCol()[ie].ExpPos.z());
      m_Iter0_Expsh_E.push_back(tmp_barcol.getEMCandidateCol()[ie].ExpEshower);
    }

    for(int is=0; is<m_layer.barShowerXCol.size(); is++){
      m_Iter0_shower0E.push_back(m_layer.barShowerXCol[is].getE());
      m_Iter0_shower0X.push_back(m_layer.barShowerXCol[is].getPos().x());
      m_Iter0_shower0Y.push_back(m_layer.barShowerXCol[is].getPos().y());
      m_Iter0_shower0Z.push_back(m_layer.barShowerXCol[is].getPos().z());
      m_Iter0_shower0T1.push_back(m_layer.barShowerXCol[is].getT1());
      m_Iter0_shower0T2.push_back(m_layer.barShowerXCol[is].getT2());
    }
    for(int is=0; is<m_layer.barShowerYCol.size(); is++){
      m_Iter0_shower1E.push_back(m_layer.barShowerYCol[is].getE());
      m_Iter0_shower1X.push_back(m_layer.barShowerYCol[is].getPos().x());
      m_Iter0_shower1Y.push_back(m_layer.barShowerYCol[is].getPos().y());
      m_Iter0_shower1Z.push_back(m_layer.barShowerYCol[is].getPos().z());
      m_Iter0_shower1T1.push_back(m_layer.barShowerYCol[is].getT1());
      m_Iter0_shower1T2.push_back(m_layer.barShowerYCol[is].getT2());
    }   
  t_dataColIter0->Fill();
  }
*/

  ClearIter0();
  std::vector<CRDEcalEDM::CRDCaloLayer> m_layerCol_iter0 = m_DataCol.LayerCol_iter0;
  for(int il=0; il<m_layerCol_iter0.size(); il++){
    CRDEcalEDM::CRDCaloLayer m_layer = m_layerCol_iter0[il];
    m_Iter0_ly_NshowerX.push_back( m_layer.barShowerXCol.size() );
    m_Iter0_ly_NshowerY.push_back( m_layer.barShowerYCol.size() );
    m_Iter0_ly_NclusterX.push_back( m_layer.barClusXCol.size() );
    m_Iter0_ly_NclusterY.push_back( m_layer.barClusYCol.size() );

    CRDEcalEDM::CRDCaloBar tmp_bar;
    if(m_layer.barXCol.size()!=0) tmp_bar = m_layer.barXCol[0];
    else if(m_layer.barYCol.size()!=0 ) tmp_bar = m_layer.barYCol[0];
    else continue;
    m_Iter0_ly_dlayer.push_back(tmp_bar.getDlayer());
    m_Iter0_ly_part.push_back(tmp_bar.getPart());
    m_Iter0_ly_stave.push_back(tmp_bar.getStave());
    m_Iter0_ly_module.push_back(tmp_bar.getModule());
    m_Iter0_ly_x.push_back(tmp_bar.getPosition().x());
    m_Iter0_ly_y.push_back(tmp_bar.getPosition().y());
    m_Iter0_ly_z.push_back(tmp_bar.getPosition().z());
  }

  std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_2dshCol_iter0 = m_DataCol.Shower2DCol_iter0;
  m_Iter0_N2dshower = m_2dshCol_iter0.size();
  for(int i2s=0; i2s<m_Iter0_N2dshower; i2s++){
    CRDEcalEDM::CRDCaloHit2DShower m_2dshower = m_2dshCol_iter0[i2s];
    m_Iter0_2ds_dlayer.push_back(m_2dshower.getDlayer());
    m_Iter0_2ds_part.push_back(m_2dshower.getPart());
    m_Iter0_2ds_stave.push_back(m_2dshower.getStave());
    m_Iter0_2ds_x.push_back(m_2dshower.getPos().x());
    m_Iter0_2ds_y.push_back(m_2dshower.getPos().y());
    m_Iter0_2ds_z.push_back(m_2dshower.getPos().z());
    m_Iter0_2ds_E.push_back(m_2dshower.getShowerE());
  }

  m_Iter0_Ngoodclus = m_DataCol.GoodClus3DCol_iter0.size();
  m_Iter0_Nbadclus = m_DataCol.BadClus3DCol_iter0.size();
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_3dclusCol_iter0 = m_DataCol.GoodClus3DCol_iter0;
  for(int icl=0; icl<m_3dclusCol_iter0.size(); icl++){
    CRDEcalEDM::CRDCaloHit3DShower m_clus = m_3dclusCol_iter0[icl];
    m_Iter0_clus_Nly.push_back(m_clus.get2DShowers().size());
    m_Iter0_clus_x.push_back(m_clus.getShowerCenter().x());
    m_Iter0_clus_y.push_back(m_clus.getShowerCenter().y());
    m_Iter0_clus_z.push_back(m_clus.getShowerCenter().z());
    m_Iter0_clus_E.push_back(m_clus.getShowerE());
    m_Iter0_clus_px.push_back(m_clus.getAxis().x());
    m_Iter0_clus_py.push_back(m_clus.getAxis().y());
    m_Iter0_clus_pz.push_back(m_clus.getAxis().z());
  }
  m_3dclusCol_iter0.clear(); 
  m_3dclusCol_iter0 = m_DataCol.BadClus3DCol_iter0;
  for(int icl=0; icl<m_3dclusCol_iter0.size(); icl++){
    CRDEcalEDM::CRDCaloHit3DShower m_clus = m_3dclusCol_iter0[icl];
    m_Iter0_clus_Nly.push_back(m_clus.get2DShowers().size());
    m_Iter0_clus_x.push_back(m_clus.getShowerCenter().x());
    m_Iter0_clus_y.push_back(m_clus.getShowerCenter().y());
    m_Iter0_clus_z.push_back(m_clus.getShowerCenter().z());
    m_Iter0_clus_E.push_back(m_clus.getShowerE());
    m_Iter0_clus_px.push_back(m_clus.getAxis().x());
    m_Iter0_clus_py.push_back(m_clus.getAxis().y());
    m_Iter0_clus_pz.push_back(m_clus.getAxis().z());
  }
  t_dataColIter0->Fill();


  //Iteration1 info
/*  std::vector<CRDEcalEDM::CRDCaloBlock> m_staveCol_iter1 = m_DataCol.BlockVec_iter1;
  std::vector<CRDEcalEDM::CRDCaloLayer> m_layerCol_iter1 = m_DataCol.LayerCol_iter1;
  for(int il=0;il<m_layerCol_iter1.size();il++){
    ClearIter1();
    CRDEcalEDM::CRDCaloLayer m_layer = m_layerCol_iter1[il];

    CRDEcalEDM::CRDCaloBlock tmp_barcol; tmp_barcol.Clear();
    for(int ib=0;ib<m_staveCol_iter1.size(); ib++){
      if(m_layer.getDlayer()==m_staveCol_iter1[ib].getDlayer()){ tmp_barcol = m_staveCol_iter1[ib]; break; }
    }
    if(tmp_barcol.getBarXCol().size()==0 && tmp_barcol.getBarYCol().size()==0) continue;

    m_Iter1_ly_NshowerX.push_back( m_layer.barShowerXCol.size() );
    m_Iter1_ly_NshowerY.push_back( m_layer.barShowerYCol.size() );
    m_Iter1_ly_NclusterX.push_back( m_layer.barClusXCol.size() );
    m_Iter1_ly_NclusterY.push_back( m_layer.barClusYCol.size() );
    m_Iter1_Nexpsh.push_back( tmp_barcol.getEMCandidateCol().size() );

    for(int ib=0; ib<m_layer.barXCol.size(); ib++){
      m_Iter1_barx.push_back(m_layer.barXCol[ib].getPosition().x());
      m_Iter1_bary.push_back(m_layer.barXCol[ib].getPosition().y());
      m_Iter1_barz.push_back(m_layer.barXCol[ib].getPosition().z());
      m_Iter1_barE.push_back(m_layer.barXCol[ib].getEnergy());
      m_Iter1_barslayer.push_back(m_layer.barXCol[ib].getSlayer());
    }
    for(int ib=0; ib<m_layer.barYCol.size(); ib++){
      m_Iter1_barx.push_back(m_layer.barYCol[ib].getPosition().x());
      m_Iter1_bary.push_back(m_layer.barYCol[ib].getPosition().y());
      m_Iter1_barz.push_back(m_layer.barYCol[ib].getPosition().z());
      m_Iter1_barE.push_back(m_layer.barYCol[ib].getEnergy());
      m_Iter1_barslayer.push_back(m_layer.barYCol[ib].getSlayer());
    }

    for(int ie=0; ie<tmp_barcol.getEMCandidateCol().size(); ie++){
      m_Iter1_Expsh_x.push_back(tmp_barcol.getEMCandidateCol()[ie].ExpPos.x());
      m_Iter1_Expsh_y.push_back(tmp_barcol.getEMCandidateCol()[ie].ExpPos.y());
      m_Iter1_Expsh_z.push_back(tmp_barcol.getEMCandidateCol()[ie].ExpPos.z());
      m_Iter1_Expsh_E.push_back(tmp_barcol.getEMCandidateCol()[ie].ExpEshower);
    }

    for(int is=0; is<m_layer.barShowerXCol.size(); is++){
      m_Iter1_shower0E.push_back(m_layer.barShowerXCol[is].getE());
      m_Iter1_shower0X.push_back(m_layer.barShowerXCol[is].getPos().x());
      m_Iter1_shower0Y.push_back(m_layer.barShowerXCol[is].getPos().y());
      m_Iter1_shower0Z.push_back(m_layer.barShowerXCol[is].getPos().z());
      m_Iter1_shower0T1.push_back(m_layer.barShowerXCol[is].getT1());
      m_Iter1_shower0T2.push_back(m_layer.barShowerXCol[is].getT2());
    }
    for(int is=0; is<m_layer.barShowerYCol.size(); is++){
      m_Iter1_shower1E.push_back(m_layer.barShowerYCol[is].getE());
      m_Iter1_shower1X.push_back(m_layer.barShowerYCol[is].getPos().x());
      m_Iter1_shower1Y.push_back(m_layer.barShowerYCol[is].getPos().y());
      m_Iter1_shower1Z.push_back(m_layer.barShowerYCol[is].getPos().z());
      m_Iter1_shower1T1.push_back(m_layer.barShowerYCol[is].getT1());
      m_Iter1_shower1T2.push_back(m_layer.barShowerYCol[is].getT2());
    }
  t_dataColIter1->Fill();
  }
*/
  ClearIter1();
  std::vector<CRDEcalEDM::CRDCaloLayer> m_layerCol_iter1 = m_DataCol.LayerCol_iter1;
  for(int il=0; il<m_layerCol_iter1.size(); il++){
    CRDEcalEDM::CRDCaloLayer m_layer = m_layerCol_iter1[il];
    m_Iter1_ly_NshowerX.push_back( m_layer.barShowerXCol.size() );
    m_Iter1_ly_NshowerY.push_back( m_layer.barShowerYCol.size() );
    m_Iter1_ly_NclusterX.push_back( m_layer.barClusXCol.size() );
    m_Iter1_ly_NclusterY.push_back( m_layer.barClusYCol.size() );

    CRDEcalEDM::CRDCaloBar tmp_bar; 
    if(m_layer.barXCol.size()!=0) tmp_bar = m_layer.barXCol[0];
    else if(m_layer.barYCol.size()!=0 ) tmp_bar = m_layer.barYCol[0];
    else continue;
    m_Iter1_ly_dlayer.push_back(tmp_bar.getDlayer());
    m_Iter1_ly_part.push_back(tmp_bar.getPart());
    m_Iter1_ly_stave.push_back(tmp_bar.getStave());
    m_Iter1_ly_module.push_back(tmp_bar.getModule());
    m_Iter1_ly_x.push_back(tmp_bar.getPosition().x());
    m_Iter1_ly_y.push_back(tmp_bar.getPosition().y());
    m_Iter1_ly_z.push_back(tmp_bar.getPosition().z());
  }

  std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_2dshCol_iter1 = m_DataCol.Shower2DCol_iter1;
  m_Iter1_N2dshower = m_2dshCol_iter1.size();
  for(int i2s=0; i2s<m_Iter1_N2dshower; i2s++){
    CRDEcalEDM::CRDCaloHit2DShower m_2dshower = m_2dshCol_iter1[i2s];
    m_Iter1_2ds_dlayer.push_back(m_2dshower.getDlayer());
    m_Iter1_2ds_part.push_back(m_2dshower.getPart());
    m_Iter1_2ds_stave.push_back(m_2dshower.getStave());
    m_Iter1_2ds_x.push_back(m_2dshower.getPos().x());
    m_Iter1_2ds_y.push_back(m_2dshower.getPos().y());
    m_Iter1_2ds_z.push_back(m_2dshower.getPos().z());
    m_Iter1_2ds_E.push_back(m_2dshower.getShowerE());
  }

  m_Iter1_Ngoodclus = m_DataCol.GoodClus3DCol_iter1.size();
  m_Iter1_Nbadclus = m_DataCol.BadClus3DCol_iter1.size();
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_3dclusCol_iter1 = m_DataCol.GoodClus3DCol_iter1;
  for(int icl=0; icl<m_3dclusCol_iter1.size(); icl++){
    CRDEcalEDM::CRDCaloHit3DShower m_clus = m_3dclusCol_iter1[icl];
    m_Iter1_clus_Nly.push_back(m_clus.get2DShowers().size());
    m_Iter1_clus_x.push_back(m_clus.getShowerCenter().x());
    m_Iter1_clus_y.push_back(m_clus.getShowerCenter().y());
    m_Iter1_clus_z.push_back(m_clus.getShowerCenter().z());
    m_Iter1_clus_E.push_back(m_clus.getShowerE());
    m_Iter1_clus_px.push_back(m_clus.getAxis().x());
    m_Iter1_clus_py.push_back(m_clus.getAxis().y());
    m_Iter1_clus_pz.push_back(m_clus.getAxis().z());

    for(int ig=0; ig<m_clus.get2DShowers().size(); ig++){
      m_Iter1_gclus_2dshx.push_back(m_clus.get2DShowers()[ig].getPos().x() );
      m_Iter1_gclus_2dshy.push_back(m_clus.get2DShowers()[ig].getPos().y() );
      m_Iter1_gclus_2dshz.push_back(m_clus.get2DShowers()[ig].getPos().z() );
      m_Iter1_gclus_2dshE.push_back(m_clus.get2DShowers()[ig].getShowerE() );
    }
  }
  m_3dclusCol_iter1.clear();
  m_3dclusCol_iter1 = m_DataCol.BadClus3DCol_iter1;
  for(int icl=0; icl<m_3dclusCol_iter1.size(); icl++){
    CRDEcalEDM::CRDCaloHit3DShower m_clus = m_3dclusCol_iter1[icl];
    m_Iter1_clus_Nly.push_back(m_clus.get2DShowers().size());
    m_Iter1_clus_x.push_back(m_clus.getShowerCenter().x());
    m_Iter1_clus_y.push_back(m_clus.getShowerCenter().y());
    m_Iter1_clus_z.push_back(m_clus.getShowerCenter().z());
    m_Iter1_clus_E.push_back(m_clus.getShowerE());
    m_Iter1_clus_px.push_back(m_clus.getAxis().x());
    m_Iter1_clus_py.push_back(m_clus.getAxis().y());
    m_Iter1_clus_pz.push_back(m_clus.getAxis().z());

    for(int ig=0; ig<m_clus.get2DShowers().size(); ig++){
      m_Iter1_bclus_2dshx.push_back(m_clus.get2DShowers()[ig].getPos().x() );
      m_Iter1_bclus_2dshy.push_back(m_clus.get2DShowers()[ig].getPos().y() );
      m_Iter1_bclus_2dshz.push_back(m_clus.get2DShowers()[ig].getPos().z() );
      m_Iter1_bclus_2dshE.push_back(m_clus.get2DShowers()[ig].getShowerE() );
    }
  }
  t_dataColIter1->Fill();


  ClearIter2();
  std::vector<CRDEcalEDM::CRDCaloLayer> m_layerCol_iter2 = m_DataCol.LayerCol;
  for(int il=0; il<m_layerCol_iter2.size(); il++){
    CRDEcalEDM::CRDCaloLayer m_layer = m_layerCol_iter2[il];
    m_Iter2_ly_NshowerX.push_back( m_layer.barShowerXCol.size() );
    m_Iter2_ly_NshowerY.push_back( m_layer.barShowerYCol.size() );
    m_Iter2_ly_NclusterX.push_back( m_layer.barClusXCol.size() );
    m_Iter2_ly_NclusterY.push_back( m_layer.barClusYCol.size() );

    CRDEcalEDM::CRDCaloBar tmp_bar; 
    if(m_layer.barXCol.size()!=0) tmp_bar = m_layer.barXCol[0];
    else if(m_layer.barYCol.size()!=0 ) tmp_bar = m_layer.barYCol[0];
    else continue;
    m_Iter2_ly_dlayer.push_back(tmp_bar.getDlayer());
    m_Iter2_ly_part.push_back(tmp_bar.getPart());
    m_Iter2_ly_stave.push_back(tmp_bar.getStave());
    m_Iter2_ly_module.push_back(tmp_bar.getModule());
    m_Iter2_ly_x.push_back(tmp_bar.getPosition().x());
    m_Iter2_ly_y.push_back(tmp_bar.getPosition().y());
    m_Iter2_ly_z.push_back(tmp_bar.getPosition().z());
  }

  std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_2dshCol_iter2 = m_DataCol.Shower2DCol;
  m_Iter2_N2dshower = m_2dshCol_iter2.size();
  for(int i2s=0; i2s<m_Iter2_N2dshower; i2s++){
    CRDEcalEDM::CRDCaloHit2DShower m_2dshower = m_2dshCol_iter2[i2s];
    m_Iter2_2ds_dlayer.push_back(m_2dshower.getDlayer());
    m_Iter2_2ds_part.push_back(m_2dshower.getPart());
    m_Iter2_2ds_stave.push_back(m_2dshower.getStave());
    m_Iter2_2ds_x.push_back(m_2dshower.getPos().x());
    m_Iter2_2ds_y.push_back(m_2dshower.getPos().y());
    m_Iter2_2ds_z.push_back(m_2dshower.getPos().z());
    m_Iter2_2ds_E.push_back(m_2dshower.getShowerE());
  }
  m_Iter2_Ngoodclus = m_DataCol.GoodClus3DCol.size();
  m_Iter2_Nbadclus = m_DataCol.BadClus3DCol.size();
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_3dclusCol_iter2 = m_DataCol.GoodClus3DCol;
  for(int icl=0; icl<m_3dclusCol_iter2.size(); icl++){
    CRDEcalEDM::CRDCaloHit3DShower m_clus = m_3dclusCol_iter2[icl];
    m_Iter2_clus_Nly.push_back(m_clus.get2DShowers().size());
    m_Iter2_clus_x.push_back(m_clus.getShowerCenter().x());
    m_Iter2_clus_y.push_back(m_clus.getShowerCenter().y());
    m_Iter2_clus_z.push_back(m_clus.getShowerCenter().z());
    m_Iter2_clus_E.push_back(m_clus.getShowerE());
    m_Iter2_clus_px.push_back(m_clus.getAxis().x());
    m_Iter2_clus_py.push_back(m_clus.getAxis().y());
    m_Iter2_clus_pz.push_back(m_clus.getAxis().z());

    for(int ig=0; ig<m_clus.get2DShowers().size(); ig++){
      m_Iter2_gclus_2dshx.push_back(m_clus.get2DShowers()[ig].getPos().x() );
      m_Iter2_gclus_2dshy.push_back(m_clus.get2DShowers()[ig].getPos().y() );
      m_Iter2_gclus_2dshz.push_back(m_clus.get2DShowers()[ig].getPos().z() );
      m_Iter2_gclus_2dshE.push_back(m_clus.get2DShowers()[ig].getShowerE() );
    }
  }
  m_3dclusCol_iter2.clear();
  m_3dclusCol_iter2 = m_DataCol.BadClus3DCol;
  for(int icl=0; icl<m_3dclusCol_iter2.size(); icl++){
    CRDEcalEDM::CRDCaloHit3DShower m_clus = m_3dclusCol_iter2[icl];
    m_Iter2_clus_Nly.push_back(m_clus.get2DShowers().size());
    m_Iter2_clus_x.push_back(m_clus.getShowerCenter().x());
    m_Iter2_clus_y.push_back(m_clus.getShowerCenter().y());
    m_Iter2_clus_z.push_back(m_clus.getShowerCenter().z());
    m_Iter2_clus_E.push_back(m_clus.getShowerE());
    m_Iter2_clus_px.push_back(m_clus.getAxis().x());
    m_Iter2_clus_py.push_back(m_clus.getAxis().y());
    m_Iter2_clus_pz.push_back(m_clus.getAxis().z());

    for(int ig=0; ig<m_clus.get2DShowers().size(); ig++){
      m_Iter2_bclus_2dshx.push_back(m_clus.get2DShowers()[ig].getPos().x() );
      m_Iter2_bclus_2dshy.push_back(m_clus.get2DShowers()[ig].getPos().y() );
      m_Iter2_bclus_2dshz.push_back(m_clus.get2DShowers()[ig].getPos().z() );
      m_Iter2_bclus_2dshE.push_back(m_clus.get2DShowers()[ig].getShowerE() );
    }
  }
  t_dataColIter2->Fill();
  
  //Save PFO information
  ClearRecPFO();
  std::vector< CRDEcalEDM::PFObject > m_pfoCol = m_DataCol.PFOCol;
  m_Npfo = m_pfoCol.size(); 
  for(int ipfo=0;ipfo<m_Npfo;ipfo++){
    m_recPFO_px.push_back(m_pfoCol[ipfo].getP4().Px());
    m_recPFO_py.push_back(m_pfoCol[ipfo].getP4().Py());
    m_recPFO_pz.push_back(m_pfoCol[ipfo].getP4().Pz());
    m_recPFO_En.push_back(m_pfoCol[ipfo].getEnergy()); 
  }

  std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_clus = m_DataCol.Clus3DCol;
  m_N3dclus = m_clus.size();
  for(int i=0;i<m_N3dclus;i++){ 
    m_N2dshInClus.push_back(m_clus[i].get2DShowers().size());
    m_2DShower_x.push_back(m_clus[i].getShowerCenter().x());
    m_2DShower_y.push_back(m_clus[i].getShowerCenter().y());
    m_2DShower_z.push_back(m_clus[i].getShowerCenter().z());
    m_2DShower_E.push_back(m_clus[i].getShowerE());
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
  t_dataColIter2->Write();
  t_recoPFO->Write();
  m_wfile->Close();

/*  delete m_pMCParticleCreator;
  delete m_pTrackCreator;
  delete m_pVertexCreator;
  delete m_pEcalHitsCreator;
  delete m_pHcalHitsCreator;
  delete m_pPfoCreator;

  delete m_pMCParticleCreatorSettings;
  delete m_pTrackCreatorSettings;
  delete m_pVertexCreatorSettings;
  delete m_EcalHitsCreatorSettings;
  delete m_pHcalHitsCreatorSettings;
  delete m_pPfoCreatorSettings;

  delete m_pEcalClusterRec;
  delete m_pEcalClusterRecSettings;

  delete m_wfile, t_SimBar, t_PreRec, t_recoPFO;
*/
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

void PandoraPlusPFAlg::ClearRecPFO(){
  m_recPFO_px.clear(); 
  m_recPFO_py.clear(); 
  m_recPFO_pz.clear(); 
  m_recPFO_En.clear(); 
  m_2DShower_x.clear();
  m_2DShower_y.clear();
  m_2DShower_z.clear();
  m_2DShower_E.clear();
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

void PandoraPlusPFAlg::ClearIter0(){
  m_Iter0_barx.clear();
  m_Iter0_bary.clear();
  m_Iter0_barz.clear();
  m_Iter0_barE.clear();
  m_Iter0_barslayer.clear();
  m_Iter0_shower0E.clear();
  m_Iter0_shower0X.clear();
  m_Iter0_shower0Y.clear();
  m_Iter0_shower0Z.clear();
  m_Iter0_shower0T1.clear();
  m_Iter0_shower0T2.clear();
  m_Iter0_shower1E.clear();
  m_Iter0_shower1X.clear();
  m_Iter0_shower1Y.clear();
  m_Iter0_shower1Z.clear();
  m_Iter0_shower1T1.clear();
  m_Iter0_shower1T2.clear();
  m_Iter0_ly_dlayer.clear();
  m_Iter0_ly_part.clear();
  m_Iter0_ly_stave.clear();
  m_Iter0_ly_module.clear();
  m_Iter0_ly_x.clear();
  m_Iter0_ly_y.clear();
  m_Iter0_ly_z.clear();
  m_Iter0_ly_NshowerX.clear();
  m_Iter0_ly_NshowerY.clear();
  m_Iter0_ly_NclusterX.clear();
  m_Iter0_ly_NclusterY.clear();
  m_Iter0_Nexpsh.clear();
  m_Iter0_N2dshower=-99;
  m_Iter0_Expsh_x.clear();
  m_Iter0_Expsh_y.clear();
  m_Iter0_Expsh_z.clear();
  m_Iter0_Expsh_E.clear();
  m_Iter0_2ds_dlayer.clear();
  m_Iter0_2ds_part.clear();
  m_Iter0_2ds_stave.clear();
  m_Iter0_2ds_x.clear();
  m_Iter0_2ds_y.clear();
  m_Iter0_2ds_z.clear();
  m_Iter0_2ds_E.clear();
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
}

void PandoraPlusPFAlg::ClearIter1(){
  m_Iter1_barx.clear();
  m_Iter1_bary.clear();
  m_Iter1_barz.clear();
  m_Iter1_barE.clear();
  m_Iter1_barslayer.clear();
  m_Iter1_shower0E.clear();
  m_Iter1_shower0X.clear();
  m_Iter1_shower0Y.clear();
  m_Iter1_shower0Z.clear();
  m_Iter1_shower0T1.clear();
  m_Iter1_shower0T2.clear();
  m_Iter1_shower1E.clear();
  m_Iter1_shower1X.clear();
  m_Iter1_shower1Y.clear();
  m_Iter1_shower1Z.clear();
  m_Iter1_shower1T1.clear();
  m_Iter1_shower1T2.clear();
  m_Iter1_ly_dlayer.clear();
  m_Iter1_ly_part.clear();
  m_Iter1_ly_stave.clear();
  m_Iter1_ly_module.clear();
  m_Iter1_ly_x.clear();
  m_Iter1_ly_y.clear();
  m_Iter1_ly_z.clear();
  m_Iter1_ly_NshowerX.clear();
  m_Iter1_ly_NshowerY.clear();
  m_Iter1_ly_NclusterX.clear();
  m_Iter1_ly_NclusterY.clear();
  m_Iter1_Nexpsh.clear();
  m_Iter1_N2dshower=-99;
  m_Iter1_Expsh_x.clear();
  m_Iter1_Expsh_y.clear();
  m_Iter1_Expsh_z.clear();
  m_Iter1_Expsh_E.clear();
  m_Iter1_2ds_dlayer.clear();
  m_Iter1_2ds_part.clear();
  m_Iter1_2ds_stave.clear();
  m_Iter1_2ds_x.clear();
  m_Iter1_2ds_y.clear();
  m_Iter1_2ds_z.clear();
  m_Iter1_2ds_E.clear();
  m_Iter1_Ngoodclus=-99;
  m_Iter1_Nbadclus=-99;
  m_Iter1_clus_Nly.clear();
  m_Iter1_clus_x.clear();
  m_Iter1_clus_y.clear();
  m_Iter1_clus_z.clear();
  m_Iter1_clus_E.clear();
  m_Iter1_clus_px.clear();
  m_Iter1_clus_py.clear();
  m_Iter1_clus_pz.clear();
  m_Iter1_gclus_2dshx.clear();
  m_Iter1_gclus_2dshy.clear();
  m_Iter1_gclus_2dshz.clear();
  m_Iter1_gclus_2dshE.clear();
  m_Iter1_bclus_2dshx.clear();
  m_Iter1_bclus_2dshy.clear();
  m_Iter1_bclus_2dshz.clear();
  m_Iter1_bclus_2dshE.clear();

}

void PandoraPlusPFAlg::ClearIter2(){
  m_Iter2_barx.clear();
  m_Iter2_bary.clear();
  m_Iter2_barz.clear();
  m_Iter2_barE.clear();
  m_Iter2_barslayer.clear();
  m_Iter2_shower0E.clear();
  m_Iter2_shower0X.clear();
  m_Iter2_shower0Y.clear();
  m_Iter2_shower0Z.clear();
  m_Iter2_shower0T1.clear();
  m_Iter2_shower0T2.clear();
  m_Iter2_shower1E.clear();
  m_Iter2_shower1X.clear();
  m_Iter2_shower1Y.clear();
  m_Iter2_shower1Z.clear();
  m_Iter2_shower1T1.clear();
  m_Iter2_shower1T2.clear();
  m_Iter2_ly_dlayer.clear();
  m_Iter2_ly_part.clear();
  m_Iter2_ly_stave.clear();
  m_Iter2_ly_module.clear();
  m_Iter2_ly_x.clear();
  m_Iter2_ly_y.clear();
  m_Iter2_ly_z.clear();
  m_Iter2_ly_NshowerX.clear();
  m_Iter2_ly_NshowerY.clear();
  m_Iter2_ly_NclusterX.clear();
  m_Iter2_ly_NclusterY.clear();
  m_Iter2_Nexpsh.clear();
  m_Iter2_N2dshower=-99;
  m_Iter2_Expsh_x.clear();
  m_Iter2_Expsh_y.clear();
  m_Iter2_Expsh_z.clear();
  m_Iter2_Expsh_E.clear();
  m_Iter2_2ds_dlayer.clear();
  m_Iter2_2ds_part.clear();
  m_Iter2_2ds_stave.clear();
  m_Iter2_2ds_x.clear();
  m_Iter2_2ds_y.clear();
  m_Iter2_2ds_z.clear();
  m_Iter2_2ds_E.clear();
  m_Iter2_Ngoodclus=-99;
  m_Iter2_Nbadclus=-99;
  m_Iter2_clus_Nly.clear();
  m_Iter2_clus_x.clear();
  m_Iter2_clus_y.clear();
  m_Iter2_clus_z.clear();
  m_Iter2_clus_E.clear();
  m_Iter2_clus_px.clear();
  m_Iter2_clus_py.clear();
  m_Iter2_clus_pz.clear();
  m_Iter2_gclus_2dshx.clear();
  m_Iter2_gclus_2dshy.clear();
  m_Iter2_gclus_2dshz.clear();
  m_Iter2_gclus_2dshE.clear();
  m_Iter2_bclus_2dshx.clear();
  m_Iter2_bclus_2dshy.clear();
  m_Iter2_bclus_2dshz.clear();
  m_Iter2_bclus_2dshE.clear();
}
#endif
