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
   declareProperty("MCParticle", r_MCParticleCol, "Handle of the Input MCParticle collection");
 
   
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
  t_PreRec = new TTree("RecBlock", "RecBlock");
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
  t_SimBar->Branch("simBar_block", &m_simBar_block);
  t_SimBar->Branch("simBar_slayer", &m_simBar_slayer);

  t_PreRec->Branch("PreRec_Bar0x",&m_PreRec_Bar0x);
  t_PreRec->Branch("PreRec_Bar0y",&m_PreRec_Bar0y);
  t_PreRec->Branch("PreRec_Bar0z",&m_PreRec_Bar0z);
  t_PreRec->Branch("PreRec_Bar0E",&m_PreRec_Bar0E);
  t_PreRec->Branch("PreRec_Bar1x",&m_PreRec_Bar1x);
  t_PreRec->Branch("PreRec_Bar1y",&m_PreRec_Bar1y);
  t_PreRec->Branch("PreRec_Bar1z",&m_PreRec_Bar1z);
  t_PreRec->Branch("PreRec_Bar1E",&m_PreRec_Bar1E);
  t_PreRec->Branch("PreRec_shower0X", &m_PreRec_shower0X);
  t_PreRec->Branch("PreRec_shower0Y", &m_PreRec_shower0Y);
  t_PreRec->Branch("PreRec_shower0Z", &m_PreRec_shower0Z);
  t_PreRec->Branch("PreRec_shower0E", &m_PreRec_shower0E);
  t_PreRec->Branch("PreRec_shower0T1", &m_PreRec_shower0T1);
  t_PreRec->Branch("PreRec_shower0T2", &m_PreRec_shower0T2);
  t_PreRec->Branch("PreRec_shower1X", &m_PreRec_shower1X);
  t_PreRec->Branch("PreRec_shower1Y", &m_PreRec_shower1Y);
  t_PreRec->Branch("PreRec_shower1Z", &m_PreRec_shower1Z);
  t_PreRec->Branch("PreRec_shower1E", &m_PreRec_shower1E);
  t_PreRec->Branch("PreRec_shower1T1", &m_PreRec_shower1T1);
  t_PreRec->Branch("PreRec_shower1T2", &m_PreRec_shower1T2);  
  t_PreRec->Branch("PreRec_NshowerX", &m_PreRec_NshowerX);
  t_PreRec->Branch("PreRec_NshowerY", &m_PreRec_NshowerY);
  t_PreRec->Branch("PreRec_NclusterX", &m_PreRec_NclusterX);
  t_PreRec->Branch("PreRec_NclusterY", &m_PreRec_NclusterY);

  t_recoPFO->Branch("Npfo", &m_Npfo);
  t_recoPFO->Branch("recPFO_px", &m_recPFO_px);
  t_recoPFO->Branch("recPFO_py", &m_recPFO_py);
  t_recoPFO->Branch("recPFO_pz", &m_recPFO_pz);
  t_recoPFO->Branch("recPFO_En", &m_recPFO_En);
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
  t_recoPFO->Branch("maPy", &m_mcPy);
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
  const edm4hep::MCParticleCollection* const_MCPCol =  r_MCParticleCol.get();
  if( const_MCPCol!=NULL ) m_pMCParticleCreator->GetMCParticle( m_DataCol, *const_MCPCol);
  m_pTrackCreator->     GetTracks( m_DataCol );
  m_pVertexCreator->    GetVertex( m_DataCol );
  m_pEcalHitsCreator->  GetEcalBars( m_DataCol, *m_edmsvc); 
  m_pHcalHitsCreator->  GetHcalHits( m_DataCol );


  //Link MCParticle-Track and MCParticle-Hit. 
  m_pMCParticleCreator->CreateTrackMCParticleRelation();
  m_pMCParticleCreator->CreateEcalBarMCParticleRelation();
  m_pMCParticleCreator->CreateHcalHitsMCParticleRelation();


  //Perform PFA algorithm
  //m_pEcalClusterRec->Initialize();
  m_pEcalClusterRec->RunAlgorithm( *m_pEcalClusterRecSettings, m_DataCol );
  //m_pEcalClusterRec->ClearAlgorithm();

  //m_DataCol.Print3DClus();

  std::cout<<"Reconstruction is done. Create PFO"<<std::endl;
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

  //Save output information
  ClearBar();
  std::vector<CRDEcalEDM::CRDCaloBlock>  tmp_blockvec = m_DataCol.blockVec;
  for(int ibl=0;ibl<tmp_blockvec.size();ibl++){
    CRDEcalEDM::CRDCaloBlock tmp_barcol = tmp_blockvec[ibl];
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
      m_simBar_block.push_back(hitbar.getBlock());
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
      m_simBar_block.push_back(hitbar.getBlock());
      m_simBar_slayer.push_back(hitbar.getSlayer());
    }
  }
  t_SimBar->Fill();

  std::vector<CRDEcalEDM::CRDCaloLayer>  m_layercol = m_DataCol.LayerCol;
  for(int il=0;il<m_layercol.size();il++){
  
  ClearPreRec();
  std::vector<CRDEcalEDM::CRDCaloBarShower> barShowerXCol = m_layercol[il].barShowerXCol;
  std::vector<CRDEcalEDM::CRDCaloBarShower> barShowerYCol = m_layercol[il].barShowerYCol;
  std::vector<CRDEcalEDM::CRDCaloBar> barColX = m_layercol[il].barXCol;
  std::vector<CRDEcalEDM::CRDCaloBar> barColY = m_layercol[il].barYCol;

  for(int i=0;i<barColX.size();i++){
    m_PreRec_Bar0x.push_back(barColX[i].getPosition().x());
    m_PreRec_Bar0y.push_back(barColX[i].getPosition().y());
    m_PreRec_Bar0z.push_back(barColX[i].getPosition().z());
    m_PreRec_Bar0E.push_back(barColX[i].getEnergy());
  }
  for(int i=0;i<barColY.size();i++){
    m_PreRec_Bar1x.push_back(barColY[i].getPosition().x());
    m_PreRec_Bar1y.push_back(barColY[i].getPosition().y());
    m_PreRec_Bar1z.push_back(barColY[i].getPosition().z());
    m_PreRec_Bar1E.push_back(barColY[i].getEnergy());
  }
  for(int i=0;i<barShowerXCol.size();i++){
    m_PreRec_shower0X.push_back(barShowerXCol[i].getPos().x());
    m_PreRec_shower0Y.push_back(barShowerXCol[i].getPos().y());
    m_PreRec_shower0Z.push_back(barShowerXCol[i].getPos().z());
    m_PreRec_shower0E.push_back(barShowerXCol[i].getE());
    m_PreRec_shower0T2.push_back(barShowerXCol[i].getT2());
    m_PreRec_shower0T1.push_back(barShowerXCol[i].getT1());
  }
  for(int i=0;i<barShowerYCol.size();i++){
    m_PreRec_shower1X.push_back(barShowerYCol[i].getPos().x());
    m_PreRec_shower1Y.push_back(barShowerYCol[i].getPos().y());
    m_PreRec_shower1Z.push_back(barShowerYCol[i].getPos().z());
    m_PreRec_shower1E.push_back(barShowerYCol[i].getE());
    m_PreRec_shower1T1.push_back(barShowerYCol[i].getT1());
    m_PreRec_shower1T2.push_back(barShowerYCol[i].getT2());
  }
  m_PreRec_NshowerX = barShowerXCol.size();
  m_PreRec_NshowerY = barShowerYCol.size();

  t_PreRec->Fill();
  }
  

  std::cout<<"Save PFO into ntuple"<<std::endl;
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
  m_N2dshInClus=0;
  for(int i=0;i<m_clus.size();i++){
    m_N2dshInClus += m_clus[i].get2DShowers().size();
    if(m_clus.size()==1){
	 std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_shower2dCol = m_clus[i].get2DShowers();
    for(int is=0; is<m_N2dshInClus; is++){
      m_2DShower_x.push_back( m_shower2dCol[is].getPos().x() );
      m_2DShower_y.push_back( m_shower2dCol[is].getPos().y() );
      m_2DShower_z.push_back( m_shower2dCol[is].getPos().z() );
      m_2DShower_E.push_back( m_shower2dCol[is].getShowerE() );
    }}
  }
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
  //t_PreRec->Write();
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
  m_simBar_block.clear();
  m_simBar_slayer.clear();
}

void PandoraPlusPFAlg::ClearPreRec(){
  m_PreRec_Bar0x.clear();
  m_PreRec_Bar0y.clear();
  m_PreRec_Bar0z.clear();
  m_PreRec_Bar0E.clear();
  m_PreRec_Bar1x.clear();
  m_PreRec_Bar1y.clear();
  m_PreRec_Bar1z.clear();
  m_PreRec_Bar1E.clear();
  m_PreRec_shower0X.clear();
  m_PreRec_shower0Y.clear();
  m_PreRec_shower0Z.clear();
  m_PreRec_shower0E.clear();
  m_PreRec_shower0T1.clear();
  m_PreRec_shower0T2.clear();
  m_PreRec_shower1X.clear();
  m_PreRec_shower1Y.clear();
  m_PreRec_shower1Z.clear();
  m_PreRec_shower1E.clear();
  m_PreRec_shower1T1.clear();
  m_PreRec_shower1T2.clear();
  m_PreRec_NshowerX=-999;
  m_PreRec_NshowerY=-999;
  m_PreRec_NclusterX=-999;
  m_PreRec_NclusterY=-999;
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
  m_mcPx.clear();
  m_mcPy.clear();
  m_mcPz.clear();
  m_mcEn.clear();
  m_Npfo=-99;
  m_Nmc=-99;
  m_N2dshInClus=-99;
}
#endif
