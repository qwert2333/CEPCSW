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
  t_HoughClusters = new TTree("RecClusters", "RecClusters");
  t_ShowersX = new TTree("RecShowersX", "RecShowersX");
  t_ShowersY = new TTree("RecShowersY", "RecShowersY");
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
  t_SimBar->Branch("simBar_id", &m_simBar_id);

  t_Layers->Branch("NshowerX", &m_NshowerX);
  t_Layers->Branch("NshowerY", &m_NshowerY);
  t_Layers->Branch("barShowerX_x", &m_barShowerX_x);
  t_Layers->Branch("barShowerX_y", &m_barShowerX_y);
  t_Layers->Branch("barShowerX_z", &m_barShowerX_z);
  t_Layers->Branch("barShowerX_E", &m_barShowerX_E);
  t_Layers->Branch("barShowerY_x", &m_barShowerY_x);
  t_Layers->Branch("barShowerY_y", &m_barShowerY_y);
  t_Layers->Branch("barShowerY_z", &m_barShowerY_z);
  t_Layers->Branch("barShowerY_E", &m_barShowerY_E);

  t_HoughClusters->Branch("NclusX", &m_NclusX);
  t_HoughClusters->Branch("NclusY", &m_NclusY);
  t_HoughClusters->Branch("clusX_x", &m_clusX_x);
  t_HoughClusters->Branch("clusX_y", &m_clusX_y);
  t_HoughClusters->Branch("clusX_z", &m_clusX_z);
  t_HoughClusters->Branch("clusX_E", &m_clusX_E);
  t_HoughClusters->Branch("clusX_alpha", &m_clusX_alpha);
  t_HoughClusters->Branch("clusX_rho", &m_clusX_rho);
  t_HoughClusters->Branch("clusX_inter", &m_clusX_inter);
  t_HoughClusters->Branch("clusX_px", &m_clusX_px);
  t_HoughClusters->Branch("clusX_py", &m_clusX_py);
  t_HoughClusters->Branch("clusX_pz", &m_clusX_pz);
  t_HoughClusters->Branch("clusX_Nhit", &m_clusX_Nhit);
  t_HoughClusters->Branch("clusX_Nlayer", &m_clusX_Nlayer);
  t_HoughClusters->Branch("clusY_x", &m_clusY_x);
  t_HoughClusters->Branch("clusY_y", &m_clusY_y);
  t_HoughClusters->Branch("clusY_z", &m_clusY_z);
  t_HoughClusters->Branch("clusY_E", &m_clusY_E);
  t_HoughClusters->Branch("clusY_alpha", &m_clusY_alpha);
  t_HoughClusters->Branch("clusY_rho", &m_clusY_rho);
  t_HoughClusters->Branch("clusY_inter", &m_clusY_inter);
  t_HoughClusters->Branch("clusY_px", &m_clusY_px);
  t_HoughClusters->Branch("clusY_py", &m_clusY_py);
  t_HoughClusters->Branch("clusY_pz", &m_clusY_pz);
  t_HoughClusters->Branch("clusY_Nhit", &m_clusY_Nhit);
  t_HoughClusters->Branch("clusY_Nlayer", &m_clusY_Nlayer);
  //t_HoughClusters->Branch("NshowerX", &m_NshowerX);
  //t_HoughClusters->Branch("NshowerY", &m_NshowerY);

/*  t_HoughClusters->Branch("barShowerX_x", &m_barShowerX_x);
  t_HoughClusters->Branch("barShowerX_y", &m_barShowerX_y);
  t_HoughClusters->Branch("barShowerX_z", &m_barShowerX_z);
  t_HoughClusters->Branch("barShowerX_E", &m_barShowerX_E);
  t_HoughClusters->Branch("barShowerY_x", &m_barShowerY_x);
  t_HoughClusters->Branch("barShowerY_y", &m_barShowerY_y);
  t_HoughClusters->Branch("barShowerY_z", &m_barShowerY_z);
  t_HoughClusters->Branch("barShowerY_E", &m_barShowerY_E);
*/

/*  t_Showers->Branch("Nshower2D", &m_Nshower2D);
  t_Showers->Branch("shower2D_dlayer", &m_shower2D_dlayer);
  t_Showers->Branch("shower2D_part", &m_shower2D_part);
  t_Showers->Branch("shower2D_stave", &m_shower2D_stave);
  t_Showers->Branch("shower2D_module", &m_shower2D_module);
  t_Showers->Branch("shower2D_type", &m_shower2D_type);
  t_Showers->Branch("shower2D_x", &m_shower2D_x);
  t_Showers->Branch("shower2D_y", &m_shower2D_y);
  t_Showers->Branch("shower2D_z", &m_shower2D_z);
  t_Showers->Branch("shower2D_E", &m_shower2D_E);
*/
  t_ShowersX->Branch("barShowerX_x", &m_barShowerX_x);
  t_ShowersX->Branch("barShowerX_y", &m_barShowerX_y);
  t_ShowersX->Branch("barShowerX_z", &m_barShowerX_z);
  t_ShowersX->Branch("barShowerX_E", &m_barShowerX_E);
  t_ShowersY->Branch("barShowerY_x", &m_barShowerY_x);
  t_ShowersY->Branch("barShowerY_y", &m_barShowerY_y);
  t_ShowersY->Branch("barShowerY_z", &m_barShowerY_z);
  t_ShowersY->Branch("barShowerY_E", &m_barShowerY_E);

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
      m_simBar_id.push_back(hitbar.getcellID());
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
  ClearLayer();
  tmp_stavevec = m_DataCol.BlockVec;
  for(int ibl=0;ibl<tmp_stavevec.size();ibl++){
    CRDEcalEDM::CRDCaloBlock tmp_barcol = tmp_stavevec[ibl];
    m_NshowerX = tmp_barcol.getShowerXCol().size();
    m_NshowerY = tmp_barcol.getShowerYCol().size();
    for(int is=0; is<m_NshowerX; is++){
      m_barShowerX_x.push_back( tmp_barcol.getShowerXCol()[is].getPos().x() );
      m_barShowerX_y.push_back( tmp_barcol.getShowerXCol()[is].getPos().y() );
      m_barShowerX_z.push_back( tmp_barcol.getShowerXCol()[is].getPos().z() );
      m_barShowerX_E.push_back( tmp_barcol.getShowerXCol()[is].getE() );
    }
    for(int is=0; is<m_NshowerY; is++){
      m_barShowerY_x.push_back( tmp_barcol.getShowerYCol()[is].getPos().x() );
      m_barShowerY_y.push_back( tmp_barcol.getShowerYCol()[is].getPos().y() );
      m_barShowerY_z.push_back( tmp_barcol.getShowerYCol()[is].getPos().z() );
      m_barShowerY_E.push_back( tmp_barcol.getShowerYCol()[is].getE() );
    }
  }
  t_Layers->Fill();


  //Save tower info
  std::vector<CRDEcalEDM::CRDCaloTower> m_towers = m_DataCol.TowerCol;
  ClearCluster();
  for(int it=0; it<m_towers.size(); it++){

    //Save bar showers
/*    ClearLayer();
    std::vector<CRDEcalEDM::CRDCaloBlock> m_blockcol = m_towers[it].getBlocks();

    for(int ib=0; ib<m_blockcol.size(); ib++){
    m_NshowerX = m_blockcol[ib].getShowerXCol().size();
    m_NshowerY = m_blockcol[ib].getShowerYCol().size();
    for(int is=0; is<m_NshowerX; is++){
      m_barShowerX_x.push_back( m_blockcol[ib].getShowerXCol()[is].getPos().x() );
      m_barShowerX_y.push_back( m_blockcol[ib].getShowerXCol()[is].getPos().y() );
      m_barShowerX_z.push_back( m_blockcol[ib].getShowerXCol()[is].getPos().z() );
      m_barShowerX_E.push_back( m_blockcol[ib].getShowerXCol()[is].getE() );
    }
    for(int is=0; is<m_NshowerY; is++){
      m_barShowerY_x.push_back( m_blockcol[ib].getShowerYCol()[is].getPos().x() );
      m_barShowerY_y.push_back( m_blockcol[ib].getShowerYCol()[is].getPos().y() );
      m_barShowerY_z.push_back( m_blockcol[ib].getShowerYCol()[is].getPos().z() );
      m_barShowerY_E.push_back( m_blockcol[ib].getShowerYCol()[is].getE() );
    }
    }
*/  

    //ClearCluster();
    //Save Hough clusters
    m_NclusX = m_towers[it].getLongiClusterXCol().size(); 
    m_NclusY = m_towers[it].getLongiClusterYCol().size();
    if(m_NclusX==0 && m_NclusY==0) continue;
    for(int ic=0; ic<m_NclusX; ic++){
      CRDEcalEDM::CRDCaloHitLongiCluster m_longicl = m_towers[it].getLongiClusterXCol()[ic];
      m_clusX_x.push_back(m_longicl.getPos().X()); 
      m_clusX_y.push_back(m_longicl.getPos().Y()); 
      m_clusX_z.push_back(m_longicl.getPos().Z()); 
      m_clusX_E.push_back(m_longicl.getE()); 
      m_clusX_alpha.push_back(m_longicl.getHoughAlpha()); 
      m_clusX_rho.push_back(m_longicl.getHoughRho()); 
      m_clusX_inter.push_back(m_longicl.getHoughIntercept());
      m_clusX_px.push_back(m_longicl.getAxis().X());
      m_clusX_py.push_back(m_longicl.getAxis().Y());
      m_clusX_pz.push_back(m_longicl.getAxis().Z());
      m_clusX_Nhit.push_back(m_longicl.getBarShowers().size());
      m_clusX_Nlayer.push_back(m_longicl.getEndDlayer()-m_longicl.getBeginningDlayer());
    }
    for(int ic=0; ic<m_NclusY; ic++){
      CRDEcalEDM::CRDCaloHitLongiCluster m_longicl = m_towers[it].getLongiClusterYCol()[ic];
      m_clusY_x.push_back(m_longicl.getPos().X());
      m_clusY_y.push_back(m_longicl.getPos().Y());
      m_clusY_z.push_back(m_longicl.getPos().Z());
      m_clusY_E.push_back(m_longicl.getE());
      m_clusY_alpha.push_back(m_longicl.getHoughAlpha());
      m_clusY_rho.push_back(m_longicl.getHoughRho());
      m_clusY_inter.push_back(m_longicl.getHoughIntercept());
      m_clusY_px.push_back(m_longicl.getAxis().X());
      m_clusY_py.push_back(m_longicl.getAxis().Y());
      m_clusY_pz.push_back(m_longicl.getAxis().Z());
      m_clusY_Nhit.push_back(m_longicl.getBarShowers().size());
      m_clusY_Nlayer.push_back(m_longicl.getEndDlayer()-m_longicl.getBeginningDlayer());
    }
    //t_HoughClusters->Fill();
  } 
  t_HoughClusters->Fill();

  //Save Hough Clusters
  m_towers = m_DataCol.TowerCol;
  for(int it=0; it<m_towers.size(); it++){

    ClearLayer();

    m_NclusX = m_towers[it].getLongiClusterXCol().size();
    m_NclusY = m_towers[it].getLongiClusterYCol().size();

    for(int is=0; is<m_NclusX; is++){
      m_barShowerX_x.clear(); m_barShowerX_y.clear(); m_barShowerX_z.clear(); m_barShowerX_E.clear();

      CRDEcalEDM::CRDCaloHitLongiCluster m_longicl = m_towers[it].getLongiClusterXCol()[is];
      for(int js=0; js<m_longicl.getBarShowers().size(); js++){
        m_barShowerX_x.push_back( m_longicl.getBarShowers()[js].getPos().x() );
        m_barShowerX_y.push_back( m_longicl.getBarShowers()[js].getPos().y() );
        m_barShowerX_z.push_back( m_longicl.getBarShowers()[js].getPos().z() );
        m_barShowerX_E.push_back( m_longicl.getBarShowers()[js].getE() );
      }
      t_ShowersX->Fill();
    }
    for(int is=0; is<m_NclusY; is++){
      m_barShowerY_x.clear(); m_barShowerY_y.clear(); m_barShowerY_z.clear(); m_barShowerY_E.clear();

      CRDEcalEDM::CRDCaloHitLongiCluster m_longicl = m_towers[it].getLongiClusterYCol()[is];
      for(int js=0; js<m_longicl.getBarShowers().size(); js++){
        m_barShowerY_x.push_back( m_longicl.getBarShowers()[js].getPos().x() );
        m_barShowerY_y.push_back( m_longicl.getBarShowers()[js].getPos().y() );
        m_barShowerY_z.push_back( m_longicl.getBarShowers()[js].getPos().z() );
        m_barShowerY_E.push_back( m_longicl.getBarShowers()[js].getE() );
      }
      t_ShowersY->Fill();
    }

   
  }

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
  t_Layers->Write();
  t_HoughClusters->Write();
  t_ShowersX->Write();
  t_ShowersY->Write();
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

  delete m_wfile, t_SimBar, t_recoPFO, t_Layers, t_HoughClusters, t_ShowersX, t_ShowersY;

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
  m_simBar_id.clear();
}


void PandoraPlusPFAlg::ClearLayer(){
  m_NshowerX=-99;
  m_NshowerY=-99;
  m_barShowerX_x.clear();
  m_barShowerX_y.clear();
  m_barShowerX_z.clear();
  m_barShowerX_E.clear();
  m_barShowerY_x.clear();
  m_barShowerY_y.clear();
  m_barShowerY_z.clear();
  m_barShowerY_E.clear();
}


void PandoraPlusPFAlg::ClearCluster(){
  m_NclusX=-99;
  m_NclusY=-99;
  m_clusX_x.clear();
  m_clusX_y.clear();
  m_clusX_z.clear();
  m_clusX_E.clear();
  m_clusX_px.clear();
  m_clusX_py.clear();
  m_clusX_pz.clear();
  m_clusX_alpha.clear();
  m_clusX_rho.clear();
  m_clusX_inter.clear();
  m_clusX_Nhit.clear();
  m_clusX_Nlayer.clear();
  m_clusY_x.clear();
  m_clusY_y.clear();
  m_clusY_z.clear();
  m_clusY_E.clear();
  m_clusY_px.clear();
  m_clusY_py.clear();
  m_clusY_pz.clear();
  m_clusY_alpha.clear();
  m_clusY_rho.clear();
  m_clusY_inter.clear();
  m_clusY_Nhit.clear();
  m_clusY_Nlayer.clear();

}

void PandoraPlusPFAlg::ClearShower(){
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
