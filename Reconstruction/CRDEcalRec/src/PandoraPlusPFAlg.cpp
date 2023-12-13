#ifndef PANDORAPLUS_ALG_C
#define PANDORAPLUS_ALG_C

#include "PandoraPlusPFAlg.h"

// #include <fstream>
// #include <ctime>

using namespace std;
using namespace dd4hep;

int PandoraPlus::CaloUnit::Nmodule = 8;
int PandoraPlus::CaloUnit::Npart = 4;
int PandoraPlus::CaloUnit::Nstave = 11;
int PandoraPlus::CaloUnit::Nlayer = 14;
int PandoraPlus::CaloUnit::NbarPhi = 47;
int PandoraPlus::CaloUnit::NbarZ = 60;
int PandoraPlus::CaloUnit::over_module[28] = {13,15,16,18,19,21,22,24,25,26,28,29,30,32,33,35,36,38,39,41,42,43,45,46};
int PandoraPlus::CaloUnit::over_module_set = 2;
float PandoraPlus::CaloUnit::barsize = 10.; //mm
float PandoraPlus::CaloUnit::ecal_innerR = 1860;  //mm

DECLARE_COMPONENT( PandoraPlusPFAlg )

PandoraPlusPFAlg::PandoraPlusPFAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{
 
  // Output collections
  declareProperty("RecCaloHitCollection", w_RecCaloCol, "Handle of Reconstructed CaloHit collection");
  declareProperty("ClusterCollection", w_ClusterCollection, "Handle of Cluster collection");   
  declareProperty("RecoPFOCollection", w_ReconstructedParticleCollection, "Handle of Reconstructed PFO collection");   
  declareProperty("RecoVtxCollection", w_VertexCollection, "Handle of Reconstructed vertex collection");   
  declareProperty("MCRecoPFOAssociationCollection", w_MCRecoParticleAssociationCollection, "Handle of MC-RecoPFO association collection");   

}

StatusCode PandoraPlusPFAlg::initialize()
{
  //Initialize global settings
  m_GlobalSettings.map_floatPars["BField"] = m_BField;
  m_GlobalSettings.map_floatPars["Seed"] = m_seed;
  m_GlobalSettings.map_intPars["Debug"] = m_Debug;
  m_GlobalSettings.map_intPars["SkipEvt"] = m_Nskip;


  //Initialize Creator settings
  m_pMCParticleCreatorSettings.map_stringPars["MCParticleCollections"] = name_MCParticleCol.value();

  m_pTrackCreatorSettings.map_stringVecPars["trackCollections"] = name_TrackCol.value();
  m_pTrackCreatorSettings.map_floatPars["BField"] = m_BField; 

  std::vector<std::string> name_CaloHits = name_EcalHits; 
  std::vector<std::string> name_CaloReadout = name_EcalReadout;
  name_CaloHits.insert( name_CaloHits.end(), name_HcalHits.begin(), name_HcalHits.end() );
  name_CaloReadout.insert(name_CaloReadout.end(), name_HcalReadout.begin(), name_HcalReadout.end());

  m_CaloHitsCreatorSettings.map_stringVecPars["CaloHitCollections"] = name_CaloHits;
  m_CaloHitsCreatorSettings.map_stringPars["EcalType"] = m_EcalType.value();

  //Initialize Creators
  m_pMCParticleCreator = new MCParticleCreator( m_pMCParticleCreatorSettings );
  m_pTrackCreator      = new TrackCreator( m_pTrackCreatorSettings );
  m_pCaloHitsCreator   = new CaloHitsCreator( m_CaloHitsCreatorSettings );
  m_pOutputCreator     = new OutputCreator( m_OutputCreatorSettings );

  //Readin collections

  //---MC particle---
  if(!name_MCParticleCol.empty()) r_MCParticleCol = new DataHandle<edm4hep::MCParticleCollection> (name_MCParticleCol, Gaudi::DataHandle::Reader, this);

  //---Tracks---
  for(auto& _trk : name_TrackCol) if(!_trk.empty()) r_TrackCols.push_back( new TrackType(_trk, Gaudi::DataHandle::Reader, this) );

  //---Calo Hits---
  for(auto& _ecal : name_EcalHits){
    if(!_ecal.empty()){ 
      //r_ECalHitCols.push_back( new CaloType(_ecal, Gaudi::DataHandle::Reader, this) );
      r_CaloHitCols.push_back( new CaloType(_ecal, Gaudi::DataHandle::Reader, this) );
  }}
  for(auto& _hcal : name_HcalHits){ 
    if(!_hcal.empty()){
      //r_HCalHitCols.push_back( new CaloType(_hcal, Gaudi::DataHandle::Reader, this) );
      r_CaloHitCols.push_back( new CaloType(_hcal, Gaudi::DataHandle::Reader, this) );
  }}

  //---MCParticle CaloHit Association
  if(!name_MCPTrkAssoCol.empty())      r_MCPTrkAssoCol = new DataHandle<edm4hep::MCRecoTrackParticleAssociationCollection> (name_MCPTrkAssoCol, Gaudi::DataHandle::Reader, this);

  std::vector<std::string> name_CaloAssoCol = name_EcalMCPAssociation; 
  name_CaloAssoCol.insert(name_CaloAssoCol.end(), name_HcalMCPAssociation.begin(), name_HcalMCPAssociation.end());
  if(name_CaloAssoCol.size()==name_CaloHits.size()){
    for(int iCol=0; iCol<name_CaloAssoCol.size(); iCol++){
      map_CaloMCPAssoCols[name_CaloHits[iCol]] = new CaloParticleAssoType(name_CaloAssoCol[iCol], Gaudi::DataHandle::Reader, this);
    }
  }


  //Register Algorithms
  //--- Initialize algorithm maps ---
  m_algorithmManager.RegisterAlgorithmFactory("ExampleAlg",             new ExampleAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("GlobalClusteringAlg",    new GlobalClusteringAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("HcalClusteringAlg",      new HcalClusteringAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("LocalMaxFindingAlg",     new LocalMaxFindingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("HoughClusteringAlg",     new HoughClusteringAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("TrackMatchingAlg",       new TrackMatchingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("ConeClustering2DAlg",    new ConeClustering2DAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("AxisMergingAlg",         new AxisMergingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("EnergySplittingAlg",     new EnergySplittingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("EnergyTimeMatchingAlg",  new EnergyTimeMatchingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("PFOCreatingAlg",         new PFOCreatingAlg::Factory);
  //m_algorithmManager.RegisterAlgorithmFactory("ConeClusteringAlg",      new ConeClusteringAlg::Factory);
  //m_algorithmManager.RegisterAlgorithmFactory("ConeClusteringAlgHCAL",  new ConeClusteringAlg::Factory);

  m_algorithmManager.RegisterAlgorithmFactory("TruthTrackMatchingAlg",       new TruthTrackMatchingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("TruthPatternRecAlg",          new TruthPatternRecAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("TruthEnergySplittingAlg",     new TruthEnergySplittingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("TruthMatchingAlg",            new TruthMatchingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("TruthClusteringAlg",          new TruthClusteringAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("TruthClusterMergingAlg",      new TruthClusterMergingAlg::Factory);

  //--- Create algorithm from readin settings ---
  for(int ialg=0; ialg<name_Algs.value().size(); ialg++){
    Settings m_settings; 
    for(int ipar=0; ipar<name_AlgPars.value()[ialg].size(); ipar++){
      if(type_AlgPars.value()[ialg].at(ipar)=="int")    m_settings.map_intPars[name_AlgPars.value()[ialg].at(ipar)] = std::stoi( (string)value_AlgPars.value()[ialg].at(ipar) );
      if(type_AlgPars.value()[ialg].at(ipar)=="double") m_settings.map_floatPars[name_AlgPars.value()[ialg].at(ipar)] = std::stod( (string)value_AlgPars.value()[ialg].at(ipar) );
      if(type_AlgPars.value()[ialg].at(ipar)=="string") m_settings.map_stringPars[name_AlgPars.value()[ialg].at(ipar)] = value_AlgPars.value()[ialg].at(ipar) ;
      if(type_AlgPars.value()[ialg].at(ipar)=="bool")   m_settings.map_boolPars[name_AlgPars.value()[ialg].at(ipar)] = (bool)std::stoi( (string)value_AlgPars.value()[ialg].at(ipar) );
    }

    m_algorithmManager.RegisterAlgorithm( name_Algs.value()[ialg], m_settings );
  }


  //Initialize services
  m_geosvc = service<IGeomSvc>("GeomSvc");
  if ( !m_geosvc )  throw "PandoraPlusPFAlg :Failed to find GeomSvc ...";

  for(unsigned int i=0; i<name_CaloReadout.size(); i++){
    if(name_CaloReadout[i].empty()) continue;
    dd4hep::DDSegmentation::BitFieldCoder* tmp_decoder = m_geosvc->getDecoder(name_CaloReadout[i]);
    if (!tmp_decoder) {
      error() << "Failed to get the decoder for: " << name_CaloReadout[i] << endmsg;
      return StatusCode::FAILURE;
    }
    map_readout_decoder[name_CaloHits[i]] = tmp_decoder;
  }

  rndm.SetSeed(m_seed);
  std::cout<<"PandoraPlusPFAlg::initialize"<<std::endl;


  //Output collections
  
  //Output ntuple for analysis.
  if(m_WriteAna){
    std::string s_outfile = m_filename;
    m_wfile = new TFile(s_outfile.c_str(), "recreate");
    t_SimBar = new TTree("SimBarHit", "SimBarHit");
    t_HalfCluster = new TTree("HalfCluster","HalfCluster");
    t_Cluster = new TTree("RecClusters", "RecClusters");
    t_Track = new TTree("RecTracks", "RecTracks");
    t_PFO = new TTree("PFO", "PFO");


    //Bar
    t_SimBar->Branch("totE_EcalSim", &m_totE_EcalSim);
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
    t_SimBar->Branch("simBar_bar", &m_simBar_bar);
    t_SimBar->Branch("totE_HcalSim", &m_totE_HcalSim);
    t_SimBar->Branch("HcalHit_x", &m_HcalHit_x); 
    t_SimBar->Branch("HcalHit_y", &m_HcalHit_y); 
    t_SimBar->Branch("HcalHit_z", &m_HcalHit_z); 
    t_SimBar->Branch("HcalHit_E", &m_HcalHit_E); 
    t_SimBar->Branch("HcalHit_layer", &m_HcalHit_layer); 

    t_HalfCluster->Branch("totE_HFClusV", &m_totE_HFClusV);
    t_HalfCluster->Branch("HalfClusterV_x", &m_HalfClusterV_x);
    t_HalfCluster->Branch("HalfClusterV_y", &m_HalfClusterV_y);
    t_HalfCluster->Branch("HalfClusterV_z", &m_HalfClusterV_z);
    t_HalfCluster->Branch("HalfClusterV_E", &m_HalfClusterV_E);
    t_HalfCluster->Branch("totE_HFClusU", &m_totE_HFClusU);
    t_HalfCluster->Branch("HalfClusterU_x", &m_HalfClusterU_x);
    t_HalfCluster->Branch("HalfClusterU_y", &m_HalfClusterU_y);
    t_HalfCluster->Branch("HalfClusterU_z", &m_HalfClusterU_z);
    t_HalfCluster->Branch("HalfClusterU_E", &m_HalfClusterU_E);  

    //Clusters and MCParticles
    t_Cluster->Branch("Nclus", &m_Nclus);
    t_Cluster->Branch("totE_EcalRec",  &m_totE_EcalRec);
    t_Cluster->Branch("Clus_x", &m_Clus_x);
    t_Cluster->Branch("Clus_y", &m_Clus_y);
    t_Cluster->Branch("Clus_z", &m_Clus_z);
    t_Cluster->Branch("Clus_Px", &m_Clus_Px);
    t_Cluster->Branch("Clus_Py", &m_Clus_Py);
    t_Cluster->Branch("Clus_Pz", &m_Clus_Pz);
    t_Cluster->Branch("Clus_E", &m_Clus_E);
    t_Cluster->Branch("Clus_Ptrk", &m_Clus_Ptrk);
    t_Cluster->Branch("Clus_Ntrk", &m_Clus_Ntrk);
    t_Cluster->Branch("Clus_Nhit", &m_Clus_Nhit);
    t_Cluster->Branch("Clus_truthPDG", &m_Clus_truthPDG);
    t_Cluster->Branch("Clus_truthFrac", &m_Clus_truthFrac);
    t_Cluster->Branch("Clus_startLayer", &m_Clus_startLayer);
    t_Cluster->Branch("Clus_endLayer", &m_Clus_endLayer);
    t_Cluster->Branch("Clus_maxELayer", &m_Clus_maxELayer);
    t_Cluster->Branch("Clus_maxWidthLayer", &m_Clus_maxWidthLayer);
    t_Cluster->Branch("Clus_typeU", &m_Clus_typeU);
    t_Cluster->Branch("Clus_typeV", &m_Clus_typeV);
    t_Cluster->Branch("Clus_width", &m_Clus_width);
    t_Cluster->Branch("Clus_ScndM", &m_Clus_ScndM);
    t_Cluster->Branch("Clus_E1Etot", &m_Clus_E1Etot);
    t_Cluster->Branch("Clus_E2Etot", &m_Clus_E2Etot);
    t_Cluster->Branch("Clus_E5Etot", &m_Clus_E5Etot);
    t_Cluster->Branch("Clus_EhalfEtot", &m_Clus_EhalfEtot);
    t_Cluster->Branch("Clus_EaxisEtot", &m_Clus_EaxisEtot);
    t_Cluster->Branch("Clus_hitx", &m_Clus_hitx);
    t_Cluster->Branch("Clus_hity", &m_Clus_hity);
    t_Cluster->Branch("Clus_hitz", &m_Clus_hitz);
    t_Cluster->Branch("Clus_hitE", &m_Clus_hitE);
    t_Cluster->Branch("Clus_hittag", &m_Clus_hittag);
    t_Cluster->Branch("Clus_hittag_trk", &m_Clus_hittag_trk);
    t_Cluster->Branch("Nmc", &m_Nmc);
    t_Cluster->Branch("mcPdgid",     &m_mcPdgid);
    t_Cluster->Branch("mcStatus",    &m_mcStatus);
    t_Cluster->Branch("mcPx", &m_mcPx);
    t_Cluster->Branch("mcPy", &m_mcPy);
    t_Cluster->Branch("mcPz", &m_mcPz);
    t_Cluster->Branch("mcEn", &m_mcEn);
    t_Cluster->Branch("mcMass", &m_mcMass);
    t_Cluster->Branch("mcCharge", &m_mcCharge);
    t_Cluster->Branch("mcEPx", &m_mcEPx);
    t_Cluster->Branch("mcEPy", &m_mcEPy);
    t_Cluster->Branch("mcEPz", &m_mcEPz);

    t_Cluster->Branch("totE_HcalRec",  &m_totE_HcalRec);
    t_Cluster->Branch("Hcal_clus_x", &m_Hcal_clus_x);
    t_Cluster->Branch("Hcal_clus_y", &m_Hcal_clus_y);
    t_Cluster->Branch("Hcal_clus_z", &m_Hcal_clus_z);
    t_Cluster->Branch("Hcal_clus_E", &m_Hcal_clus_E);
    t_Cluster->Branch("Hcal_hit_tag", &m_Hcal_hit_tag);
    t_Cluster->Branch("Hcal_hit_x", &m_Hcal_hit_x);
    t_Cluster->Branch("Hcal_hit_y", &m_Hcal_hit_y);
    t_Cluster->Branch("Hcal_hit_z", &m_Hcal_hit_z);
    t_Cluster->Branch("Hcal_hit_E", &m_Hcal_hit_E);   


    // Tracks
    t_Track->Branch("m_Ntrk", &m_Ntrk);
    t_Track->Branch("m_type", &m_type);
    t_Track->Branch("m_trkstate_d0", &m_trkstate_d0);
    t_Track->Branch("m_trkstate_z0", &m_trkstate_z0);
    t_Track->Branch("m_trkstate_phi", &m_trkstate_phi);
    t_Track->Branch("m_trkstate_tanL", &m_trkstate_tanL);
    t_Track->Branch("m_trkstate_kappa", &m_trkstate_kappa);
    t_Track->Branch("m_trkstate_omega", &m_trkstate_omega);
    t_Track->Branch("m_trkstate_refx", &m_trkstate_refx);
    t_Track->Branch("m_trkstate_refy", &m_trkstate_refy);
    t_Track->Branch("m_trkstate_refz", &m_trkstate_refz);
    t_Track->Branch("m_trkstate_location", &m_trkstate_location);
    t_Track->Branch("m_trkstate_tag", &m_trkstate_tag);


    //PFOs
    t_PFO->Branch("pfo_tag",     &pfo_tag);
    t_PFO->Branch("n_track",     &n_track);
    t_PFO->Branch("n_ecal_clus", &n_ecal_clus);
    t_PFO->Branch("n_hcal_clus", &n_hcal_clus);
    t_PFO->Branch("trk_pfo_tag", &m_trk_pfo_tag);
    t_PFO->Branch("trk_pt", &m_trk_pt);
    t_PFO->Branch("trk_pz", &m_trk_pz);
    t_PFO->Branch("trk_mcpid", &m_trk_mcpid);
    t_PFO->Branch("trk_mc_px", &m_trk_mc_px);
    t_PFO->Branch("trk_mc_py", &m_trk_mc_py);
    t_PFO->Branch("trk_mc_pz", &m_trk_mc_pz);
    t_PFO->Branch("trk_mc_E", &m_trk_mc_E);
    t_PFO->Branch("ecal_pfo_tag", &m_ecal_pfo_tag);
    t_PFO->Branch("ecal_clus_x", &m_ecal_clus_x);
    t_PFO->Branch("ecal_clus_y", &m_ecal_clus_y);
    t_PFO->Branch("ecal_clus_z", &m_ecal_clus_z);
    t_PFO->Branch("ecal_clus_E", &m_ecal_clus_E);
    t_PFO->Branch("ecal_clus_mcpid", &m_ecal_clus_mcpid);
    t_PFO->Branch("ecal_clus_mc_px", &m_ecal_clus_mc_px);
    t_PFO->Branch("ecal_clus_mc_py", &m_ecal_clus_mc_py);
    t_PFO->Branch("ecal_clus_mc_pz", &m_ecal_clus_mc_pz);
    t_PFO->Branch("ecal_clus_mc_E",  &m_ecal_clus_mc_E);
    t_PFO->Branch("hcal_pfo_tag", &m_hcal_pfo_tag);
    t_PFO->Branch("hcal_clus_x", &m_hcal_clus_x);
    t_PFO->Branch("hcal_clus_y", &m_hcal_clus_y);
    t_PFO->Branch("hcal_clus_z", &m_hcal_clus_z);
    t_PFO->Branch("hcal_clus_E", &m_hcal_clus_E);
    t_PFO->Branch("hcal_clus_mcpid", &m_hcal_clus_mcpid);
    t_PFO->Branch("hcal_clus_mc_px", &m_hcal_clus_mc_px);
    t_PFO->Branch("hcal_clus_mc_py", &m_hcal_clus_mc_py);
    t_PFO->Branch("hcal_clus_mc_pz", &m_hcal_clus_mc_pz);
    t_PFO->Branch("hcal_clus_mc_E",  &m_hcal_clus_mc_E);

  }

  return GaudiAlgorithm::initialize();
}

StatusCode PandoraPlusPFAlg::execute()
{
// clock_t yyy_start, yyy_endrec, yyy_endfill;
// yyy_start = clock(); // 记录开始时间

  if(_nEvt==0) std::cout<<"PandoraPlusPFAlg::execute Start"<<std::endl;
  std::cout<<"Processing event: "<<_nEvt<<std::endl;

  if(_nEvt<m_Nskip){ _nEvt++;  return GaudiAlgorithm::initialize(); }

  //InitializeForNewEvent(); 
  PandoraPlusDataCol     m_DataCol;
  m_DataCol.Clear();

  //Readin collections 
std::cout<<"Readin MCParticle"<<std::endl;
  m_pMCParticleCreator->CreateMCParticle( m_DataCol, *r_MCParticleCol );
cout<<"Readin Tracks"<<endl;
  m_pTrackCreator->CreateTracks( m_DataCol, r_TrackCols, r_MCPTrkAssoCol );
cout<<"Readin CaloHits"<<endl;
  m_pCaloHitsCreator->CreateCaloHits( m_DataCol, r_CaloHitCols, map_readout_decoder, map_CaloMCPAssoCols );

  //Perform PFA algorithm
cout<<"Run Algorithms"<<endl;
  m_algorithmManager.RunAlgorithm( m_DataCol );

cout<<"Create output"<<endl;
cout<<"  Creating CaloHits"<<endl;
  m_pOutputCreator->CreateRecCaloHits( m_DataCol, w_RecCaloCol );
cout<<"  Creating clusters"<<endl;
  m_pOutputCreator->CreateCluster( m_DataCol, w_ClusterCollection );
cout<<"  Creating PFOs"<<endl;
  m_pOutputCreator->CreatePFO( m_DataCol, w_ReconstructedParticleCollection );


// yyy_endrec = clock();  // 重建结束的时间

cout<<"Write tuples"<<endl;
  //---------------------Write Ana tuples-------------------------
  //Save Raw bars information
  ClearBar();
  m_totE_EcalSim = 0.;
  for(int ibar=0;ibar<m_DataCol.map_BarCol["BarCol"].size();ibar++){
    auto p_hitbar = m_DataCol.map_BarCol["BarCol"][ibar].get();
    m_simBar_x.push_back(p_hitbar->getPosition().x());
    m_simBar_y.push_back(p_hitbar->getPosition().y());
    m_simBar_z.push_back(p_hitbar->getPosition().z());
    m_simBar_Q1.push_back(p_hitbar->getQ1());
    m_simBar_Q2.push_back(p_hitbar->getQ2());
    m_simBar_T1.push_back(p_hitbar->getT1());
    m_simBar_T2.push_back(p_hitbar->getT2());
    m_simBar_module.push_back(p_hitbar->getModule());
    m_simBar_dlayer.push_back(p_hitbar->getDlayer());
    m_simBar_part.push_back(p_hitbar->getPart());
    m_simBar_stave.push_back(p_hitbar->getStave());
    m_simBar_slayer.push_back(p_hitbar->getSlayer());
    m_simBar_bar.push_back(p_hitbar->getBar());
    m_totE_EcalSim += (p_hitbar->getQ1()+p_hitbar->getQ2())/2.; 
  }

  std::vector<PandoraPlus::CaloHit*> m_hcalHitsCol; m_hcalHitsCol.clear();
  for(int ih=0; ih<m_DataCol.map_CaloHit["HCALBarrel"].size(); ih++)
    m_hcalHitsCol.push_back( m_DataCol.map_CaloHit["HCALBarrel"][ih].get() );

  m_totE_HcalSim = 0.;
  for(int ihit=0; ihit<m_hcalHitsCol.size(); ihit++){
    m_HcalHit_x.push_back( m_hcalHitsCol[ihit]->getPosition().x() );
    m_HcalHit_y.push_back( m_hcalHitsCol[ihit]->getPosition().y() );
    m_HcalHit_z.push_back( m_hcalHitsCol[ihit]->getPosition().z() );
    m_HcalHit_E.push_back( m_hcalHitsCol[ihit]->getEnergy() );
    m_HcalHit_layer.push_back( m_hcalHitsCol[ihit]->getLayer() );
    m_totE_HcalSim += m_hcalHitsCol[ihit]->getEnergy(); 
  }
  t_SimBar->Fill();

  //----ECAL HalfClusters
  ClearHalfCluster();
  m_totE_HFClusU = 0.;
  m_totE_HFClusV = 0.;
  std::vector<PandoraPlus::CaloHalfCluster*> m_halfclusterV; m_halfclusterV.clear();
  std::vector<PandoraPlus::CaloHalfCluster*> m_halfclusterU; m_halfclusterU.clear();
  for(int i=0; i<m_DataCol.map_HalfCluster["TruthESClusterU"].size(); i++){
    m_halfclusterU.push_back( m_DataCol.map_HalfCluster["TruthESClusterU"][i].get() );
  }
  for(int i=0; i<m_DataCol.map_HalfCluster["TruthESClusterV"].size(); i++){
    m_halfclusterV.push_back( m_DataCol.map_HalfCluster["TruthESClusterV"][i].get() );
  }
  for(int i=0; i<m_halfclusterV.size(); i++){  // loop half cluster V
    m_HalfClusterV_x.push_back(m_halfclusterV[i]->getPos().x());
    m_HalfClusterV_y.push_back(m_halfclusterV[i]->getPos().x());
    m_HalfClusterV_z.push_back(m_halfclusterV[i]->getPos().x());
    m_HalfClusterV_E.push_back(m_halfclusterV[i]->getEnergy());    
    m_totE_HFClusV += m_halfclusterV[i]->getEnergy();
  }  
  for(int i=0; i<m_halfclusterU.size(); i++){  // loop half cluster U
    m_HalfClusterU_x.push_back(m_halfclusterU[i]->getPos().x());
    m_HalfClusterU_y.push_back(m_halfclusterU[i]->getPos().x());
    m_HalfClusterU_z.push_back(m_halfclusterU[i]->getPos().x());
    m_HalfClusterU_E.push_back(m_halfclusterU[i]->getEnergy());
    m_totE_HFClusU += m_halfclusterU[i]->getEnergy();
  }
  t_HalfCluster->Fill();


  //Save Cluster and MCP info
  ClearCluster();
  std::vector<PandoraPlus::Calo3DCluster*> m_clusvec;
  for(int i=0; i<m_DataCol.map_CaloCluster["TruthMergedEcalCluster"].size(); i++)
    m_clusvec.push_back( m_DataCol.map_CaloCluster["TruthMergedEcalCluster"][i].get() );

  m_Nclus = m_clusvec.size();
  m_totE_EcalRec = 0;
  for(int ic=0; ic<m_clusvec.size(); ic++){
    //For ECAL cluster
    m_Clus_x.push_back( m_clusvec[ic]->getShowerCenter().x() );
    m_Clus_y.push_back( m_clusvec[ic]->getShowerCenter().y() );
    m_Clus_z.push_back( m_clusvec[ic]->getShowerCenter().z() );
    m_Clus_E.push_back( m_clusvec[ic]->getLongiE() );
    m_Clus_Px.push_back( m_clusvec[ic]->getAxis().x() );
    m_Clus_Py.push_back( m_clusvec[ic]->getAxis().y() );
    m_Clus_Pz.push_back( m_clusvec[ic]->getAxis().z() );
    m_Clus_Nhit.push_back( m_clusvec[ic]->getCluster().size() );
    m_Clus_Ntrk.push_back( m_clusvec[ic]->getAssociatedTracks().size() );
    if( m_clusvec[ic]->getAssociatedTracks().size()>=1 ) m_Clus_Ptrk.push_back( m_clusvec[ic]->getAssociatedTracks()[0]->getMomentum() );
    else m_Clus_Ptrk.push_back(-99.);
    float maxFrac=-99;
    int maxFracPDG=0;
    for(auto imcp:m_clusvec[ic]->getLinkedMCP()){
      if(imcp.second>maxFrac){
        maxFrac = imcp.second;
        maxFracPDG = imcp.first.getPDG();
      }
    }
    m_Clus_truthFrac.push_back(maxFrac);
    m_Clus_truthPDG.push_back(maxFracPDG);

    m_totE_EcalRec += m_clusvec[ic]->getLongiE();

    for(int ihit=0; ihit<m_clusvec[ic]->getCluster().size(); ihit++){
      m_Clus_hitx.push_back( m_clusvec[ic]->getCluster()[ihit]->getPos().x() );
      m_Clus_hity.push_back( m_clusvec[ic]->getCluster()[ihit]->getPos().y() );
      m_Clus_hitz.push_back( m_clusvec[ic]->getCluster()[ihit]->getPos().z() );
      m_Clus_hitE.push_back( m_clusvec[ic]->getCluster()[ihit]->getEnergy() );
      m_Clus_hittag.push_back( ic );
      m_Clus_hittag_trk.push_back( m_clusvec[ic]->getAssociatedTracks().size() );
    }

    //Cluster properties for photon ID
    const CaloHalfCluster* p_HFClusterU = m_clusvec[ic]->getHalfClusterUCol("LinkedLongiCluster")[0];
    const CaloHalfCluster* p_HFClusterV = m_clusvec[ic]->getHalfClusterVCol("LinkedLongiCluster")[0];
    int startLayer = min(p_HFClusterU->getBeginningDlayer(), p_HFClusterV->getBeginningDlayer());
    int endLayer = max(p_HFClusterU->getEndDlayer(), p_HFClusterV->getEndDlayer());
    int maxELayer = -99;
    int maxWidthLayer = -99;
    double maxE = -99;
    double maxWidth = -99; 
    double maxScndM = -99;
    double E1=0;
    double E2=0; 
    double E5=0;
    double Ehalf=0;

    int icount = 0;
    for(int il=startLayer; il<endLayer; il++){
      std::vector<const Calo1DCluster*> m_showersU = p_HFClusterU->getClusterInLayer(il);
      std::vector<const Calo1DCluster*> m_showersV = p_HFClusterV->getClusterInLayer(il);
      if(m_showersU.size()==0 && m_showersV.size()==0) continue;

      double showerE = 0;
      double widthU = 0;
      double widthV = 0;
      double scndMomentU = 0;
      double scndMomentV = 0;
      for(auto ish:m_showersU){
        showerE += ish->getEnergy();
        widthU += ish->getWidth();
        scndMomentU += ish->getScndMoment();
      }
      for(auto ish:m_showersV){
        showerE += ish->getEnergy();
        widthV += ish->getWidth();
        scndMomentV += ish->getScndMoment();
      } 
      double width = sqrt(widthU*widthU + widthV*widthV);
      double scndM = sqrt(scndMomentU*scndMomentU + scndMomentV*scndMomentV);
      if(showerE>maxE){ maxE = showerE; maxELayer = il; } 
      if(width>maxWidth){ maxWidth = width; maxWidthLayer = il; }
      if(scndM>maxScndM) maxScndM = scndM;

      if(icount==0) E1 = showerE;
      if(icount<2) E2 += showerE;
      if(icount<5) E5 += showerE;
      if(il<(startLayer+endLayer)/2) Ehalf += showerE;
      icount++;
    }
    m_Clus_startLayer.push_back(startLayer);
    m_Clus_endLayer.push_back(endLayer);
    m_Clus_maxELayer.push_back(maxELayer);
    m_Clus_maxWidthLayer.push_back(maxWidthLayer);
    m_Clus_typeU.push_back(p_HFClusterU->getType());
    m_Clus_typeV.push_back(p_HFClusterV->getType());
    m_Clus_width.push_back(maxWidth);
    m_Clus_ScndM.push_back(maxScndM);
    m_Clus_E1Etot.push_back(E1/m_clusvec[ic]->getLongiE());
    m_Clus_E2Etot.push_back(E2/m_clusvec[ic]->getLongiE());
    m_Clus_E5Etot.push_back(E5/m_clusvec[ic]->getLongiE()); 
    m_Clus_EhalfEtot.push_back(Ehalf/m_clusvec[ic]->getLongiE());

  }


  //For MCParticle
  std::vector<edm4hep::MCParticle> m_MCPCol = m_DataCol.collectionMap_MC[name_MCParticleCol.value()];  
  m_Nmc = m_MCPCol.size(); 
  for(int imc=0; imc<m_MCPCol.size(); imc++){
    m_mcPdgid.push_back( m_MCPCol[imc].getPDG() );
    m_mcStatus.push_back( m_MCPCol[imc].getGeneratorStatus() );
    m_mcPx.push_back( m_MCPCol[imc].getMomentum()[0] );
    m_mcPy.push_back( m_MCPCol[imc].getMomentum()[1] );
    m_mcPz.push_back( m_MCPCol[imc].getMomentum()[2] );
    m_mcEn.push_back( m_MCPCol[imc].getEnergy() );
    m_mcMass.push_back( m_MCPCol[imc].getMass() );
    m_mcCharge.push_back( m_MCPCol[imc].getCharge() );
    m_mcEPx.push_back( m_MCPCol[imc].getEndpoint()[0] );
    m_mcEPy.push_back( m_MCPCol[imc].getEndpoint()[1] );
    m_mcEPz.push_back( m_MCPCol[imc].getEndpoint()[2] );
  }


  //For HCAL cluster
  std::vector<PandoraPlus::Calo3DCluster*> m_hcalclusvec;
  for(int i=0; i<m_DataCol.map_CaloCluster["TruthMergedHcalCluster"].size(); i++)
    m_hcalclusvec.push_back( m_DataCol.map_CaloCluster["TruthMergedHcalCluster"][i].get() );
 
  m_totE_HcalRec = 0.; 
  for(int ic=0; ic<m_hcalclusvec.size(); ic++)
  {
    m_Hcal_clus_x.push_back( m_hcalclusvec[ic]->getHitCenter().x() );
    m_Hcal_clus_y.push_back( m_hcalclusvec[ic]->getHitCenter().y() );
    m_Hcal_clus_z.push_back( m_hcalclusvec[ic]->getHitCenter().z() );
    m_Hcal_clus_E.push_back( m_hcalclusvec[ic]->getHitsE() );
    m_totE_HcalRec += m_hcalclusvec[ic]->getHitsE();
    for(int ihit=0; ihit<m_hcalclusvec[ic]->getCaloHits().size(); ihit++)
    {
      m_Hcal_hit_x.push_back( m_hcalclusvec[ic]->getCaloHits()[ihit]->getPosition().x() );
      m_Hcal_hit_y.push_back( m_hcalclusvec[ic]->getCaloHits()[ihit]->getPosition().y() );
      m_Hcal_hit_z.push_back( m_hcalclusvec[ic]->getCaloHits()[ihit]->getPosition().z() );
      m_Hcal_hit_E.push_back( m_hcalclusvec[ic]->getCaloHits()[ihit]->getEnergy() );
      m_Hcal_hit_tag.push_back( ic );
    }
  }
  t_Cluster->Fill();

  // ---------------------------------------------
  // Save Track info
  ClearTrack();
  std::vector<PandoraPlus::Track*> m_trkCol; 
  for(int it=0; it<m_DataCol.TrackCol.size(); it++)
    m_trkCol.push_back( m_DataCol.TrackCol[it].get() );

  m_Ntrk = m_trkCol.size();
  for(int itrk=0; itrk<m_Ntrk; itrk++){
    m_type.push_back(m_trkCol[itrk]->getType());
    std::vector<TrackState> AllTrackStates = m_trkCol[itrk]->getAllTrackStates();
    for(int istate=0; istate<AllTrackStates.size(); istate++){
      m_trkstate_d0.push_back( AllTrackStates[istate].D0 );
      m_trkstate_z0.push_back( AllTrackStates[istate].Z0 );
      m_trkstate_phi.push_back( AllTrackStates[istate].phi0 );
      m_trkstate_tanL.push_back( AllTrackStates[istate].tanLambda );
      m_trkstate_kappa.push_back( AllTrackStates[istate].Kappa);
      m_trkstate_omega.push_back( AllTrackStates[istate].Omega );
      m_trkstate_refx.push_back( AllTrackStates[istate].referencePoint.X() );
      m_trkstate_refy.push_back( AllTrackStates[istate].referencePoint.Y() );
      m_trkstate_refz.push_back( AllTrackStates[istate].referencePoint.Z() );
      m_trkstate_location.push_back( AllTrackStates[istate].location );
      m_trkstate_tag.push_back(itrk);
  }}
  t_Track->Fill();


  //Save PFO info
  ClearPFO();
  std::vector<PandoraPlus::PFObject*> m_pfobjects; m_pfobjects.clear();
  for(int ip=0; ip<m_DataCol.map_PFObjects["TruthCombPFO"].size(); ip++)
    m_pfobjects.push_back(m_DataCol.map_PFObjects["TruthCombPFO"][ip].get());

  for(int ip=0; ip<m_pfobjects.size(); ip++){
    std::vector<const Track*> t_tracks = m_pfobjects[ip]->getTracks();
    std::vector<const Calo3DCluster*> t_ecal_clusters = m_pfobjects[ip]->getECALClusters();
    std::vector<const Calo3DCluster*> t_hcal_clusters =  m_pfobjects[ip]->getHCALClusters();

    pfo_tag.push_back(ip);
    n_track.push_back(t_tracks.size());
    n_ecal_clus.push_back(t_ecal_clusters.size());
    n_hcal_clus.push_back(t_hcal_clusters.size());

    for(int it=0; it<t_tracks.size(); it++){
      m_trk_pfo_tag.push_back(ip);
      m_trk_pt.push_back(t_tracks[it]->getPt());
      m_trk_pz.push_back(t_tracks[it]->getPz());
      auto mcp = t_tracks[it]->getLeadingMCP();
      m_trk_mcpid.push_back(mcp.getPDG());
      m_trk_mc_px.push_back(mcp.getMomentum().x);
      m_trk_mc_py.push_back(mcp.getMomentum().y);
      m_trk_mc_pz.push_back(mcp.getMomentum().z);
      m_trk_mc_E.push_back(mcp.getEnergy());
    }

    for(int ie=0; ie<t_ecal_clusters.size(); ie++){
      m_ecal_pfo_tag.push_back(ip);
      m_ecal_clus_x.push_back(t_ecal_clusters[ie]->getShowerCenter().x());
      m_ecal_clus_y.push_back(t_ecal_clusters[ie]->getShowerCenter().y());
      m_ecal_clus_z.push_back(t_ecal_clusters[ie]->getShowerCenter().z());
      m_ecal_clus_E.push_back(t_ecal_clusters[ie]->getLongiE());
      auto mcp = t_ecal_clusters[ie]->getLeadingMCP();
      m_ecal_clus_mcpid.push_back(mcp.getPDG());
      m_ecal_clus_mc_px.push_back(mcp.getMomentum().x);
      m_ecal_clus_mc_py.push_back(mcp.getMomentum().y);
      m_ecal_clus_mc_pz.push_back(mcp.getMomentum().z);
      m_ecal_clus_mc_E.push_back(mcp.getEnergy());
    }
    for(int ih=0; ih<t_hcal_clusters.size(); ih++){
      m_hcal_pfo_tag.push_back(ip);
      m_hcal_clus_x.push_back(t_hcal_clusters[ih]->getHitCenter().x());
      m_hcal_clus_y.push_back(t_hcal_clusters[ih]->getHitCenter().y());
      m_hcal_clus_z.push_back(t_hcal_clusters[ih]->getHitCenter().z());
      m_hcal_clus_E.push_back(t_hcal_clusters[ih]->getHitsE());
      auto mcp = t_hcal_clusters[ih]->getLeadingMCP();
      m_hcal_clus_mcpid.push_back(mcp.getPDG());
      m_hcal_clus_mc_px.push_back(mcp.getMomentum().x);
      m_hcal_clus_mc_py.push_back(mcp.getMomentum().y);
      m_hcal_clus_mc_pz.push_back(mcp.getMomentum().z);
      m_hcal_clus_mc_E.push_back(mcp.getEnergy());
    }

  }
  t_PFO->Fill();


  //Clean Events
  //system("/cefs/higgs/songwz/winter22/CEPCSW/workarea/memory/memory_test.sh before_clean");
  m_DataCol.Clear();
  //system("/cefs/higgs/songwz/winter22/CEPCSW/workarea/memory/memory_test.sh event_end");


  std::cout<<"Event: "<<_nEvt<<" is done"<<std::endl;
  _nEvt ++ ;
  return StatusCode::SUCCESS;
}

StatusCode PandoraPlusPFAlg::finalize()
{
  m_wfile->cd();
  t_SimBar->Write();
  t_HalfCluster->Write();
  t_Cluster->Write();
  t_Track->Write();
  t_PFO->Write();
  m_wfile->Close();
  delete m_wfile, t_SimBar, t_HalfCluster, t_Cluster, t_Track, t_PFO;

  delete m_pMCParticleCreator;
  delete m_pTrackCreator; 
  delete m_pCaloHitsCreator;
  delete m_pOutputCreator;

  delete r_MCParticleCol;
  for(auto iter : r_TrackCols) delete iter; 
  //for(auto iter : r_ECalHitCols) delete iter; 
  //for(auto iter : r_HCalHitCols) delete iter; 
  for(auto iter : r_CaloHitCols) delete iter; 
  for(auto iter : map_CaloMCPAssoCols) delete iter.second;
  r_TrackCols.clear();
  //r_ECalHitCols.clear();
  //r_HCalHitCols.clear();
  r_CaloHitCols.clear();


  info() << "Processed " << _nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}

void PandoraPlusPFAlg::ClearBar(){
  m_totE_EcalSim = -99;
  m_totE_HcalSim = -99;
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
  m_simBar_bar.clear();

  m_HcalHit_x.clear();
  m_HcalHit_y.clear();
  m_HcalHit_z.clear();
  m_HcalHit_E.clear();
  m_HcalHit_layer.clear();
}


void PandoraPlusPFAlg::ClearHalfCluster(){
  m_totE_HFClusU = -99.;
  m_totE_HFClusV = -99.;
  m_HalfClusterV_x.clear();
  m_HalfClusterV_y.clear();
  m_HalfClusterV_z.clear();
  m_HalfClusterV_E.clear();
  m_HalfClusterU_x.clear();
  m_HalfClusterU_y.clear();
  m_HalfClusterU_z.clear();
  m_HalfClusterU_E.clear();
}


void PandoraPlusPFAlg::ClearCluster(){
  m_Nclus = -99;
  m_Nmc = -99;
  m_totE_EcalRec = -99;
  m_totE_HcalRec = -99;
  m_Clus_x.clear();
  m_Clus_y.clear();
  m_Clus_z.clear();
  m_Clus_E.clear();
  m_Clus_Px.clear();
  m_Clus_Py.clear();
  m_Clus_Pz.clear();
  m_Clus_Ptrk.clear();
  m_Clus_Ntrk.clear();
  m_Clus_Nhit.clear();
  m_Clus_truthFrac.clear();
  m_Clus_truthPDG.clear();
  m_Clus_startLayer.clear(); 
  m_Clus_endLayer.clear(); 
  m_Clus_maxELayer.clear();
  m_Clus_maxWidthLayer.clear();
  m_Clus_typeU.clear();
  m_Clus_typeV.clear();
  m_Clus_width.clear();
  m_Clus_ScndM.clear();
  m_Clus_E1Etot.clear();
  m_Clus_E2Etot.clear();
  m_Clus_E5Etot.clear();
  m_Clus_EhalfEtot.clear();
  m_Clus_EaxisEtot.clear();
  m_Clus_hitx.clear();
  m_Clus_hity.clear();
  m_Clus_hitz.clear();
  m_Clus_hitE.clear();
  m_Clus_hittag.clear();
  m_Clus_hittag_trk.clear();
  m_mcPdgid.clear();
  m_mcStatus.clear();
  m_mcPx.clear();
  m_mcPy.clear();
  m_mcPz.clear();
  m_mcEn.clear();
  m_mcMass.clear();
  m_mcCharge.clear();
  m_mcEPx.clear();
  m_mcEPy.clear();
  m_mcEPz.clear();

  m_Hcal_clus_x.clear();
  m_Hcal_clus_y.clear();
  m_Hcal_clus_z.clear();
  m_Hcal_clus_E.clear();
  m_Hcal_hit_tag.clear();
  m_Hcal_hit_x.clear();
  m_Hcal_hit_y.clear();
  m_Hcal_hit_z.clear();
  m_Hcal_hit_E.clear();
}


void PandoraPlusPFAlg::ClearTrack(){
  m_Ntrk=-99;
  m_type.clear();
  m_trkstate_d0.clear(); 
  m_trkstate_z0.clear();
  m_trkstate_phi.clear();
  m_trkstate_tanL.clear();
  m_trkstate_kappa.clear();
  m_trkstate_omega.clear();
  m_trkstate_refx.clear();
  m_trkstate_refy.clear();
  m_trkstate_refz.clear();
  m_trkstate_location.clear();
  m_trkstate_tag.clear();
}


void PandoraPlusPFAlg::ClearPFO(){
  pfo_tag.clear();
  n_track.clear();
  n_ecal_clus.clear();
  n_hcal_clus.clear();
  m_trk_pfo_tag.clear();
  m_trk_mcpid.clear();
  m_trk_pt.clear();
  m_trk_pz.clear();
  m_trk_mc_px.clear();
  m_trk_mc_py.clear();
  m_trk_mc_pz.clear();
  m_trk_mc_E.clear();
  m_ecal_pfo_tag.clear();
  m_ecal_clus_x.clear();
  m_ecal_clus_y.clear();
  m_ecal_clus_z.clear();
  m_ecal_clus_E.clear();
  m_ecal_clus_mcpid.clear();
  m_ecal_clus_mc_px.clear();
  m_ecal_clus_mc_py.clear();
  m_ecal_clus_mc_pz.clear();
  m_ecal_clus_mc_E.clear();
  m_hcal_pfo_tag.clear();
  m_hcal_clus_x.clear();
  m_hcal_clus_y.clear();
  m_hcal_clus_z.clear();
  m_hcal_clus_E.clear();
  m_hcal_clus_mcpid.clear();
  m_hcal_clus_mc_px.clear();
  m_hcal_clus_mc_py.clear();
  m_hcal_clus_mc_pz.clear();
  m_hcal_clus_mc_E.clear();
}


#endif
