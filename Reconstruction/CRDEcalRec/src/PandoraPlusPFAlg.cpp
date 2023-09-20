#ifndef PANDORAPLUS_ALG_C
#define PANDORAPLUS_ALG_C

#include "PandoraPlusPFAlg.h"

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
  if(!name_MCPRecoCaloAssoCol.empty()) r_MCPRecoCaloAssoCol = new DataHandle<edm4hep::MCRecoCaloParticleAssociationCollection> (name_MCPRecoCaloAssoCol, Gaudi::DataHandle::Reader, this);


  //Register Algorithms
  //--- Initialize algorithm maps ---
  m_algorithmManager.RegisterAlgorithmFactory("ExampleAlg",             new ExampleAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("GlobalClusteringAlg",    new GlobalClusteringAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("LocalMaxFindingAlg",     new LocalMaxFindingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("HoughClusteringAlg",     new HoughClusteringAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("TrackMatchingAlg",       new TrackMatchingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("ConeClustering2DAlg",    new ConeClustering2DAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("AxisMergingAlg",         new AxisMergingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("EnergySplittingAlg",     new EnergySplittingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("EnergyTimeMatchingAlg",  new EnergyTimeMatchingAlg::Factory);
  //m_algorithmManager.RegisterAlgorithmFactory("ConeClusteringAlg",      new ConeClusteringAlg::Factory);
  //m_algorithmManager.RegisterAlgorithmFactory("ConeClusteringAlgHCAL",  new ConeClusteringAlg::Factory);


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
    t_Layers = new TTree("RecLayers","RecLayers");
    t_Hough = new TTree("Hough","Hough");   // yyy: Hough
    t_Tower = new TTree("RecTowers", "RecTowers");
    t_Cluster = new TTree("RecClusters", "RecClusters");
    t_Track = new TTree("RecTracks", "RecTracks");
    t_axis = new TTree("Axis","Axis"); 

    //Bar
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
    t_SimBar->Branch("HcalHit_x", &m_HcalHit_x); 
    t_SimBar->Branch("HcalHit_y", &m_HcalHit_y); 
    t_SimBar->Branch("HcalHit_z", &m_HcalHit_z); 
    t_SimBar->Branch("HcalHit_E", &m_HcalHit_E); 
    t_SimBar->Branch("HcalHit_layer", &m_HcalHit_layer); 

  
    //BarShowers
    t_Layers->Branch("NshowerU", &m_NshowerU);
    t_Layers->Branch("NshowerV", &m_NshowerV);
    t_Layers->Branch("barShowerU_tag", &m_barShowerU_tag);
    t_Layers->Branch("barShowerU_x", &m_barShowerU_x);
    t_Layers->Branch("barShowerU_y", &m_barShowerU_y);
    t_Layers->Branch("barShowerU_z", &m_barShowerU_z);
    t_Layers->Branch("barShowerU_E", &m_barShowerU_E);
    t_Layers->Branch("barShowerU_T", &m_barShowerU_T);
    //t_Layers->Branch("barShowerU_module", &m_barShowerU_module);
    //t_Layers->Branch("barShowerU_stave", &m_barShowerU_stave);
    //t_Layers->Branch("barShowerU_part", &m_barShowerU_part);
    //t_Layers->Branch("barShowerU_dlayer", &m_barShowerU_dlayer);
    //t_Layers->Branch("barShowerU_slayer", &m_barShowerU_slayer);
    //t_Layers->Branch("barShowerU_bar", &m_barShowerU_bar);
    t_Layers->Branch("barShowerV_tag", &m_barShowerV_tag);
    t_Layers->Branch("barShowerV_x", &m_barShowerV_x);
    t_Layers->Branch("barShowerV_y", &m_barShowerV_y);
    t_Layers->Branch("barShowerV_z", &m_barShowerV_z);
    t_Layers->Branch("barShowerV_E", &m_barShowerV_E);
    t_Layers->Branch("barShowerV_T", &m_barShowerV_T);
    //t_Layers->Branch("barShowerV_module", &m_barShowerV_module);
    //t_Layers->Branch("barShowerV_stave", &m_barShowerV_stave);
    //t_Layers->Branch("barShowerV_part", &m_barShowerV_part);
    //t_Layers->Branch("barShowerV_dlayer", &m_barShowerV_dlayer);
    //t_Layers->Branch("barShowerV_slayer", &m_barShowerV_slayer);
    //t_Layers->Branch("barShowerV_bar", &m_barShowerV_bar);


    // yyy: Hough cluster
    t_Hough->Branch("NHoughAxisV", &m_NHoughAxisV);
    t_Hough->Branch("halfclusV_tag", &m_halfclusV_tag);
    t_Hough->Branch("houghV_x", &m_houghV_x);
    t_Hough->Branch("houghV_y", &m_houghV_y);
    t_Hough->Branch("houghV_z", &m_houghV_z);
    t_Hough->Branch("houghV_E", &m_houghV_E);
    //t_Hough->Branch("houghV_module", &m_houghV_module);
    //t_Hough->Branch("houghV_part", &m_houghV_part);
    //t_Hough->Branch("houghV_stave", &m_houghV_stave);
    //t_Hough->Branch("houghV_dlayer", &m_houghV_dlayer);
    //t_Hough->Branch("houghV_slayer", &m_houghV_slayer);
    t_Hough->Branch("NHoughAxisU", &m_NHoughAxisU);
    t_Hough->Branch("halfclusU_tag", &m_halfclusU_tag);
    t_Hough->Branch("houghU_x", &m_houghU_x);
    t_Hough->Branch("houghU_y", &m_houghU_y);
    t_Hough->Branch("houghU_z", &m_houghU_z);
    t_Hough->Branch("houghU_E", &m_houghU_E);
    //t_Hough->Branch("houghU_module", &m_houghU_module);
    //t_Hough->Branch("houghU_part", &m_houghU_part);
    //t_Hough->Branch("houghU_stave", &m_houghU_stave);
    //t_Hough->Branch("houghU_dlayer", &m_houghU_dlayer);
    //t_Hough->Branch("houghU_slayer", &m_houghU_slayer);


    //Clusters and MCParticles
    t_Cluster->Branch("Nclus", &m_Nclus);
    t_Cluster->Branch("totE",  &m_totE);
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


    // Axis
    t_axis->Branch("NaxisU", &m_NaxisU);
    t_axis->Branch("NaxisV", &m_NaxisV);
    t_axis->Branch("axisU_x", &m_axisU_x);
    t_axis->Branch("axisU_y", &m_axisU_y);
    t_axis->Branch("axisU_z", &m_axisU_z);
    t_axis->Branch("axisU_E", &m_axisU_E);
    t_axis->Branch("axisU_Nhit", &m_axisU_Nhit);
    t_axis->Branch("axisU_type", &m_axisU_type);
    t_axis->Branch("axisU_Nmcp", &m_axisU_Nmcp);
    t_axis->Branch("axisU_truthFracMax", &m_axisU_truthFracMax);
    t_axis->Branch("axisV_x", &m_axisV_x);
    t_axis->Branch("axisV_y", &m_axisV_y);
    t_axis->Branch("axisV_z", &m_axisV_z);
    t_axis->Branch("axisV_E", &m_axisV_E);
    t_axis->Branch("axisV_Nhit", &m_axisV_Nhit);
    t_axis->Branch("axisV_type", &m_axisV_type);
    t_axis->Branch("axisV_Nmcp", &m_axisV_Nmcp);
    t_axis->Branch("axisV_truthFracMax", &m_axisV_truthFracMax);
    t_axis->Branch("axisUhit_tag", &m_axisUhit_tag);
    t_axis->Branch("axisUhit_type", &m_axisUhit_type);
    t_axis->Branch("axisUhit_x", &m_axisUhit_x);
    t_axis->Branch("axisUhit_y", &m_axisUhit_y);
    t_axis->Branch("axisUhit_z", &m_axisUhit_z);
    t_axis->Branch("axisUhit_E", &m_axisUhit_E);
    t_axis->Branch("axisVhit_tag", &m_axisVhit_tag);
    t_axis->Branch("axisVhit_type", &m_axisVhit_type);
    t_axis->Branch("axisVhit_x", &m_axisVhit_x);
    t_axis->Branch("axisVhit_y", &m_axisVhit_y);
    t_axis->Branch("axisVhit_z", &m_axisVhit_z);
    t_axis->Branch("axisVhit_E", &m_axisVhit_E);

    //Tower
    t_Tower->Branch("module", &m_module );
    t_Tower->Branch("part", &m_part );
    t_Tower->Branch("stave", &m_stave );
    t_Tower->Branch("HFClusU_Nhit", &m_HFClusU_Nhit );
    t_Tower->Branch("HFClusU_type", &m_HFClusU_type );
    t_Tower->Branch("HFClusU_x",    &m_HFClusU_x );
    t_Tower->Branch("HFClusU_y",    &m_HFClusU_y );
    t_Tower->Branch("HFClusU_z",    &m_HFClusU_z );
    t_Tower->Branch("HFClusU_E",    &m_HFClusU_E );
    t_Tower->Branch("HFClusV_Nhit", &m_HFClusV_Nhit );
    t_Tower->Branch("HFClusV_type", &m_HFClusV_type );
    t_Tower->Branch("HFClusV_x",    &m_HFClusV_x );
    t_Tower->Branch("HFClusV_y",    &m_HFClusV_y );
    t_Tower->Branch("HFClusV_z",    &m_HFClusV_z );
    t_Tower->Branch("HFClusV_E",    &m_HFClusV_E );
    t_Tower->Branch("HFClusUhit_type", &m_HFClusUhit_type );
    t_Tower->Branch("HFClusUhit_tag",  &m_HFClusUhit_tag );
    t_Tower->Branch("HFClusUhit_x",    &m_HFClusUhit_x );
    t_Tower->Branch("HFClusUhit_y",    &m_HFClusUhit_y );
    t_Tower->Branch("HFClusUhit_z",    &m_HFClusUhit_z );
    t_Tower->Branch("HFClusUhit_E",    &m_HFClusUhit_E );
    t_Tower->Branch("HFClusVhit_type", &m_HFClusVhit_type );
    t_Tower->Branch("HFClusVhit_tag",  &m_HFClusVhit_tag );
    t_Tower->Branch("HFClusVhit_x",    &m_HFClusVhit_x );
    t_Tower->Branch("HFClusVhit_y",    &m_HFClusVhit_y );
    t_Tower->Branch("HFClusVhit_z",    &m_HFClusVhit_z );
    t_Tower->Branch("HFClusVhit_E",    &m_HFClusVhit_E );

  }

  return GaudiAlgorithm::initialize();
}

StatusCode PandoraPlusPFAlg::execute()
{
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
  m_pTrackCreator->CreateTracks( m_DataCol, r_TrackCols );
cout<<"Readin CaloHits"<<endl;
  m_pCaloHitsCreator->CreateCaloHits( m_DataCol, r_CaloHitCols, map_readout_decoder, r_MCPRecoCaloAssoCol );

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



cout<<"Write tuples"<<endl;
  //---------------------Write Ana tuples-------------------------
  //Save Raw bars information
  ClearBar();

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
  }

  std::vector<PandoraPlus::CaloHit*> m_hcalHitsCol; m_hcalHitsCol.clear();
  for(int ih=0; ih<m_DataCol.map_CaloHit["HCALBarrel"].size(); ih++)
    m_hcalHitsCol.push_back( m_DataCol.map_CaloHit["HCALBarrel"][ih].get() );

  for(int ihit=0; ihit<m_hcalHitsCol.size(); ihit++){
    m_HcalHit_x.push_back( m_hcalHitsCol[ihit]->getPosition().x() );
    m_HcalHit_y.push_back( m_hcalHitsCol[ihit]->getPosition().y() );
    m_HcalHit_z.push_back( m_hcalHitsCol[ihit]->getPosition().z() );
    m_HcalHit_E.push_back( m_hcalHitsCol[ihit]->getEnergy() );
    m_HcalHit_layer.push_back( m_hcalHitsCol[ihit]->getLayer() );
  }

  t_SimBar->Fill();
  // ------------------------------------

  //Save Layer info (local Max)
  ClearLayer();
  std::vector<PandoraPlus::CaloHalfCluster*> m_halfclusters; m_halfclusters.clear();
  for(int i=0; i<m_DataCol.map_HalfCluster["HalfClusterColU"].size(); i++)
    m_halfclusters.push_back( m_DataCol.map_HalfCluster["HalfClusterColU"][i].get() );

  for(int ic=0;ic<m_halfclusters.size(); ic++){
    std::vector<const Calo1DCluster*> tmp_shower = m_halfclusters[ic]->getLocalMaxCol("AllLocalMax");
    m_NshowerU = tmp_shower.size(); 
    for(int is=0; is<tmp_shower.size(); is++)
    {
      m_barShowerU_tag.push_back( ic );
      m_barShowerU_x.push_back( tmp_shower[is]->getPos().x() );
      m_barShowerU_y.push_back( tmp_shower[is]->getPos().y() );
      m_barShowerU_z.push_back( tmp_shower[is]->getPos().z() );
      m_barShowerU_E.push_back( tmp_shower[is]->getEnergy() );
      m_barShowerU_T.push_back( (tmp_shower[is]->getT1()+tmp_shower[is]->getT2())/2. );
      //m_barShowerU_module.push_back( tmp_shower[is]->getSeeds().at(iseed)->getModule() );
      //m_barShowerU_stave.push_back( tmp_shower[is]->getSeeds().at(iseed)->getStave() );
      //m_barShowerU_part.push_back( tmp_shower[is]->getSeeds().at(iseed)->getPart() );
      //m_barShowerU_dlayer.push_back( tmp_shower[is]->getSeeds().at(iseed)->getDlayer() );
      //m_barShowerU_slayer.push_back( tmp_shower[is]->getSeeds().at(iseed)->getSlayer() );
      //m_barShowerU_bar.push_back( tmp_shower[is]->getSeeds().at(iseed)->getBar() );
      
    }
  }

  m_halfclusters.clear();
  for(int i=0; i<m_DataCol.map_HalfCluster["HalfClusterColV"].size(); i++)
    m_halfclusters.push_back( m_DataCol.map_HalfCluster["HalfClusterColV"][i].get() );

  for(int ic=0;ic<m_halfclusters.size();ic++){
    std::vector<const Calo1DCluster*> tmp_shower = m_halfclusters[ic]->getLocalMaxCol("AllLocalMax");
    m_NshowerV = tmp_shower.size(); 
    for(int is=0; is<m_NshowerV; is++)
    {
      m_barShowerV_tag.push_back( ic );
      m_barShowerV_x.push_back( tmp_shower[is]->getPos().x() );
      m_barShowerV_y.push_back( tmp_shower[is]->getPos().y() );
      m_barShowerV_z.push_back( tmp_shower[is]->getPos().z() );
      m_barShowerV_E.push_back( tmp_shower[is]->getEnergy() );
      m_barShowerV_T.push_back( (tmp_shower[is]->getT1()+tmp_shower[is]->getT2())/2. );
      //m_barShowerV_module.push_back( tmp_shower[is]->getSeeds().at(iseed)->getModule() );
      //m_barShowerV_stave.push_back( tmp_shower[is]->getSeeds().at(iseed)->getStave() );
      //m_barShowerV_part.push_back( tmp_shower[is]->getSeeds().at(iseed)->getPart() );
      //m_barShowerV_dlayer.push_back( tmp_shower[is]->getSeeds().at(iseed)->getDlayer() );
      //m_barShowerV_slayer.push_back( tmp_shower[is]->getSeeds().at(iseed)->getSlayer() );
      //m_barShowerV_bar.push_back( tmp_shower[is]->getSeeds().at(iseed)->getBar() );
    }
  }
  t_Layers->Fill();


  // ------------------------------------
  // yyy: check Hough cluster
  ClearHough();
  std::vector<PandoraPlus::CaloHalfCluster*> m_halfclusterV; m_halfclusterV.clear();
  std::vector<PandoraPlus::CaloHalfCluster*> m_halfclusterU; m_halfclusterU.clear();
  std::vector<const PandoraPlus::CaloHalfCluster*> m_houghaxisU; m_houghaxisU.clear(); 
  std::vector<const PandoraPlus::CaloHalfCluster*> m_houghaxisV; m_houghaxisV.clear(); 
  for(int i=0; i<m_DataCol.map_HalfCluster["HalfClusterColU"].size(); i++){
    m_halfclusterU.push_back( m_DataCol.map_HalfCluster["HalfClusterColU"][i].get() );
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_clus = m_DataCol.map_HalfCluster["HalfClusterColU"][i].get()->getHalfClusterCol("HoughAxis"); 
    m_houghaxisU.insert( m_houghaxisU.end(), tmp_clus.begin(), tmp_clus.end() );
  }
  for(int i=0; i<m_DataCol.map_HalfCluster["HalfClusterColV"].size(); i++){
    m_halfclusterV.push_back( m_DataCol.map_HalfCluster["HalfClusterColV"][i].get() );
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_clus = m_DataCol.map_HalfCluster["HalfClusterColV"][i].get()->getHalfClusterCol("HoughAxis");
    m_houghaxisV.insert( m_houghaxisV.end(), tmp_clus.begin(), tmp_clus.end() );
  }

  m_NHoughAxisU = m_houghaxisU.size();
  m_NHoughAxisV = m_houghaxisV.size();
//  for(int ihc=0; ihc<m_houghaxisV.size(); ihc++){
//    m_houghV_x.push_back( m_houghaxisV[ihc]->getPos().x() );
//    m_houghV_y.push_back( m_houghaxisV[ihc]->getPos().y() );
//    m_houghV_z.push_back( m_houghaxisV[ihc]->getPos().z() );
//    m_houghV_E.push_back( m_houghaxisV[ihc]->getEnergy()  );
//  }
//  for(int ihc=0; ihc<m_houghaxisU.size(); ihc++){
//    m_houghU_x.push_back( m_houghaxisU[ihc]->getPos().x() );
//    m_houghU_y.push_back( m_houghaxisU[ihc]->getPos().y() );
//    m_houghU_z.push_back( m_houghaxisU[ihc]->getPos().z() );
//    m_houghU_E.push_back( m_houghaxisU[ihc]->getEnergy()  );
//  }


  for(int ilc=0; ilc<m_houghaxisV.size(); ilc++){
    for(int ilm=0; ilm<m_houghaxisV[ilc]->getCluster().size(); ilm++){
      m_halfclusV_tag.push_back( m_houghaxisV[ilc]->getEnergy() );
      m_houghV_x.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getPos().x() );
      m_houghV_y.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getPos().y() );
      m_houghV_z.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getPos().z() );
      m_houghV_E.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getEnergy() );
      //m_houghV_module.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getTowerID()[0][0] );
      //m_houghV_part.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getTowerID()[0][1] );
      //m_houghV_stave.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getTowerID()[0][2] );
      //m_houghV_dlayer.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getDlayer() );
      //m_houghV_slayer.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getSlayer() );
    }
  }

  for(int ilc=0; ilc<m_houghaxisU.size(); ilc++){
    for(int ilm=0; ilm<m_houghaxisU[ilc]->getCluster().size(); ilm++){
      m_halfclusU_tag.push_back( m_houghaxisU[ilc]->getEnergy() );
      m_houghU_x.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getPos().x() );
      m_houghU_y.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getPos().y() );
      m_houghU_z.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getPos().z() );
      m_houghU_E.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getEnergy() );
      //m_houghU_module.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getTowerID()[0][0] );
      //m_houghU_part.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getTowerID()[0][1] );
      //m_houghU_stave.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getTowerID()[0][2] );
      //m_houghU_dlayer.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getDlayer() );
      //m_houghU_slayer.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getSlayer() );
    }
  }
  t_Hough->Fill();


  // ---------------------------------------------
  // Save axis info (Splitted showers in tower)
  ClearAxis();
  std::vector<const PandoraPlus::CaloHalfCluster*> m_allaxisU; m_allaxisU.clear();
  std::vector<const PandoraPlus::CaloHalfCluster*> m_allaxisV; m_allaxisV.clear();

  std::vector<PandoraPlus::Calo3DCluster*> m_towerCol;
  int ntower = m_DataCol.map_CaloCluster["ESTower"].size();
  for(int it=0; it<ntower; it++)
    m_towerCol.push_back( m_DataCol.map_CaloCluster["ESTower"][it].get() );

/*  for(int it=0; it<m_towerCol.size(); it++){
    std::vector<const PandoraPlus::CaloHalfCluster*> m_HFClusUCol = m_towerCol.at(it)->getHalfClusterUCol("ESHalfClusterU");
    std::vector<const PandoraPlus::CaloHalfCluster*> m_HFClusVCol = m_towerCol.at(it)->getHalfClusterVCol("ESHalfClusterV");

    m_allaxisU.insert(m_allaxisU.end(), m_HFClusUCol.begin(), m_HFClusUCol.end());
    m_allaxisV.insert(m_allaxisV.end(), m_HFClusVCol.begin(), m_HFClusVCol.end());
  }

  for(int i=0; i<m_halfclusterU.size(); i++){
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_clus = m_halfclusterU[i]->getHalfClusterCol("MergedAxis");
    //std::vector<const PandoraPlus::CaloHalfCluster*> tmp_clus = m_halfclusterU[i]->getAllHalfClusterCol();
    m_allaxisU.insert( m_allaxisU.end(), tmp_clus.begin(), tmp_clus.end() );
  }
  for(int i=0; i<m_halfclusterV.size(); i++){
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_clus = m_halfclusterV[i]->getHalfClusterCol("MergedAxis");
    //std::vector<const PandoraPlus::CaloHalfCluster*> tmp_clus = m_halfclusterV[i]->getAllHalfClusterCol();
    m_allaxisV.insert( m_allaxisV.end(), tmp_clus.begin(), tmp_clus.end() );
  }
*/
  for(int i=0; i<m_DataCol.map_HalfCluster["ESHalfClusterU"].size(); i++){
    m_allaxisU.push_back( m_DataCol.map_HalfCluster["ESHalfClusterU"][i].get());
  }
  for(int i=0; i<m_DataCol.map_HalfCluster["ESHalfClusterV"].size(); i++){
    m_allaxisV.push_back( m_DataCol.map_HalfCluster["ESHalfClusterV"][i].get());
  }

  m_NaxisU = m_allaxisU.size();
  m_NaxisV = m_allaxisV.size();
  for(int iax=0; iax<m_allaxisU.size(); iax++){
    m_axisU_Nhit.push_back( m_allaxisU[iax]->getCluster().size() );
    m_axisU_type.push_back( m_allaxisU[iax]->getType() );
    m_axisU_x.push_back( m_allaxisU[iax]->getPos().x() );
    m_axisU_y.push_back( m_allaxisU[iax]->getPos().y() );
    m_axisU_z.push_back( m_allaxisU[iax]->getPos().z() );
    m_axisU_E.push_back( m_allaxisU[iax]->getEnergy() );
    for(int ihit=0; ihit<m_allaxisU[iax]->getCluster().size(); ihit++){
      m_axisUhit_tag.push_back( m_allaxisU[iax]->getEnergy() );
      m_axisUhit_type.push_back( m_allaxisU[iax]->getType() );
      m_axisUhit_x.push_back( m_allaxisU[iax]->getCluster()[ihit]->getPos().x() );
      m_axisUhit_y.push_back( m_allaxisU[iax]->getCluster()[ihit]->getPos().y() );
      m_axisUhit_z.push_back( m_allaxisU[iax]->getCluster()[ihit]->getPos().z() );
      m_axisUhit_E.push_back( m_allaxisU[iax]->getCluster()[ihit]->getEnergy() );
    }

    std::map<edm4hep::MCParticle, float> map_truthP_totE; map_truthP_totE.clear();
    for(auto ish : m_allaxisU[iax]->getCluster()){
      for(auto ibar : ish->getBars()){
        for(auto ipair : ibar->getLinkedMCP()) map_truthP_totE[ipair.first] += ibar->getEnergy()*ipair.second;
      }
    }
    float maxFrac = -99;
    int Nmcp=0;
    for(auto imcp : map_truthP_totE){
      if(imcp.second/m_allaxisU[iax]->getEnergy()>maxFrac) maxFrac=imcp.second/m_allaxisU[iax]->getEnergy();
      if(imcp.second/m_allaxisU[iax]->getEnergy()>0.05) Nmcp++;
    }
    map_truthP_totE.clear();
    m_axisU_Nmcp.push_back(Nmcp);
    m_axisU_truthFracMax.push_back(maxFrac); 
  }

  for(int iax=0; iax<m_allaxisV.size(); iax++){
    m_axisV_Nhit.push_back( m_allaxisV[iax]->getCluster().size() );
    m_axisV_type.push_back( m_allaxisV[iax]->getType() );
    m_axisV_x.push_back( m_allaxisV[iax]->getPos().x() );
    m_axisV_y.push_back( m_allaxisV[iax]->getPos().y() );
    m_axisV_z.push_back( m_allaxisV[iax]->getPos().z() );
    m_axisV_E.push_back( m_allaxisV[iax]->getEnergy() );
    for(int ihit=0; ihit<m_allaxisV[iax]->getCluster().size(); ihit++){
      m_axisVhit_tag.push_back( m_allaxisV[iax]->getEnergy() );
      m_axisVhit_type.push_back( m_allaxisV[iax]->getType() );
      m_axisVhit_x.push_back( m_allaxisV[iax]->getCluster()[ihit]->getPos().x() );
      m_axisVhit_y.push_back( m_allaxisV[iax]->getCluster()[ihit]->getPos().y() );
      m_axisVhit_z.push_back( m_allaxisV[iax]->getCluster()[ihit]->getPos().z() );
      m_axisVhit_E.push_back( m_allaxisV[iax]->getCluster()[ihit]->getEnergy() );
    }

    std::map<edm4hep::MCParticle, float> map_truthP_totE; map_truthP_totE.clear();
    for(auto ish : m_allaxisV[iax]->getCluster()){
      for(auto ibar : ish->getBars()){
        for(auto ipair : ibar->getLinkedMCP()) map_truthP_totE[ipair.first] += ibar->getEnergy()*ipair.second;
      }
    }
    float maxFrac = -99;
    int Nmcp = 0;
    for(auto imcp : map_truthP_totE){
      if(imcp.second/m_allaxisV[iax]->getEnergy()>maxFrac) maxFrac=imcp.second/m_allaxisV[iax]->getEnergy();
      if(imcp.second/m_allaxisV[iax]->getEnergy()>0.05) Nmcp++;
    }
    map_truthP_totE.clear();
    m_axisV_truthFracMax.push_back(maxFrac);
    m_axisV_Nmcp.push_back(Nmcp);
  }
  t_axis->Fill();


  // --------------------------------------------
  // Save splitted half cluster in tower
  ntower = m_DataCol.map_CaloCluster["ESTower"].size();
    ClearTower();
  for(int it=0; it<ntower; it++){

    const PandoraPlus::Calo3DCluster* p_tower = m_DataCol.map_CaloCluster["ESTower"][it].get();
    if(p_tower->getTowerID().size()<1){ cout<<"ERROR: no towerID! "<<endl; continue; }
    m_module = p_tower->getTowerID()[0][0];
    m_part = p_tower->getTowerID()[0][1];
    m_stave = p_tower->getTowerID()[0][2];

    std::vector<const PandoraPlus::CaloHalfCluster*> m_HFClusUCol = m_towerCol.at(it)->getHalfClusterUCol("ESHalfClusterU");
    std::vector<const PandoraPlus::CaloHalfCluster*> m_HFClusVCol = m_towerCol.at(it)->getHalfClusterVCol("ESHalfClusterV");
    for(int ih=0; ih<m_HFClusUCol.size(); ih++){
      m_HFClusU_Nhit.push_back( m_HFClusUCol[ih]->getCluster().size() );
      m_HFClusU_type.push_back( m_HFClusUCol[ih]->getType() );
      m_HFClusU_x.push_back( m_HFClusUCol[ih]->getPos().x() );
      m_HFClusU_y.push_back( m_HFClusUCol[ih]->getPos().y() );
      m_HFClusU_z.push_back( m_HFClusUCol[ih]->getPos().z() );
      m_HFClusU_E.push_back( m_HFClusUCol[ih]->getEnergy() );
      for(int ihit=0; ihit<m_HFClusUCol[ih]->getCluster().size(); ihit++){
        m_HFClusUhit_tag.push_back( m_HFClusUCol[ih]->getEnergy() );
        m_HFClusUhit_type.push_back( m_HFClusUCol[ih]->getType() );
        m_HFClusUhit_x.push_back( m_HFClusUCol[ih]->getCluster()[ihit]->getPos().x() );
        m_HFClusUhit_y.push_back( m_HFClusUCol[ih]->getCluster()[ihit]->getPos().y() );
        m_HFClusUhit_z.push_back( m_HFClusUCol[ih]->getCluster()[ihit]->getPos().z() );
        m_HFClusUhit_E.push_back( m_HFClusUCol[ih]->getCluster()[ihit]->getEnergy() );
      }
    }

    for(int ih=0; ih<m_HFClusVCol.size(); ih++){
      m_HFClusV_Nhit.push_back( m_HFClusVCol[ih]->getCluster().size() );
      m_HFClusV_type.push_back( m_HFClusVCol[ih]->getType() );
      m_HFClusV_x.push_back( m_HFClusVCol[ih]->getPos().x() );
      m_HFClusV_y.push_back( m_HFClusVCol[ih]->getPos().y() );
      m_HFClusV_z.push_back( m_HFClusVCol[ih]->getPos().z() );
      m_HFClusV_E.push_back( m_HFClusVCol[ih]->getEnergy() );
      for(int ihit=0; ihit<m_HFClusVCol[ih]->getCluster().size(); ihit++){
        m_HFClusVhit_tag.push_back( m_HFClusVCol[ih]->getEnergy() );
        m_HFClusVhit_type.push_back( m_HFClusVCol[ih]->getType() );
        m_HFClusVhit_x.push_back( m_HFClusVCol[ih]->getCluster()[ihit]->getPos().x() );
        m_HFClusVhit_y.push_back( m_HFClusVCol[ih]->getCluster()[ihit]->getPos().y() );
        m_HFClusVhit_z.push_back( m_HFClusVCol[ih]->getCluster()[ihit]->getPos().z() );
        m_HFClusVhit_E.push_back( m_HFClusVCol[ih]->getCluster()[ihit]->getEnergy() );
      }
    }

  }
    t_Tower->Fill();

cout<<"  Write 3DClusters and MCP "<<endl;
  // ---------------------------------------------
  //Save Cluster and MCP info
  ClearCluster();
  std::vector<PandoraPlus::Calo3DCluster*> m_clusvec;
  //std::vector<PandoraPlus::Calo3DCluster*> m_clusvec = m_DataCol.map_CaloCluster["HCALCluster"];
  for(int i=0; i<m_DataCol.map_CaloCluster["EcalCluster"].size(); i++)
    m_clusvec.push_back( m_DataCol.map_CaloCluster["EcalCluster"][i].get() );

  m_Nclus = m_clusvec.size();
  m_totE = 0;
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

    m_totE += m_clusvec[ic]->getLongiE();

    for(int ihit=0; ihit<m_clusvec[ic]->getCluster().size(); ihit++){
      m_Clus_hitx.push_back( m_clusvec[ic]->getCluster()[ihit]->getPos().x() );
      m_Clus_hity.push_back( m_clusvec[ic]->getCluster()[ihit]->getPos().y() );
      m_Clus_hitz.push_back( m_clusvec[ic]->getCluster()[ihit]->getPos().z() );
      m_Clus_hitE.push_back( m_clusvec[ic]->getCluster()[ihit]->getEnergy() );
      m_Clus_hittag.push_back( m_clusvec[ic]->getEnergy() );
      m_Clus_hittag_trk.push_back( m_clusvec[ic]->getAssociatedTracks().size() );
    }

    //For HCAL cluster
//    m_Clus_x.push_back( m_clusvec[ic]->getHitCenter().x() );
//  m_Clus_y.push_back( m_clusvec[ic]->getHitCenter().y() );
//  m_Clus_z.push_back( m_clusvec[ic]->getHitCenter().z() );
//  m_Clus_E.push_back( m_clusvec[ic]->getHitsE() );
//  m_Clus_Nhit.push_back( m_clusvec[ic]->getCaloHits().size() );
//  for(int ihit=0; ihit<m_clusvec[ic]->getCaloHits().size(); ihit++){
//    m_Clus_hitx.push_back( m_clusvec[ic]->getCaloHits()[ihit]->getPosition().x() );
//    m_Clus_hity.push_back( m_clusvec[ic]->getCaloHits()[ihit]->getPosition().y() );
//    m_Clus_hitz.push_back( m_clusvec[ic]->getCaloHits()[ihit]->getPosition().z() );
//    m_Clus_hitE.push_back( m_clusvec[ic]->getCaloHits()[ihit]->getEnergy() );
//    m_Clus_hittag.push_back( m_clusvec[ic]->getHitsE() );
//  }

  }

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
    std::vector<TrackState> AllTrackStates = m_trkCol[itrk]->getTrackStates("Input");
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
  t_Layers->Write();
  t_Hough->Write(); // yyy
  t_Cluster->Write();
  t_Tower->Write();
  t_Track->Write();
  t_axis->Write();
  m_wfile->Close();
  delete m_wfile, t_SimBar, t_Layers, t_Hough, t_Cluster, t_Tower, t_Track, t_axis;

  delete m_pMCParticleCreator;
  delete m_pTrackCreator; 
  delete m_pCaloHitsCreator;
  delete m_pOutputCreator;

  delete r_MCParticleCol;
  delete r_MCPRecoCaloAssoCol; 
  for(auto iter : r_TrackCols) delete iter; 
  //for(auto iter : r_ECalHitCols) delete iter; 
  //for(auto iter : r_HCalHitCols) delete iter; 
  for(auto iter : r_CaloHitCols) delete iter; 
  r_TrackCols.clear();
  //r_ECalHitCols.clear();
  //r_HCalHitCols.clear();
  r_CaloHitCols.clear();


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
  m_simBar_bar.clear();

  m_HcalHit_x.clear();
  m_HcalHit_y.clear();
  m_HcalHit_z.clear();
  m_HcalHit_E.clear();
  m_HcalHit_layer.clear();
}

void PandoraPlusPFAlg::ClearLayer(){
  m_NshowerU=-99;
  m_NshowerV=-99;
  m_barShowerU_tag.clear();
  m_barShowerU_x.clear();
  m_barShowerU_y.clear();
  m_barShowerU_z.clear();
  m_barShowerU_E.clear();
  m_barShowerU_T.clear();
  m_barShowerU_module.clear();
  m_barShowerU_part.clear();
  m_barShowerU_stave.clear();
  m_barShowerU_dlayer.clear();
  m_barShowerU_slayer.clear();
  m_barShowerU_bar.clear();
  m_barShowerV_tag.clear();
  m_barShowerV_x.clear();
  m_barShowerV_y.clear();
  m_barShowerV_z.clear();
  m_barShowerV_E.clear();
  m_barShowerV_T.clear();
  m_barShowerV_module.clear();
  m_barShowerV_part.clear();
  m_barShowerV_stave.clear();
  m_barShowerV_dlayer.clear();
  m_barShowerV_slayer.clear();
  m_barShowerV_bar.clear();
}

// yyy
void PandoraPlusPFAlg::ClearHough(){
  m_NHoughAxisU = -99;
  m_NHoughAxisV = -99;
  m_NConeAxisU = -99;
  m_NConeAxisV = -99;
  m_NMergedAxisU = -99;
  m_NMergedAxisV = -99;
  m_halfclusV_tag.clear();
  m_houghV_x.clear();
  m_houghV_y.clear();
  m_houghV_z.clear();
  m_houghV_E.clear();
  m_halfclusU_tag.clear();
  m_houghU_x.clear();
  m_houghU_y.clear();
  m_houghU_z.clear();
  m_houghU_E.clear();
}


void PandoraPlusPFAlg::ClearCluster(){
  m_Nclus = -99;
  m_Nmc = -99;
  m_totE = -99;
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


void PandoraPlusPFAlg::ClearAxis(){
  m_NaxisU = -99;
  m_NaxisV = -99;
  m_axisU_Nhit.clear(); 
  m_axisU_type.clear(); 
  m_axisU_Nmcp.clear();
  m_axisV_Nhit.clear(); 
  m_axisV_type.clear();
  m_axisV_Nmcp.clear();
  m_axisU_x.clear();
  m_axisU_y.clear();
  m_axisU_z.clear();
  m_axisU_E.clear();
  m_axisU_truthFracMax.clear();
  m_axisV_x.clear();
  m_axisV_y.clear();
  m_axisV_z.clear();
  m_axisV_E.clear();
  m_axisV_truthFracMax.clear();
  m_axisUhit_tag.clear();
  m_axisUhit_type.clear();
  m_axisUhit_x.clear(); 
  m_axisUhit_y.clear(); 
  m_axisUhit_z.clear(); 
  m_axisUhit_E.clear();
  m_axisVhit_tag.clear();
  m_axisVhit_type.clear();
  m_axisVhit_x.clear();
  m_axisVhit_y.clear();
  m_axisVhit_z.clear();
  m_axisVhit_E.clear();
  
}

void PandoraPlusPFAlg::ClearTower(){
  m_module = -1;
  m_part = -1;
  m_stave = -1;
  m_HFClusU_Nhit.clear();
  m_HFClusU_type.clear();
  m_HFClusV_Nhit.clear();
  m_HFClusV_type.clear();
  m_HFClusU_x.clear();
  m_HFClusU_y.clear();
  m_HFClusU_z.clear();
  m_HFClusU_E.clear();
  m_HFClusV_x.clear();
  m_HFClusV_y.clear();
  m_HFClusV_z.clear();
  m_HFClusV_E.clear();
  m_HFClusUhit_tag.clear();
  m_HFClusUhit_type.clear();
  m_HFClusUhit_x.clear();
  m_HFClusUhit_y.clear();
  m_HFClusUhit_z.clear();
  m_HFClusUhit_E.clear();
  m_HFClusVhit_tag.clear();
  m_HFClusVhit_type.clear();
  m_HFClusVhit_x.clear();
  m_HFClusVhit_y.clear();
  m_HFClusVhit_z.clear();
  m_HFClusVhit_E.clear();

}

#endif
