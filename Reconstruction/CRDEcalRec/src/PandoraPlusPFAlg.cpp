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
    t_MCParticle = new TTree("MCParticle", "MCParticle");
    t_SimBar = new TTree("SimBarHit", "SimBarHit");
    t_LocalMax = new TTree("LocalMax", "LocalMax");
    t_Layers = new TTree("RecLayers","RecLayers");
    t_axis = new TTree("Axis","Axis");
    t_HalfCluster = new TTree("HalfCluster","HalfCluster");
    t_Tower = new TTree("Tower", "Tower");
    t_Cluster = new TTree("RecClusters", "RecClusters");

    //MC particle 
    t_MCParticle->Branch("mcPdgid",     &m_mcPdgid);
    t_MCParticle->Branch("mcStatus",    &m_mcStatus);
    t_MCParticle->Branch("mcPx", &m_mcPx);
    t_MCParticle->Branch("mcPy", &m_mcPy);
    t_MCParticle->Branch("mcPz", &m_mcPz);
    t_MCParticle->Branch("mcEn", &m_mcEn);
    t_MCParticle->Branch("mcMass", &m_mcMass);
    t_MCParticle->Branch("mcCharge", &m_mcCharge);
    t_MCParticle->Branch("mcEPx", &m_mcEPx);
    t_MCParticle->Branch("mcEPy", &m_mcEPy);
    t_MCParticle->Branch("mcEPz", &m_mcEPz);


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
    t_SimBar->Branch("simBar_truthMC_tag", &m_simBar_truthMC_tag);
    t_SimBar->Branch("simBar_truthMC_pid", &m_simBar_truthMC_pid);
    t_SimBar->Branch("simBar_truthMC_px", &m_simBar_truthMC_px);
    t_SimBar->Branch("simBar_truthMC_py", &m_simBar_truthMC_py);
    t_SimBar->Branch("simBar_truthMC_pz", &m_simBar_truthMC_pz);
    t_SimBar->Branch("simBar_truthMC_E", &m_simBar_truthMC_E);
    t_SimBar->Branch("simBar_truthMC_EPx", &m_simBar_truthMC_EPx);
    t_SimBar->Branch("simBar_truthMC_EPy", &m_simBar_truthMC_EPy);
    t_SimBar->Branch("simBar_truthMC_EPz", &m_simBar_truthMC_EPz);
    t_SimBar->Branch("simBar_truthMC_weight", &m_simBar_truthMC_weight);
    t_SimBar->Branch("totE_HcalSim", &m_totE_HcalSim);
    t_SimBar->Branch("HcalHit_x", &m_HcalHit_x); 
    t_SimBar->Branch("HcalHit_y", &m_HcalHit_y); 
    t_SimBar->Branch("HcalHit_z", &m_HcalHit_z); 
    t_SimBar->Branch("HcalHit_E", &m_HcalHit_E); 
    t_SimBar->Branch("HcalHit_layer", &m_HcalHit_layer); 
    t_SimBar->Branch("HcalHit_truthMC_tag", &m_HcalHit_truthMC_tag);
    t_SimBar->Branch("HcalHit_truthMC_pid", &m_HcalHit_truthMC_pid);
    t_SimBar->Branch("HcalHit_truthMC_px", &m_HcalHit_truthMC_px);
    t_SimBar->Branch("HcalHit_truthMC_py", &m_HcalHit_truthMC_py);
    t_SimBar->Branch("HcalHit_truthMC_pz", &m_HcalHit_truthMC_pz);
    t_SimBar->Branch("HcalHit_truthMC_E", &m_HcalHit_truthMC_E);
    t_SimBar->Branch("HcalHit_truthMC_EPx", &m_HcalHit_truthMC_EPx);
    t_SimBar->Branch("HcalHit_truthMC_EPy", &m_HcalHit_truthMC_EPy);
    t_SimBar->Branch("HcalHit_truthMC_EPz", &m_HcalHit_truthMC_EPz);
    t_SimBar->Branch("HcalHit_truthMC_weight", &m_HcalHit_truthMC_weight);    

    //ECAL local max
    t_LocalMax->Branch("NlocalMaxU", &m_NlmU);
    t_LocalMax->Branch("NlocalMaxV", &m_NlmV);
    t_LocalMax->Branch("localMaxU_tag", &m_localMaxU_tag);
    t_LocalMax->Branch("localMaxU_x", &m_localMaxU_x);
    t_LocalMax->Branch("localMaxU_y", &m_localMaxU_y);
    t_LocalMax->Branch("localMaxU_z", &m_localMaxU_z);
    t_LocalMax->Branch("localMaxU_E", &m_localMaxU_E);
    t_LocalMax->Branch("localMaxU_mc_tag", &m_localMaxU_mc_tag);
    t_LocalMax->Branch("localMaxU_mc_pdg", &m_localMaxU_mc_pdg);
    t_LocalMax->Branch("localMaxU_mc_px", &m_localMaxU_mc_px);
    t_LocalMax->Branch("localMaxU_mc_py", &m_localMaxU_mc_py);
    t_LocalMax->Branch("localMaxU_mc_pz", &m_localMaxU_mc_pz);
    t_LocalMax->Branch("localMaxU_mc_weight", &m_localMaxU_mc_weight);
    t_LocalMax->Branch("localMaxV_tag", &m_localMaxV_tag);
    t_LocalMax->Branch("localMaxV_x", &m_localMaxV_x);
    t_LocalMax->Branch("localMaxV_y", &m_localMaxV_y);
    t_LocalMax->Branch("localMaxV_z", &m_localMaxV_z);
    t_LocalMax->Branch("localMaxV_E", &m_localMaxV_E);
    t_LocalMax->Branch("localMaxV_mc_tag", &m_localMaxV_mc_tag);
    t_LocalMax->Branch("localMaxV_mc_pdg", &m_localMaxV_mc_pdg);
    t_LocalMax->Branch("localMaxV_mc_px", &m_localMaxV_mc_px);
    t_LocalMax->Branch("localMaxV_mc_py", &m_localMaxV_mc_py);
    t_LocalMax->Branch("localMaxV_mc_pz", &m_localMaxV_mc_pz);
    t_LocalMax->Branch("localMaxV_mc_weight", &m_localMaxV_mc_weight);

    //1D Showers
    t_Layers->Branch("NshowerU", &m_NshowerU);
    t_Layers->Branch("NshowerV", &m_NshowerV);
    t_Layers->Branch("barShowerU_tag", &m_barShowerU_tag);
    t_Layers->Branch("barShowerU_x", &m_barShowerU_x);
    t_Layers->Branch("barShowerU_y", &m_barShowerU_y);
    t_Layers->Branch("barShowerU_z", &m_barShowerU_z);
    t_Layers->Branch("barShowerU_E", &m_barShowerU_E);
    t_Layers->Branch("barShowerU_mc_tag", &m_barShowerU_mc_tag);
    t_Layers->Branch("barShowerU_mc_pdg", &m_barShowerU_mc_pdg);
    t_Layers->Branch("barShowerU_mc_px", &m_barShowerU_mc_px);
    t_Layers->Branch("barShowerU_mc_py", &m_barShowerU_mc_py);
    t_Layers->Branch("barShowerU_mc_pz", &m_barShowerU_mc_pz);
    t_Layers->Branch("barShowerU_mc_weight", &m_barShowerU_mc_weight);
    t_Layers->Branch("barShowerV_tag", &m_barShowerV_tag);
    t_Layers->Branch("barShowerV_x", &m_barShowerV_x);
    t_Layers->Branch("barShowerV_y", &m_barShowerV_y);
    t_Layers->Branch("barShowerV_z", &m_barShowerV_z);
    t_Layers->Branch("barShowerV_E", &m_barShowerV_E);
    t_Layers->Branch("barShowerV_mc_tag", &m_barShowerV_mc_tag);
    t_Layers->Branch("barShowerV_mc_pdg", &m_barShowerV_mc_pdg);
    t_Layers->Branch("barShowerV_mc_px", &m_barShowerV_mc_px);
    t_Layers->Branch("barShowerV_mc_py", &m_barShowerV_mc_py);
    t_Layers->Branch("barShowerV_mc_pz", &m_barShowerV_mc_pz);
    t_Layers->Branch("barShowerV_mc_weight", &m_barShowerV_mc_weight);

    //Axis
    t_axis->Branch("axisU_tag", &m_axisU_tag);
    t_axis->Branch("axisU_type", &m_axisU_type);
    t_axis->Branch("axisU_x", &m_axisU_x);
    t_axis->Branch("axisU_y", &m_axisU_y);
    t_axis->Branch("axisU_z", &m_axisU_z);
    t_axis->Branch("axisU_E", &m_axisU_E);
    t_axis->Branch("axisU_truth_tag", &m_axisU_truth_tag);
    t_axis->Branch("axisU_truth_MC_px", &m_axisU_truth_MC_px);
    t_axis->Branch("axisU_truth_MC_py", &m_axisU_truth_MC_py);
    t_axis->Branch("axisU_truth_MC_pz", &m_axisU_truth_MC_pz);
    t_axis->Branch("axisU_truth_MC_E", &m_axisU_truth_MC_E);
    t_axis->Branch("axisU_truth_MC_weight", &m_axisU_truth_MC_weight);
    t_axis->Branch("axisU_hit_tag", &m_axisU_hit_tag);
    t_axis->Branch("axisU_hit_x", &m_axisU_hit_x);
    t_axis->Branch("axisU_hit_y", &m_axisU_hit_y);
    t_axis->Branch("axisU_hit_z", &m_axisU_hit_z);
    t_axis->Branch("axisU_hit_E", &m_axisU_hit_E);
    t_axis->Branch("axisV_tag", &m_axisV_tag);
    t_axis->Branch("axisV_type", &m_axisV_type);
    t_axis->Branch("axisV_x", &m_axisV_x);
    t_axis->Branch("axisV_y", &m_axisV_y);
    t_axis->Branch("axisV_z", &m_axisV_z);
    t_axis->Branch("axisV_E", &m_axisV_E);
    t_axis->Branch("axisV_truth_tag", &m_axisV_truth_tag);
    t_axis->Branch("axisV_truth_MC_px", &m_axisV_truth_MC_px);
    t_axis->Branch("axisV_truth_MC_py", &m_axisV_truth_MC_py);
    t_axis->Branch("axisV_truth_MC_pz", &m_axisV_truth_MC_pz);
    t_axis->Branch("axisV_truth_MC_E", &m_axisV_truth_MC_E);
    t_axis->Branch("axisV_truth_MC_weight", &m_axisV_truth_MC_weight);
    t_axis->Branch("axisV_hit_tag", &m_axisV_hit_tag);
    t_axis->Branch("axisV_hit_x", &m_axisV_hit_x);
    t_axis->Branch("axisV_hit_y", &m_axisV_hit_y);
    t_axis->Branch("axisV_hit_z", &m_axisV_hit_z);
    t_axis->Branch("axisV_hit_E", &m_axisV_hit_E);    

    t_axis->Branch("emptyAxisU_tag", &m_emptyAxisU_tag);
    t_axis->Branch("emptyAxisU_x", &m_emptyAxisU_x);
    t_axis->Branch("emptyAxisU_y", &m_emptyAxisU_y);
    t_axis->Branch("emptyAxisU_z", &m_emptyAxisU_z);
    t_axis->Branch("emptyAxisU_E", &m_emptyAxisU_E);
    t_axis->Branch("emptyAxisV_tag", &m_emptyAxisV_tag);
    t_axis->Branch("emptyAxisV_x", &m_emptyAxisV_x);
    t_axis->Branch("emptyAxisV_y", &m_emptyAxisV_y);
    t_axis->Branch("emptyAxisV_z", &m_emptyAxisV_z);
    t_axis->Branch("emptyAxisV_E", &m_emptyAxisV_E);

    //Half clusters
    t_HalfCluster->Branch("totE_HFClusV", &m_totE_HFClusV);
    t_HalfCluster->Branch("HalfClusterV_x", &m_HalfClusterV_x);
    t_HalfCluster->Branch("HalfClusterV_y", &m_HalfClusterV_y);
    t_HalfCluster->Branch("HalfClusterV_z", &m_HalfClusterV_z);
    t_HalfCluster->Branch("HalfClusterV_E", &m_HalfClusterV_E);
    t_HalfCluster->Branch("HalfClusterV_tag", &m_HalfClusterV_tag);
    t_HalfCluster->Branch("HalfClusterV_type", &m_HalfClusterV_type);
    t_HalfCluster->Branch("HalfClusterV_hit_x", &m_HalfClusterV_hit_x);
    t_HalfCluster->Branch("HalfClusterV_hit_y", &m_HalfClusterV_hit_y);
    t_HalfCluster->Branch("HalfClusterV_hit_z", &m_HalfClusterV_hit_z);
    t_HalfCluster->Branch("HalfClusterV_hit_E", &m_HalfClusterV_hit_E);
    t_HalfCluster->Branch("HalfClusterV_hit_tag", &m_HalfClusterV_hit_tag);
    t_HalfCluster->Branch("HalfClusterV_truth_tag", &m_HalfClusterV_truth_tag);
    t_HalfCluster->Branch("HalfClusterV_truthMC_px", &m_HalfClusterV_truthMC_px);
    t_HalfCluster->Branch("HalfClusterV_truthMC_py", &m_HalfClusterV_truthMC_py);
    t_HalfCluster->Branch("HalfClusterV_truthMC_pz", &m_HalfClusterV_truthMC_pz);
    t_HalfCluster->Branch("HalfClusterV_truthMC_E", &m_HalfClusterV_truthMC_E);
    t_HalfCluster->Branch("HalfClusterV_truthMC_weight", &m_HalfClusterV_truthMC_weight);

    t_HalfCluster->Branch("totE_HFClusU", &m_totE_HFClusU);
    t_HalfCluster->Branch("HalfClusterU_x", &m_HalfClusterU_x);
    t_HalfCluster->Branch("HalfClusterU_y", &m_HalfClusterU_y);
    t_HalfCluster->Branch("HalfClusterU_z", &m_HalfClusterU_z);
    t_HalfCluster->Branch("HalfClusterU_E", &m_HalfClusterU_E);  
    t_HalfCluster->Branch("HalfClusterU_tag", &m_HalfClusterU_tag);  
    t_HalfCluster->Branch("HalfClusterU_type", &m_HalfClusterU_type);  
    t_HalfCluster->Branch("HalfClusterU_hit_x", &m_HalfClusterU_hit_x);
    t_HalfCluster->Branch("HalfClusterU_hit_y", &m_HalfClusterU_hit_y);
    t_HalfCluster->Branch("HalfClusterU_hit_z", &m_HalfClusterU_hit_z);
    t_HalfCluster->Branch("HalfClusterU_hit_E", &m_HalfClusterU_hit_E);
    t_HalfCluster->Branch("HalfClusterU_hit_tag", &m_HalfClusterU_hit_tag);
    t_HalfCluster->Branch("HalfClusterU_truth_tag", &m_HalfClusterU_truth_tag);
    t_HalfCluster->Branch("HalfClusterU_truthMC_px", &m_HalfClusterU_truthMC_px);
    t_HalfCluster->Branch("HalfClusterU_truthMC_py", &m_HalfClusterU_truthMC_py);
    t_HalfCluster->Branch("HalfClusterU_truthMC_pz", &m_HalfClusterU_truthMC_pz);
    t_HalfCluster->Branch("HalfClusterU_truthMC_E", &m_HalfClusterU_truthMC_E);
    t_HalfCluster->Branch("HalfClusterU_truthMC_weight", &m_HalfClusterU_truthMC_weight);


    //Tower
    t_Tower->Branch("towerID", towerID, "towerID[3]/I");
    t_Tower->Branch("NclusU", &m_NclusU);
    t_Tower->Branch("NclusV", &m_NclusV);
    t_Tower->Branch("totEn", &m_totEn);
    t_Tower->Branch("HalfClusterU_x", &m_HalfClusterU_x);
    t_Tower->Branch("HalfClusterU_y", &m_HalfClusterU_y);
    t_Tower->Branch("HalfClusterU_z", &m_HalfClusterU_z);
    t_Tower->Branch("HalfClusterU_E", &m_HalfClusterU_E);
    t_Tower->Branch("HalfClusterU_tag", &m_HalfClusterU_tag);
    t_Tower->Branch("HalfClusterU_type", &m_HalfClusterU_type);
    t_Tower->Branch("HalfClusterU_nTrk", &m_HalfClusterU_nTrk);
    t_Tower->Branch("HalfClusterV_x", &m_HalfClusterV_x);
    t_Tower->Branch("HalfClusterV_y", &m_HalfClusterV_y);
    t_Tower->Branch("HalfClusterV_z", &m_HalfClusterV_z);
    t_Tower->Branch("HalfClusterV_E", &m_HalfClusterV_E);
    t_Tower->Branch("HalfClusterV_tag", &m_HalfClusterV_tag);
    t_Tower->Branch("HalfClusterV_type", &m_HalfClusterV_type);
    t_Tower->Branch("HalfClusterV_nTrk", &m_HalfClusterV_nTrk);

    //Clusters
    t_Cluster->Branch("totE_Ecal", &m_totE_Ecal);
    t_Cluster->Branch("totE_Hcal", &m_totE_Hcal);
    t_Cluster->Branch("Nclus_Ecal", &m_Nclus_Ecal);
    t_Cluster->Branch("Nclus_Hcal", &m_Nclus_Hcal);
    t_Cluster->Branch("EcalClus_x", &m_EcalClus_x);
    t_Cluster->Branch("EcalClus_y", &m_EcalClus_y);
    t_Cluster->Branch("EcalClus_z", &m_EcalClus_z);
    t_Cluster->Branch("EcalClus_E", &m_EcalClus_E);
    t_Cluster->Branch("EcalClus_nTrk", &m_EcalClus_nTrk);
    t_Cluster->Branch("EcalClus_ptrk", &m_EcalClus_pTrk);
    t_Cluster->Branch("EcalClus_typeU", &m_EcalClus_typeU);
    t_Cluster->Branch("EcalClus_typeV", &m_EcalClus_typeV);
    t_Cluster->Branch("EcalClus_hitU_x", &m_EcalClus_hitU_x);
    t_Cluster->Branch("EcalClus_hitU_y", &m_EcalClus_hitU_y);
    t_Cluster->Branch("EcalClus_hitU_z", &m_EcalClus_hitU_z);
    t_Cluster->Branch("EcalClus_hitU_E", &m_EcalClus_hitU_E);
    t_Cluster->Branch("EcalClus_hitU_tag", &m_EcalClus_hitU_tag);
    t_Cluster->Branch("EcalClus_hitV_x", &m_EcalClus_hitV_x);
    t_Cluster->Branch("EcalClus_hitV_y", &m_EcalClus_hitV_y);
    t_Cluster->Branch("EcalClus_hitV_z", &m_EcalClus_hitV_z);
    t_Cluster->Branch("EcalClus_hitV_E", &m_EcalClus_hitV_E);
    t_Cluster->Branch("EcalClus_hitV_tag", &m_EcalClus_hitV_tag);
    t_Cluster->Branch("EcalClus_truthMC_tag", &m_EcalClus_truthMC_tag);
    t_Cluster->Branch("EcalClus_truthMC_pid", &m_EcalClus_truthMC_pid);
    t_Cluster->Branch("EcalClus_truthMC_px", &m_EcalClus_truthMC_px);
    t_Cluster->Branch("EcalClus_truthMC_py", &m_EcalClus_truthMC_py);
    t_Cluster->Branch("EcalClus_truthMC_pz", &m_EcalClus_truthMC_pz);
    t_Cluster->Branch("EcalClus_truthMC_E", &m_EcalClus_truthMC_E);
    t_Cluster->Branch("EcalClus_truthMC_EPx", &m_EcalClus_truthMC_EPx);
    t_Cluster->Branch("EcalClus_truthMC_EPy", &m_EcalClus_truthMC_EPy);
    t_Cluster->Branch("EcalClus_truthMC_EPz", &m_EcalClus_truthMC_EPz);
    t_Cluster->Branch("EcalClus_truthMC_weight", &m_EcalClus_truthMC_weight);
    t_Cluster->Branch("HcalClus_x", &m_HcalClus_x);
    t_Cluster->Branch("HcalClus_y", &m_HcalClus_y);
    t_Cluster->Branch("HcalClus_z", &m_HcalClus_z);
    t_Cluster->Branch("HcalClus_E", &m_HcalClus_E);
    t_Cluster->Branch("HcalClus_nTrk", &m_HcalClus_nTrk);
    t_Cluster->Branch("HcalClus_ptrk", &m_HcalClus_pTrk);
    t_Cluster->Branch("HcalClus_truthMC_tag", &m_HcalClus_truthMC_tag);
    t_Cluster->Branch("HcalClus_truthMC_pid", &m_HcalClus_truthMC_pid);
    t_Cluster->Branch("HcalClus_truthMC_px", &m_HcalClus_truthMC_px);
    t_Cluster->Branch("HcalClus_truthMC_py", &m_HcalClus_truthMC_py);
    t_Cluster->Branch("HcalClus_truthMC_pz", &m_HcalClus_truthMC_pz);
    t_Cluster->Branch("HcalClus_truthMC_E", &m_HcalClus_truthMC_E);
    t_Cluster->Branch("HcalClus_truthMC_EPx", &m_HcalClus_truthMC_EPx);
    t_Cluster->Branch("HcalClus_truthMC_EPy", &m_HcalClus_truthMC_EPy);
    t_Cluster->Branch("HcalClus_truthMC_EPz", &m_HcalClus_truthMC_EPz);
    t_Cluster->Branch("HcalClus_truthMC_weight", &m_HcalClus_truthMC_weight);
    t_Cluster->Branch("SimpleHcalClus_x", &m_SimpleHcalClus_x);
    t_Cluster->Branch("SimpleHcalClus_y", &m_SimpleHcalClus_y);
    t_Cluster->Branch("SimpleHcalClus_z", &m_SimpleHcalClus_z);
    t_Cluster->Branch("SimpleHcalClus_E", &m_SimpleHcalClus_E);
    t_Cluster->Branch("SimpleHcalClus_nTrk", &m_SimpleHcalClus_nTrk);
    t_Cluster->Branch("SimpleHcalClus_ptrk", &m_SimpleHcalClus_pTrk);
    t_Cluster->Branch("SimpleHcalClus_truthMC_tag", &m_SimpleHcalClus_truthMC_tag);
    t_Cluster->Branch("SimpleHcalClus_truthMC_pid", &m_SimpleHcalClus_truthMC_pid);
    t_Cluster->Branch("SimpleHcalClus_truthMC_px", &m_SimpleHcalClus_truthMC_px);
    t_Cluster->Branch("SimpleHcalClus_truthMC_py", &m_SimpleHcalClus_truthMC_py);
    t_Cluster->Branch("SimpleHcalClus_truthMC_pz", &m_SimpleHcalClus_truthMC_pz);
    t_Cluster->Branch("SimpleHcalClus_truthMC_E", &m_SimpleHcalClus_truthMC_E);
    t_Cluster->Branch("SimpleHcalClus_truthMC_EPx", &m_SimpleHcalClus_truthMC_EPx);
    t_Cluster->Branch("SimpleHcalClus_truthMC_EPy", &m_SimpleHcalClus_truthMC_EPy);
    t_Cluster->Branch("SimpleHcalClus_truthMC_EPz", &m_SimpleHcalClus_truthMC_EPz);
    t_Cluster->Branch("SimpleHcalClus_truthMC_weight", &m_SimpleHcalClus_truthMC_weight);
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
  m_pMCParticleCreator->CreateMCParticle( m_DataCol, *r_MCParticleCol );

  if(m_useTruthTrk) m_pTrackCreator->CreateTracksFromMCParticle(m_DataCol, *r_MCParticleCol);
  else m_pTrackCreator->CreateTracks( m_DataCol, r_TrackCols, r_MCPTrkAssoCol );

  m_pCaloHitsCreator->CreateCaloHits( m_DataCol, r_CaloHitCols, map_readout_decoder, map_CaloMCPAssoCols );

  //Perform PFA algorithm
  m_algorithmManager.RunAlgorithm( m_DataCol );

  m_pOutputCreator->CreateRecCaloHits( m_DataCol, w_RecCaloCol );
  m_pOutputCreator->CreateCluster( m_DataCol, w_ClusterCollection );
  m_pOutputCreator->CreatePFO( m_DataCol, w_ReconstructedParticleCollection );


// yyy_endrec = clock();  // 重建结束的时间

cout<<"Write tuples"<<endl;
  //---------------------Write Ana tuples-------------------------
  // MC particles  
  ClearMCParticle();
  std::vector<edm4hep::MCParticle> m_MCPCol = m_DataCol.collectionMap_MC[name_MCParticleCol.value()];
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
  t_MCParticle->Fill();


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

    auto truthMap = p_hitbar->getLinkedMCP();
    for(auto iter: truthMap){
      m_simBar_truthMC_tag.push_back(ibar);
      m_simBar_truthMC_pid.push_back(iter.first.getPDG());
      m_simBar_truthMC_px.push_back(iter.first.getMomentum().x);
      m_simBar_truthMC_py.push_back(iter.first.getMomentum().y);
      m_simBar_truthMC_pz.push_back(iter.first.getMomentum().z);
      m_simBar_truthMC_E.push_back(iter.first.getEnergy());
      m_simBar_truthMC_EPx.push_back(iter.first.getEndpoint().x);
      m_simBar_truthMC_EPy.push_back(iter.first.getEndpoint().y);
      m_simBar_truthMC_EPz.push_back(iter.first.getEndpoint().z);
      m_simBar_truthMC_weight.push_back(iter.second);
    }
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

    auto truthMap = m_hcalHitsCol[ihit]->getLinkedMCP();
    for(auto iter: truthMap){
      m_HcalHit_truthMC_tag.push_back(ihit);
      m_HcalHit_truthMC_pid.push_back(iter.first.getPDG());
      m_HcalHit_truthMC_px.push_back(iter.first.getMomentum().x);
      m_HcalHit_truthMC_py.push_back(iter.first.getMomentum().y);
      m_HcalHit_truthMC_pz.push_back(iter.first.getMomentum().z);
      m_HcalHit_truthMC_E.push_back(iter.first.getEnergy());
      m_HcalHit_truthMC_EPx.push_back(iter.first.getEndpoint().x);
      m_HcalHit_truthMC_EPy.push_back(iter.first.getEndpoint().y);
      m_HcalHit_truthMC_EPz.push_back(iter.first.getEndpoint().z);
      m_HcalHit_truthMC_weight.push_back(iter.second);
    }    

  }
  t_SimBar->Fill();

  //Save localMax
  ClearLocalMax();
  std::vector<PandoraPlus::CaloHalfCluster*> m_halfclusters; m_halfclusters.clear();
  for(int i=0; i<m_DataCol.map_HalfCluster["HalfClusterColU"].size(); i++)
    m_halfclusters.push_back( m_DataCol.map_HalfCluster["HalfClusterColU"][i].get() );

  std::vector<const Calo1DCluster*> m_local_max; m_local_max.clear();
  for(int ic=0;ic<m_halfclusters.size(); ic++){
    std::vector<const Calo1DCluster*> tmp_shower = m_halfclusters[ic]->getLocalMaxCol("AllLocalMax");
    m_local_max.insert(m_local_max.end(), tmp_shower.begin(), tmp_shower.end());
  }
  for(int il=0; il<m_local_max.size(); il++){
    m_localMaxU_tag.push_back( il );
    m_localMaxU_x.push_back( m_local_max[il]->getPos().x() );
    m_localMaxU_y.push_back( m_local_max[il]->getPos().y() );
    m_localMaxU_z.push_back( m_local_max[il]->getPos().z() );
    m_localMaxU_E.push_back( m_local_max[il]->getEnergy() );

    auto truthMap = m_local_max[il]->getLinkedMCP();
    for(auto iter: truthMap){
      m_localMaxU_mc_tag.push_back(il);
      m_localMaxU_mc_pdg.push_back(iter.first.getPDG());
      m_localMaxU_mc_px.push_back(iter.first.getMomentum().x);
      m_localMaxU_mc_py.push_back(iter.first.getMomentum().y);
      m_localMaxU_mc_pz.push_back(iter.first.getMomentum().z);
      m_localMaxU_mc_weight.push_back(iter.second);
    }
  }
  m_halfclusters.clear();
  m_local_max.clear(); 
  for(int i=0; i<m_DataCol.map_HalfCluster["HalfClusterColV"].size(); i++)
    m_halfclusters.push_back( m_DataCol.map_HalfCluster["HalfClusterColV"][i].get() );

  for(int ic=0;ic<m_halfclusters.size(); ic++){
    std::vector<const Calo1DCluster*> tmp_shower = m_halfclusters[ic]->getLocalMaxCol("AllLocalMax");
    m_local_max.insert(m_local_max.end(), tmp_shower.begin(), tmp_shower.end());
  }
  for(int il=0; il<m_local_max.size(); il++){
    m_localMaxV_tag.push_back( il );
    m_localMaxV_x.push_back( m_local_max[il]->getPos().x() );
    m_localMaxV_y.push_back( m_local_max[il]->getPos().y() );
    m_localMaxV_z.push_back( m_local_max[il]->getPos().z() );
    m_localMaxV_E.push_back( m_local_max[il]->getEnergy() );

    auto truthMap = m_local_max[il]->getLinkedMCP();
    for(auto iter: truthMap){
      m_localMaxV_mc_tag.push_back(il);
      m_localMaxV_mc_pdg.push_back(iter.first.getPDG());
      m_localMaxV_mc_px.push_back(iter.first.getMomentum().x);
      m_localMaxV_mc_py.push_back(iter.first.getMomentum().y);
      m_localMaxV_mc_pz.push_back(iter.first.getMomentum().z);
      m_localMaxV_mc_weight.push_back(iter.second);
    }
  }
  t_LocalMax->Fill();

  //Save 1DCluster
  ClearLayer();
  m_halfclusters.clear();
  for(int i=0; i<m_DataCol.map_HalfCluster["ESHalfClusterU"].size(); i++)
    m_halfclusters.push_back( m_DataCol.map_HalfCluster["ESHalfClusterU"][i].get() );

  m_local_max.clear();
  for(int ic=0;ic<m_halfclusters.size(); ic++){
    //std::vector<const Calo1DCluster*> tmp_shower = m_halfclusters[ic]->getLocalMaxCol("AllLocalMax");
    std::vector<const CaloHalfCluster*> m_axis = m_halfclusters[ic]->getHalfClusterCol("MergedAxis");
    if(m_axis.size()>0){
      std::vector<const Calo1DCluster*> tmp_shower = m_axis[0]->getCluster();
      m_local_max.insert(m_local_max.end(), tmp_shower.begin(), tmp_shower.end());
    }
  }
  for(int il=0; il<m_local_max.size(); il++){
    m_barShowerU_tag.push_back( il );
    m_barShowerU_x.push_back( m_local_max[il]->getPos().x() );
    m_barShowerU_y.push_back( m_local_max[il]->getPos().y() );
    m_barShowerU_z.push_back( m_local_max[il]->getPos().z() );
    m_barShowerU_E.push_back( m_local_max[il]->getEnergy() );  

    auto truthMap = m_local_max[il]->getLinkedMCP();
    for(auto iter: truthMap){
      m_barShowerU_mc_tag.push_back(il);
      m_barShowerU_mc_pdg.push_back(iter.first.getPDG());
      m_barShowerU_mc_px.push_back(iter.first.getMomentum().x);
      m_barShowerU_mc_py.push_back(iter.first.getMomentum().y);
      m_barShowerU_mc_pz.push_back(iter.first.getMomentum().z);
      m_barShowerU_mc_weight.push_back(iter.second);
    }
  }
  m_halfclusters.clear();
  m_local_max.clear();
  for(int i=0; i<m_DataCol.map_HalfCluster["ESHalfClusterV"].size(); i++)
    m_halfclusters.push_back( m_DataCol.map_HalfCluster["ESHalfClusterV"][i].get() );

  for(int ic=0;ic<m_halfclusters.size(); ic++){
    //std::vector<const Calo1DCluster*> tmp_shower = m_halfclusters[ic]->getLocalMaxCol("AllLocalMax");
    std::vector<const CaloHalfCluster*> m_axis = m_halfclusters[ic]->getHalfClusterCol("MergedAxis");
    if(m_axis.size()>0){
      std::vector<const Calo1DCluster*> tmp_shower = m_axis[0]->getCluster();
      m_local_max.insert(m_local_max.end(), tmp_shower.begin(), tmp_shower.end());
    }
  }
  for(int il=0; il<m_local_max.size(); il++){
    m_barShowerV_tag.push_back( il );
    m_barShowerV_x.push_back( m_local_max[il]->getPos().x() );
    m_barShowerV_y.push_back( m_local_max[il]->getPos().y() );
    m_barShowerV_z.push_back( m_local_max[il]->getPos().z() );
    m_barShowerV_E.push_back( m_local_max[il]->getEnergy() );

    auto truthMap = m_local_max[il]->getLinkedMCP();
    for(auto iter: truthMap){
      m_barShowerV_mc_tag.push_back(il);
      m_barShowerV_mc_pdg.push_back(iter.first.getPDG());
      m_barShowerV_mc_px.push_back(iter.first.getMomentum().x);
      m_barShowerV_mc_py.push_back(iter.first.getMomentum().y);
      m_barShowerV_mc_pz.push_back(iter.first.getMomentum().z);
      m_barShowerV_mc_weight.push_back(iter.second);
    }
  }
  t_Layers->Fill();


  //Axis
  ClearAxis();
  std::vector<const PandoraPlus::CaloHalfCluster*> m_halfclusterV; m_halfclusterV.clear();
  std::vector<const PandoraPlus::CaloHalfCluster*> m_halfclusterU; m_halfclusterU.clear();
  for(int i=0; i<m_DataCol.map_HalfCluster["HalfClusterColU"].size(); i++){
    m_halfclusterU.push_back( m_DataCol.map_HalfCluster["HalfClusterColU"][i].get() );
  }
  for(int i=0; i<m_DataCol.map_HalfCluster["HalfClusterColV"].size(); i++){
    m_halfclusterV.push_back( m_DataCol.map_HalfCluster["HalfClusterColV"][i].get() );
  }
  int axisU_index=0;
  int axisV_index=0;
  for(int i=0; i<m_halfclusterU.size(); i++){  // loop half cluster U
    std::vector<const PandoraPlus::CaloHalfCluster*> m_mergedaxisU = m_halfclusterU[i]->getHalfClusterCol("MergedAxis");
    for(int ita=0; ita<m_mergedaxisU.size(); ita++){ // loop  axis U
      // General information of the axis
      m_axisU_tag.push_back(axisU_index);
      m_axisU_type.push_back(m_mergedaxisU[ita]->getType());
      m_axisU_x.push_back(m_mergedaxisU[ita]->getPos().x());
      m_axisU_y.push_back(m_mergedaxisU[ita]->getPos().y());
      m_axisU_z.push_back(m_mergedaxisU[ita]->getPos().z());
      m_axisU_E.push_back(m_mergedaxisU[ita]->getEnergy());

      // MC truth information of the Axis
      auto truthMap = m_mergedaxisU[ita]->getLinkedMCP();
      for(auto iter: truthMap){
        m_axisU_truth_tag.push_back(axisU_index);
        m_axisU_truth_MC_px.push_back(iter.first.getMomentum().x);
        m_axisU_truth_MC_py.push_back(iter.first.getMomentum().y);
        m_axisU_truth_MC_pz.push_back(iter.first.getMomentum().z);
        m_axisU_truth_MC_E.push_back(iter.first.getEnergy());
        m_axisU_truth_MC_weight.push_back(iter.second);
      }
      // Hits on axis
      for(int ilm=0; ilm<m_mergedaxisU[ita]->getCluster().size(); ilm++){ // loop local max
        m_axisU_hit_tag.push_back(axisU_index);
        m_axisU_hit_x.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getPos().x() );
        m_axisU_hit_y.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getPos().y() );
        m_axisU_hit_z.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getPos().z() );
        m_axisU_hit_E.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getEnergy()  );
      }

      axisU_index++;
    }
  }
  for(int i=0; i<m_halfclusterV.size(); i++){  // loop half cluster V
    std::vector<const PandoraPlus::CaloHalfCluster*> m_mergedaxisV = m_halfclusterV[i]->getHalfClusterCol("MergedAxis");
    for(int ita=0; ita<m_mergedaxisV.size(); ita++){ // loop  axis V
      // General information of the axis
      m_axisV_tag.push_back(axisV_index);
      m_axisV_type.push_back(m_mergedaxisV[ita]->getType());
      m_axisV_x.push_back(m_mergedaxisV[ita]->getPos().x());
      m_axisV_y.push_back(m_mergedaxisV[ita]->getPos().y());
      m_axisV_z.push_back(m_mergedaxisV[ita]->getPos().z());
      m_axisV_E.push_back(m_mergedaxisV[ita]->getEnergy());

      // MC truth information of the Axis
      auto truthMap = m_mergedaxisV[ita]->getLinkedMCP();
      for(auto iter: truthMap){
        m_axisV_truth_tag.push_back(axisV_index);
        m_axisV_truth_MC_px.push_back(iter.first.getMomentum().x);
        m_axisV_truth_MC_py.push_back(iter.first.getMomentum().y);
        m_axisV_truth_MC_pz.push_back(iter.first.getMomentum().z);
        m_axisV_truth_MC_E.push_back(iter.first.getEnergy());
        m_axisV_truth_MC_weight.push_back(iter.second);
      }
      // Hits on axis
      for(int ilm=0; ilm<m_mergedaxisV[ita]->getCluster().size(); ilm++){
        m_axisV_hit_tag.push_back(axisV_index);
        m_axisV_hit_x.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getPos().x());
        m_axisV_hit_y.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getPos().y());
        m_axisV_hit_z.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getPos().z());
        m_axisV_hit_E.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getEnergy());
      }

      axisV_index++;
    }
  }
  m_halfclusterU.clear();
  m_halfclusterV.clear();
  for(int i=0; i<m_DataCol.map_HalfCluster["emptyHalfClusterU"].size(); i++){
    m_halfclusterU.push_back( m_DataCol.map_HalfCluster["emptyHalfClusterU"][i].get() );
  }
  for(int i=0; i<m_DataCol.map_HalfCluster["emptyHalfClusterV"].size(); i++){
    m_halfclusterV.push_back( m_DataCol.map_HalfCluster["emptyHalfClusterV"][i].get() );
  }
  for(int i=0; i<m_halfclusterU.size(); i++){
    m_emptyAxisU_tag.push_back(m_halfclusterU[i]->getType());
    m_emptyAxisU_x.push_back(m_halfclusterU[i]->getPos().x());
    m_emptyAxisU_y.push_back(m_halfclusterU[i]->getPos().y());
    m_emptyAxisU_z.push_back(m_halfclusterU[i]->getPos().z());
    m_emptyAxisU_E.push_back(m_halfclusterU[i]->getEnergy());
  }
  for(int i=0; i<m_halfclusterV.size(); i++){
    m_emptyAxisV_tag.push_back(m_halfclusterV[i]->getType());
    m_emptyAxisV_x.push_back(m_halfclusterV[i]->getPos().x());
    m_emptyAxisV_y.push_back(m_halfclusterV[i]->getPos().y());
    m_emptyAxisV_z.push_back(m_halfclusterV[i]->getPos().z());
    m_emptyAxisV_E.push_back(m_halfclusterV[i]->getEnergy());
  }

  t_axis->Fill();


  //Half cluster
  ClearHalfCluster();
  m_halfclusterV.clear();
  m_halfclusterU.clear();
  m_totE_HFClusV = 0;
  m_totE_HFClusU = 0;
  for(int i=0; i<m_DataCol.map_HalfCluster["ESHalfClusterU"].size(); i++){
    m_halfclusterU.push_back( m_DataCol.map_HalfCluster["ESHalfClusterU"][i]->getHalfClusterCol("MergedAxis")[0] );
  }
  for(int i=0; i<m_DataCol.map_HalfCluster["ESHalfClusterV"].size(); i++){
    m_halfclusterV.push_back( m_DataCol.map_HalfCluster["ESHalfClusterV"][i]->getHalfClusterCol("MergedAxis")[0] );
  }  
  for(int i=0; i<m_halfclusterV.size(); i++){
    m_HalfClusterV_x.push_back(m_halfclusterV[i]->getPos().x());
    m_HalfClusterV_y.push_back(m_halfclusterV[i]->getPos().y());
    m_HalfClusterV_z.push_back(m_halfclusterV[i]->getPos().z());
    m_HalfClusterV_E.push_back(m_halfclusterV[i]->getEnergy());
    m_HalfClusterV_tag.push_back(i);
    m_HalfClusterV_type.push_back(m_halfclusterV[i]->getType());
    m_totE_HFClusV += m_halfclusterV[i]->getEnergy();    

      // MC truth information of the HFCluster
      auto truthMap = m_halfclusterV[i]->getLinkedMCP();
      for(auto iter: truthMap){
        m_HalfClusterV_truth_tag.push_back(i);
        m_HalfClusterV_truthMC_px.push_back(iter.first.getMomentum().x);
        m_HalfClusterV_truthMC_py.push_back(iter.first.getMomentum().y);
        m_HalfClusterV_truthMC_pz.push_back(iter.first.getMomentum().z);
        m_HalfClusterV_truthMC_E.push_back(iter.first.getEnergy());
        m_HalfClusterV_truthMC_weight.push_back(iter.second);
      }

      // 1DClusters (hits)
      for(int ilm=0; ilm<m_halfclusterV[i]->getCluster().size(); ilm++){ 
        m_HalfClusterV_hit_tag.push_back(i);
        m_HalfClusterV_hit_x.push_back( m_halfclusterV[i]->getCluster()[ilm]->getPos().x() );
        m_HalfClusterV_hit_y.push_back( m_halfclusterV[i]->getCluster()[ilm]->getPos().y() );
        m_HalfClusterV_hit_z.push_back( m_halfclusterV[i]->getCluster()[ilm]->getPos().z() );
        m_HalfClusterV_hit_E.push_back( m_halfclusterV[i]->getCluster()[ilm]->getEnergy()  );
      }
  }
  for(int i=0; i<m_halfclusterU.size(); i++){
    m_HalfClusterU_x.push_back(m_halfclusterU[i]->getPos().x());
    m_HalfClusterU_y.push_back(m_halfclusterU[i]->getPos().y());
    m_HalfClusterU_z.push_back(m_halfclusterU[i]->getPos().z());
    m_HalfClusterU_E.push_back(m_halfclusterU[i]->getEnergy());
    m_HalfClusterU_tag.push_back(i);
    m_HalfClusterU_type.push_back(m_halfclusterU[i]->getType());
    m_totE_HFClusU += m_halfclusterU[i]->getEnergy();

      // MC truth information of the HFCluster
      auto truthMap = m_halfclusterU[i]->getLinkedMCP();
      for(auto iter: truthMap){
        m_HalfClusterU_truth_tag.push_back(i);
        m_HalfClusterU_truthMC_px.push_back(iter.first.getMomentum().x);
        m_HalfClusterU_truthMC_py.push_back(iter.first.getMomentum().y);
        m_HalfClusterU_truthMC_pz.push_back(iter.first.getMomentum().z);
        m_HalfClusterU_truthMC_E.push_back(iter.first.getEnergy());
        m_HalfClusterU_truthMC_weight.push_back(iter.second);
      }

      // 1DClusters (hits)
      for(int ilm=0; ilm<m_halfclusterU[i]->getCluster().size(); ilm++){
        m_HalfClusterU_hit_tag.push_back(i);
        m_HalfClusterU_hit_x.push_back( m_halfclusterU[i]->getCluster()[ilm]->getPos().x() );
        m_HalfClusterU_hit_y.push_back( m_halfclusterU[i]->getCluster()[ilm]->getPos().y() );
        m_HalfClusterU_hit_z.push_back( m_halfclusterU[i]->getCluster()[ilm]->getPos().z() );
        m_HalfClusterU_hit_E.push_back( m_halfclusterU[i]->getCluster()[ilm]->getEnergy()  );
      }
  }
  t_HalfCluster->Fill();


  //Tower
  ClearTower();
  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> m_tower = m_DataCol.map_CaloCluster["ESTower"];
  for(int it=0; it<m_tower.size(); it++){
    ClearTower();
    towerID[0] = m_tower[it]->getTowerID()[0][0];
    towerID[1] = m_tower[it]->getTowerID()[0][1];
    towerID[2] = m_tower[it]->getTowerID()[0][2];

    std::vector<const CaloHalfCluster*> m_HFClusU = m_tower[it]->getHalfClusterUCol("ESHalfClusterU");
    std::vector<const CaloHalfCluster*> m_HFClusV = m_tower[it]->getHalfClusterVCol("ESHalfClusterV");

    m_NclusU = m_HFClusU.size();
    m_NclusV = m_HFClusV.size();
    m_totEn = m_tower[it]->getEnergy();

    for(int ic=0; ic<m_NclusU; ic++){
      m_HalfClusterU_tag.push_back(it);
      m_HalfClusterU_x.push_back(m_HFClusU[ic]->getPos().x());
      m_HalfClusterU_y.push_back(m_HFClusU[ic]->getPos().y());
      m_HalfClusterU_z.push_back(m_HFClusU[ic]->getPos().z());
      m_HalfClusterU_E.push_back(m_HFClusU[ic]->getEnergy());
      m_HalfClusterU_type.push_back(m_HFClusU[ic]->getType());
      m_HalfClusterU_nTrk.push_back(m_HFClusU[ic]->getAssociatedTracks().size());
    }

    for(int ic=0; ic<m_NclusV; ic++){
      m_HalfClusterV_tag.push_back(it);
      m_HalfClusterV_x.push_back(m_HFClusV[ic]->getPos().x());
      m_HalfClusterV_y.push_back(m_HFClusV[ic]->getPos().y());
      m_HalfClusterV_z.push_back(m_HFClusV[ic]->getPos().z());
      m_HalfClusterV_E.push_back(m_HFClusV[ic]->getEnergy());
      m_HalfClusterV_type.push_back(m_HFClusV[ic]->getType());
      m_HalfClusterV_nTrk.push_back(m_HFClusV[ic]->getAssociatedTracks().size());
    }

    t_Tower->Fill();
  }


  //3D cluster
  ClearCluster();
  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> m_EcalClusterCol = m_DataCol.map_CaloCluster["EcalCluster"];
  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> m_HcalClusterCol = m_DataCol.map_CaloCluster["HCALCluster"];
  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> m_SimpleHcalClusterCol = m_DataCol.map_CaloCluster["SimpleHCALCluster"];
  m_totE_Ecal = 0.;
  m_totE_Hcal = 0.;
  m_Nclus_Ecal = m_EcalClusterCol.size();
  m_Nclus_Hcal = m_HcalClusterCol.size();
  for(int icl=0; icl<m_EcalClusterCol.size(); icl++){
    m_EcalClus_x.push_back(m_EcalClusterCol[icl]->getShowerCenter().x());
    m_EcalClus_y.push_back(m_EcalClusterCol[icl]->getShowerCenter().y());
    m_EcalClus_z.push_back(m_EcalClusterCol[icl]->getShowerCenter().z());
    m_EcalClus_E.push_back(m_EcalClusterCol[icl]->getLongiE());
    m_EcalClus_nTrk.push_back(m_EcalClusterCol[icl]->getAssociatedTracks().size());
    if(m_EcalClusterCol[icl]->getAssociatedTracks().size()==1)
      m_EcalClus_pTrk.push_back(m_EcalClusterCol[icl]->getAssociatedTracks()[0]->getMomentum());    
    else
      m_EcalClus_pTrk.push_back(-99);

    m_EcalClus_typeU.push_back(m_EcalClusterCol[icl]->getHalfClusterUCol("LinkedLongiCluster")[0]->getType());
    m_EcalClus_typeV.push_back(m_EcalClusterCol[icl]->getHalfClusterVCol("LinkedLongiCluster")[0]->getType());

    for(int ihit=0; ihit<m_EcalClusterCol[icl]->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("MergedAxis")[0]->getCluster().size(); ihit++){
      auto shower = m_EcalClusterCol[icl]->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("MergedAxis")[0]->getCluster()[ihit];
      m_EcalClus_hitU_tag.push_back(icl);
      m_EcalClus_hitU_x.push_back(shower->getPos().x());
      m_EcalClus_hitU_y.push_back(shower->getPos().y());
      m_EcalClus_hitU_z.push_back(shower->getPos().z());
      m_EcalClus_hitU_E.push_back(shower->getEnergy());
    }
    for(int ihit=0; ihit<m_EcalClusterCol[icl]->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("MergedAxis")[0]->getCluster().size(); ihit++){
      auto shower = m_EcalClusterCol[icl]->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("MergedAxis")[0]->getCluster()[ihit];
      m_EcalClus_hitV_tag.push_back(icl);
      m_EcalClus_hitV_x.push_back(shower->getPos().x());
      m_EcalClus_hitV_y.push_back(shower->getPos().y());
      m_EcalClus_hitV_z.push_back(shower->getPos().z());
      m_EcalClus_hitV_E.push_back(shower->getEnergy());
    }

    m_totE_Ecal += m_EcalClusterCol[icl]->getLongiE();
    auto truthMap = m_EcalClusterCol[icl]->getLinkedMCP();
    for(auto iter: truthMap){
      m_EcalClus_truthMC_tag.push_back(icl);
      m_EcalClus_truthMC_pid.push_back(iter.first.getPDG() );
      m_EcalClus_truthMC_px.push_back(iter.first.getMomentum().x);
      m_EcalClus_truthMC_py.push_back(iter.first.getMomentum().y);
      m_EcalClus_truthMC_pz.push_back(iter.first.getMomentum().z);
      m_EcalClus_truthMC_E.push_back(iter.first.getEnergy());
      m_EcalClus_truthMC_EPx.push_back(iter.first.getEndpoint().x);
      m_EcalClus_truthMC_EPy.push_back(iter.first.getEndpoint().y);
      m_EcalClus_truthMC_EPz.push_back(iter.first.getEndpoint().z);
      m_EcalClus_truthMC_weight.push_back(iter.second);
    }
  }
  for(int icl=0; icl<m_HcalClusterCol.size(); icl++){
    m_HcalClus_x.push_back(m_HcalClusterCol[icl]->getHitCenter().x());
    m_HcalClus_y.push_back(m_HcalClusterCol[icl]->getHitCenter().y());
    m_HcalClus_z.push_back(m_HcalClusterCol[icl]->getHitCenter().z());
    m_HcalClus_E.push_back(m_HcalClusterCol[icl]->getHitsE());
    m_HcalClus_nTrk.push_back(m_HcalClusterCol[icl]->getAssociatedTracks().size());
    if(m_HcalClusterCol[icl]->getAssociatedTracks().size()==1)
      m_HcalClus_pTrk.push_back(m_HcalClusterCol[icl]->getAssociatedTracks()[0]->getMomentum());
    else
      m_HcalClus_pTrk.push_back(-99);

    m_totE_Hcal += m_HcalClusterCol[icl]->getHitsE();
    auto truthMap = m_HcalClusterCol[icl]->getLinkedMCP();
    for(auto iter: truthMap){
      m_HcalClus_truthMC_tag.push_back(icl);
      m_HcalClus_truthMC_pid.push_back(iter.first.getPDG() );
      m_HcalClus_truthMC_px.push_back(iter.first.getMomentum().x);
      m_HcalClus_truthMC_py.push_back(iter.first.getMomentum().y);
      m_HcalClus_truthMC_pz.push_back(iter.first.getMomentum().z);
      m_HcalClus_truthMC_E.push_back(iter.first.getEnergy());
      m_HcalClus_truthMC_EPx.push_back(iter.first.getEndpoint().x);
      m_HcalClus_truthMC_EPy.push_back(iter.first.getEndpoint().y);
      m_HcalClus_truthMC_EPz.push_back(iter.first.getEndpoint().z);
      m_HcalClus_truthMC_weight.push_back(iter.second);
    }
  }
  for(int icl=0; icl<m_SimpleHcalClusterCol.size(); icl++){
    m_SimpleHcalClus_x.push_back(m_SimpleHcalClusterCol[icl]->getHitCenter().x());
    m_SimpleHcalClus_y.push_back(m_SimpleHcalClusterCol[icl]->getHitCenter().y());
    m_SimpleHcalClus_z.push_back(m_SimpleHcalClusterCol[icl]->getHitCenter().z());
    m_SimpleHcalClus_E.push_back(m_SimpleHcalClusterCol[icl]->getHitsE());
    m_SimpleHcalClus_nTrk.push_back(m_SimpleHcalClusterCol[icl]->getAssociatedTracks().size());
    if(m_SimpleHcalClusterCol[icl]->getAssociatedTracks().size()==1)
      m_SimpleHcalClus_pTrk.push_back(m_SimpleHcalClusterCol[icl]->getAssociatedTracks()[0]->getMomentum());
    else
      m_SimpleHcalClus_pTrk.push_back(-99);

    auto truthMap = m_SimpleHcalClusterCol[icl]->getLinkedMCP();
    for(auto iter: truthMap){
      m_SimpleHcalClus_truthMC_tag.push_back(icl);
      m_SimpleHcalClus_truthMC_pid.push_back(iter.first.getPDG() );
      m_SimpleHcalClus_truthMC_px.push_back(iter.first.getMomentum().x);
      m_SimpleHcalClus_truthMC_py.push_back(iter.first.getMomentum().y);
      m_SimpleHcalClus_truthMC_pz.push_back(iter.first.getMomentum().z);
      m_SimpleHcalClus_truthMC_E.push_back(iter.first.getEnergy());
      m_SimpleHcalClus_truthMC_EPx.push_back(iter.first.getEndpoint().x);
      m_SimpleHcalClus_truthMC_EPy.push_back(iter.first.getEndpoint().y);
      m_SimpleHcalClus_truthMC_EPz.push_back(iter.first.getEndpoint().z);
      m_SimpleHcalClus_truthMC_weight.push_back(iter.second);
    }
  }
  t_Cluster->Fill();

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
  t_MCParticle->Write();
  t_SimBar->Write();
  t_LocalMax->Write();
  t_Layers->Write();
  t_axis->Write();
  t_HalfCluster->Write();
  t_Tower->Write();
  t_Cluster->Write();
  m_wfile->Close();
  delete m_wfile, t_MCParticle, t_SimBar, t_LocalMax, t_Layers, t_axis, t_HalfCluster, t_Cluster; 

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

void PandoraPlusPFAlg::ClearMCParticle(){
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
  m_simBar_truthMC_tag.clear(); 
  m_simBar_truthMC_pid.clear();
  m_simBar_truthMC_px.clear();
  m_simBar_truthMC_py.clear();
  m_simBar_truthMC_pz.clear();
  m_simBar_truthMC_E.clear();
  m_simBar_truthMC_EPx.clear();
  m_simBar_truthMC_EPy.clear();
  m_simBar_truthMC_EPz.clear();
  m_simBar_truthMC_weight.clear();

  m_HcalHit_x.clear();
  m_HcalHit_y.clear();
  m_HcalHit_z.clear();
  m_HcalHit_E.clear();
  m_HcalHit_layer.clear();
  m_HcalHit_truthMC_tag.clear();
  m_HcalHit_truthMC_pid.clear();
  m_HcalHit_truthMC_px.clear();
  m_HcalHit_truthMC_py.clear();
  m_HcalHit_truthMC_pz.clear();
  m_HcalHit_truthMC_E.clear();
  m_HcalHit_truthMC_EPx.clear();
  m_HcalHit_truthMC_EPy.clear();
  m_HcalHit_truthMC_EPz.clear();
  m_HcalHit_truthMC_weight.clear();
}

void PandoraPlusPFAlg::ClearLocalMax(){
  m_NlmU=-99;
  m_NlmV=-99;
  m_localMaxU_tag.clear();
  m_localMaxU_x.clear();
  m_localMaxU_y.clear();
  m_localMaxU_z.clear();
  m_localMaxU_E.clear();
  m_localMaxU_mc_tag.clear();
  m_localMaxU_mc_pdg.clear();
  m_localMaxU_mc_px.clear();
  m_localMaxU_mc_py.clear();
  m_localMaxU_mc_pz.clear();
  m_localMaxU_mc_weight.clear();
  m_localMaxV_tag.clear();
  m_localMaxV_x.clear();
  m_localMaxV_y.clear();
  m_localMaxV_z.clear();
  m_localMaxV_E.clear();
  m_localMaxV_mc_tag.clear();
  m_localMaxV_mc_pdg.clear();
  m_localMaxV_mc_px.clear();
  m_localMaxV_mc_py.clear();
  m_localMaxV_mc_pz.clear();
  m_localMaxV_mc_weight.clear();
}

void PandoraPlusPFAlg::ClearLayer(){
  m_NshowerU=-99;
  m_NshowerV=-99;
  m_barShowerU_tag.clear();
  m_barShowerU_x.clear();
  m_barShowerU_y.clear();
  m_barShowerU_z.clear();
  m_barShowerU_E.clear();
  m_barShowerU_mc_tag.clear();
  m_barShowerU_mc_pdg.clear();
  m_barShowerU_mc_px.clear();
  m_barShowerU_mc_py.clear();
  m_barShowerU_mc_pz.clear();
  m_barShowerU_mc_weight.clear();
  m_barShowerV_tag.clear();
  m_barShowerV_x.clear();
  m_barShowerV_y.clear();
  m_barShowerV_z.clear();
  m_barShowerV_E.clear();
  m_barShowerV_mc_tag.clear();
  m_barShowerV_mc_pdg.clear();
  m_barShowerV_mc_px.clear();
  m_barShowerV_mc_py.clear();
  m_barShowerV_mc_pz.clear();
  m_barShowerV_mc_weight.clear();
}

void PandoraPlusPFAlg::ClearAxis(){
  m_axisU_tag.clear();
  m_axisU_type.clear();
  m_axisU_x.clear();
  m_axisU_y.clear();
  m_axisU_z.clear();
  m_axisU_E.clear();
  m_axisU_truth_tag.clear();
  m_axisU_truth_MC_px.clear();
  m_axisU_truth_MC_py.clear();
  m_axisU_truth_MC_pz.clear();
  m_axisU_truth_MC_E.clear();
  m_axisU_truth_MC_weight.clear();
  m_axisU_hit_tag.clear();
  m_axisU_hit_x.clear();
  m_axisU_hit_y.clear();
  m_axisU_hit_z.clear();
  m_axisU_hit_E.clear();
  m_axisV_tag.clear();
  m_axisV_type.clear();
  m_axisV_x.clear();
  m_axisV_y.clear();
  m_axisV_z.clear();
  m_axisV_E.clear();
  m_axisV_truth_tag.clear();
  m_axisV_truth_MC_px.clear();
  m_axisV_truth_MC_py.clear();
  m_axisV_truth_MC_pz.clear();
  m_axisV_truth_MC_E.clear();
  m_axisV_truth_MC_weight.clear();
  m_axisV_hit_tag.clear();
  m_axisV_hit_x.clear();
  m_axisV_hit_y.clear();
  m_axisV_hit_z.clear();
  m_axisV_hit_E.clear();

  m_emptyAxisU_tag.clear();
  m_emptyAxisU_x.clear();
  m_emptyAxisU_y.clear();
  m_emptyAxisU_z.clear();
  m_emptyAxisU_E.clear();
  m_emptyAxisV_tag.clear();
  m_emptyAxisV_x.clear();
  m_emptyAxisV_y.clear();
  m_emptyAxisV_z.clear();
  m_emptyAxisV_E.clear();
}

void PandoraPlusPFAlg::ClearHalfCluster(){
  m_totE_HFClusU = -99.;
  m_totE_HFClusV = -99.;
  m_HalfClusterV_x.clear();
  m_HalfClusterV_y.clear();
  m_HalfClusterV_z.clear();
  m_HalfClusterV_E.clear();
  m_HalfClusterV_tag.clear();
  m_HalfClusterV_type.clear();
  m_HalfClusterV_nTrk.clear();
  m_HalfClusterV_hit_x.clear();
  m_HalfClusterV_hit_y.clear();
  m_HalfClusterV_hit_z.clear();
  m_HalfClusterV_hit_E.clear();
  m_HalfClusterV_hit_tag.clear();
  m_HalfClusterV_truth_tag.clear();
  m_HalfClusterV_truthMC_px.clear();
  m_HalfClusterV_truthMC_py.clear();
  m_HalfClusterV_truthMC_pz.clear();
  m_HalfClusterV_truthMC_E.clear();
  m_HalfClusterV_truthMC_weight.clear();
  m_HalfClusterU_x.clear();
  m_HalfClusterU_y.clear();
  m_HalfClusterU_z.clear();
  m_HalfClusterU_E.clear();
  m_HalfClusterU_tag.clear();
  m_HalfClusterU_type.clear();
  m_HalfClusterU_nTrk.clear();
  m_HalfClusterU_hit_x.clear();
  m_HalfClusterU_hit_y.clear();
  m_HalfClusterU_hit_z.clear();
  m_HalfClusterU_hit_E.clear();
  m_HalfClusterU_hit_tag.clear();
  m_HalfClusterU_truth_tag.clear();
  m_HalfClusterU_truthMC_px.clear();
  m_HalfClusterU_truthMC_py.clear();
  m_HalfClusterU_truthMC_pz.clear();
  m_HalfClusterU_truthMC_E.clear();
  m_HalfClusterU_truthMC_weight.clear();
}

void PandoraPlusPFAlg::ClearTower(){
  towerID[0] = 0;
  towerID[1] = 0;
  towerID[2] = 0;
  m_NclusU = -99;
  m_NclusV = -99;
  m_totEn = -99;
  m_HalfClusterV_x.clear();
  m_HalfClusterV_y.clear();
  m_HalfClusterV_z.clear();
  m_HalfClusterV_E.clear();
  m_HalfClusterV_tag.clear();
  m_HalfClusterV_type.clear();
  m_HalfClusterV_nTrk.clear();
  m_HalfClusterU_x.clear();
  m_HalfClusterU_y.clear();
  m_HalfClusterU_z.clear();
  m_HalfClusterU_E.clear();
  m_HalfClusterU_tag.clear();
  m_HalfClusterU_type.clear();
  m_HalfClusterU_nTrk.clear();

}

void PandoraPlusPFAlg::ClearCluster(){
  m_totE_Ecal = -99.;
  m_totE_Hcal = -99.;
  m_Nclus_Ecal = -99;
  m_Nclus_Hcal = -99;
  m_EcalClus_x.clear();
  m_EcalClus_y.clear();
  m_EcalClus_z.clear();
  m_EcalClus_E.clear();
  m_EcalClus_nTrk.clear();
  m_EcalClus_pTrk.clear();
  m_EcalClus_typeU.clear();
  m_EcalClus_typeV.clear();
  m_EcalClus_hitU_x.clear();
  m_EcalClus_hitU_y.clear();
  m_EcalClus_hitU_z.clear();
  m_EcalClus_hitU_E.clear();
  m_EcalClus_hitU_tag.clear();
  m_EcalClus_hitV_x.clear();
  m_EcalClus_hitV_y.clear();
  m_EcalClus_hitV_z.clear();
  m_EcalClus_hitV_E.clear();
  m_EcalClus_hitV_tag.clear();
  m_EcalClus_truthMC_tag.clear();
  m_EcalClus_truthMC_pid.clear();
  m_EcalClus_truthMC_px.clear();
  m_EcalClus_truthMC_py.clear();
  m_EcalClus_truthMC_pz.clear();
  m_EcalClus_truthMC_E.clear();
  m_EcalClus_truthMC_EPx.clear();
  m_EcalClus_truthMC_EPy.clear();
  m_EcalClus_truthMC_EPz.clear();
  m_EcalClus_truthMC_weight.clear();
  m_HcalClus_x.clear();
  m_HcalClus_y.clear();
  m_HcalClus_z.clear();
  m_HcalClus_E.clear();
  m_HcalClus_nTrk.clear();
  m_HcalClus_pTrk.clear();
  m_HcalClus_truthMC_tag.clear();
  m_HcalClus_truthMC_pid.clear();
  m_HcalClus_truthMC_px.clear();
  m_HcalClus_truthMC_py.clear();
  m_HcalClus_truthMC_pz.clear();
  m_HcalClus_truthMC_E.clear();
  m_HcalClus_truthMC_EPx.clear();
  m_HcalClus_truthMC_EPy.clear();
  m_HcalClus_truthMC_EPz.clear();
  m_HcalClus_truthMC_weight.clear();
  m_SimpleHcalClus_x.clear();
  m_SimpleHcalClus_y.clear();
  m_SimpleHcalClus_z.clear();
  m_SimpleHcalClus_E.clear();
  m_SimpleHcalClus_nTrk.clear();
  m_SimpleHcalClus_pTrk.clear();
  m_SimpleHcalClus_truthMC_tag.clear();
  m_SimpleHcalClus_truthMC_pid.clear();
  m_SimpleHcalClus_truthMC_px.clear();
  m_SimpleHcalClus_truthMC_py.clear();
  m_SimpleHcalClus_truthMC_pz.clear();
  m_SimpleHcalClus_truthMC_E.clear();
  m_SimpleHcalClus_truthMC_EPx.clear();
  m_SimpleHcalClus_truthMC_EPy.clear();
  m_SimpleHcalClus_truthMC_EPz.clear();
  m_SimpleHcalClus_truthMC_weight.clear();


}


#endif
