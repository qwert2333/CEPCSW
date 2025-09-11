//=============================================================
// CyberPFA: a PFA developed for CEPC referenece detector
// Ver. CyberPFA-5.0.1(2025.01.09)
//-------------------------------------------------------------
//  Author: Fangyi Guo, Yang Zhang, Weizheng Song, Shengsen Sun
//          (IHEP, CAS)
//  Contact: guofangyi@ihep.ac.cn,
//           sunss@ihep.ac.cn
//=============================================================

#ifndef PANDORAPLUS_ALG_C
#define PANDORAPLUS_ALG_C

#include "CyberPFAlg.h"

// #include <fstream>
// #include <ctime>

using namespace std;
using namespace dd4hep;

int Cyber::CaloUnit::System_Barrel = 20;
int Cyber::CaloUnit::System_Endcap = 29;
int Cyber::CaloUnit::Nmodule = 32;
int Cyber::CaloUnit::Nstave = 15;
int Cyber::CaloUnit::Nlayer = 9;
int Cyber::CaloUnit::NbarPhi_odd[9] = {25, 25, 25, 23, 23, 23, 23, 23, 21};
int Cyber::CaloUnit::NbarPhi_even[9] = {19, 19, 21, 23, 23, 25, 27, 27, 29};
int Cyber::CaloUnit::NbarZ = 24;
//int Cyber::CaloUnit::over_module[28] = {13,15,16,18,19,21,22,24,25,26,28,29,30,32,33,35,36,38,39,41,42,43,45,46};
//int Cyber::CaloUnit::over_module_set = 2;
float Cyber::CaloUnit::barsize = 15.2; //mm
float Cyber::CaloUnit::ecal_innerR = 1830;  //mm
float Cyber::CaloUnit::ecal_endcap_deadarea = 10.5; //mm, a bit larger than real value 8.5 mm in geometry
float Cyber::CaloUnit::ecal_endcap_barsize = 16.2; //mm, a bit larger than real value 15.2 mm in geometry

DECLARE_COMPONENT( CyberPFAlg )

CyberPFAlg::CyberPFAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{
 
  // Output collections
  declareProperty("RecCaloHitCollection", w_RecEcalCol, "Handle of Reconstructed CaloHit collection");
  declareProperty("RecCaloCoreCollection", w_RecCoreCol, "Handle of Reconstructed ECAL Core collection");
  declareProperty("RecHCALHitCollection",  w_RecHcalCol, "Handle of Reconstructed HCAL CaloHit collection");
  declareProperty("RecTrackCollection", w_RecTrkCol, "Handle of Reconstructed Tracks linked to PFO"); 
  declareProperty("RecoPFOCollection", w_ReconstructedParticleCollection, "Handle of Reconstructed PFO collection");   
  // declareProperty("RecoVtxCollection", w_VertexCollection, "Handle of Reconstructed vertex collection");   
  // declareProperty("MCRecoPFOAssociationCollection", w_MCRecoParticleAssociationCollection, "Handle of MC-RecoPFO association collection");   

}

StatusCode CyberPFAlg::initialize()
{
  //Initialize global settings
  m_GlobalSettings.map_floatPars["BField"] = m_BField;
  m_GlobalSettings.map_floatPars["Seed"] = m_seed;
  m_GlobalSettings.map_intPars["Debug"] = m_Debug;
  m_GlobalSettings.map_intPars["SkipEvt"] = m_Nskip;


  //Initialize Creator settings
  m_pMCParticleCreatorSettings.map_stringPars["MCParticleCollections"] = name_MCParticleCol.value();

  m_pTrackCreatorSettings.map_stringVecPars["trackCollections"] = name_TrackCol.value();
  m_pTrackCreatorSettings.map_boolPars["DoCleanTrack"] = m_doCleanTrack.value();
  m_pTrackCreatorSettings.map_boolPars["UseTruthMatchTrk"] = m_useTruthMatchTrk.value();
  m_pTrackCreatorSettings.map_stringPars["TrackIDWeightFile"] = m_trackIDFile.value();
  m_pTrackCreatorSettings.map_stringPars["TrackIDMethod"] = m_trackIDMethod.value();
  m_pTrackCreatorSettings.map_floatPars["BField"] = m_BField; 
  m_pTrackCreatorSettings.map_floatPars["TrkMaxIP"] = 100.;
  m_pTrackCreatorSettings.map_floatPars["BDTCut"] = -0.8289;
  //m_pTrackCreatorSettings.map_floatPars["TrkEndZCut"] = 200.;
  //m_pTrackCreatorSettings.map_floatPars["TrkEndRCutMin"] = 700.;
  //m_pTrackCreatorSettings.map_floatPars["TrkEndRCutMax"] = 1600.;
  //m_pTrackCreatorSettings.map_floatPars["TrkStartRCutMin"] = 635.;
  //m_pTrackCreatorSettings.map_floatPars["TrkStartRCutMax"] = 640.;
  //m_pTrackCreatorSettings.map_floatPars["TrkLengthCut"] = 130.;
  //m_pTrackCreatorSettings.map_floatPars["BrokenTrkMinP"] = 0.5;
  //m_pTrackCreatorSettings.map_floatPars["BrokenTrkDeltaPCut"] = 0.15;
  //m_pTrackCreatorSettings.map_floatPars["BrokenTrkDistance"] = 10.;

  std::vector<std::string> name_CaloHits = name_EcalHits; 
  std::vector<std::string> name_CaloReadout = name_EcalReadout;
  name_CaloHits.insert( name_CaloHits.end(), name_HcalHits.begin(), name_HcalHits.end() );
  name_CaloReadout.insert(name_CaloReadout.end(), name_HcalReadout.begin(), name_HcalReadout.end());

  m_CaloHitsCreatorSettings.map_stringVecPars["CaloHitCollections"] = name_CaloHits;
  m_CaloHitsCreatorSettings.map_stringPars["EcalType"] = m_EcalType.value();

  m_OutputCreatorSettings.map_stringPars["OutputPFO"] = name_PFObject.value();
  m_OutputCreatorSettings.map_boolPars["UseTruthTrk"] = m_useMCPTrk.value();
  m_OutputCreatorSettings.map_boolPars["UseTruthMatchTrk"] = m_useTruthMatchTrk.value();
  m_OutputCreatorSettings.map_floatPars["BField"] = m_BField.value();
  m_OutputCreatorSettings.map_floatPars["ECALChargedCalib"] = m_EcalChargedCalib.value();
  m_OutputCreatorSettings.map_floatPars["HCALChargedCalib"] = m_HcalChargedCalib.value();
  m_OutputCreatorSettings.map_floatPars["ECALNeutralCalib"] = m_EcalNeutralCalib.value();
  m_OutputCreatorSettings.map_floatPars["HCALNeutralCalib"] = m_HcalNeutralCalib.value();

  //Initialize Creators
  m_pMCParticleCreator = new MCParticleCreator( m_pMCParticleCreatorSettings );
  m_pTrackCreator      = new TrackCreator( m_pTrackCreatorSettings );
  m_pCaloHitsCreator   = new CaloHitsCreator( m_CaloHitsCreatorSettings );
  m_pOutputCreator     = new OutputCreator( m_OutputCreatorSettings );

  //Readin collections

  //---MC particle---
  if(!name_MCParticleCol.empty()) r_MCParticleCol = new DataHandle<edm4hep::MCParticleCollection> (name_MCParticleCol, Gaudi::DataHandle::Reader, this);

  //---Tracks, dN/dx and TOF---
  for(auto& _trk : name_TrackCol) if(!_trk.empty()) r_TrackCols.push_back( new TrackType(_trk, Gaudi::DataHandle::Reader, this) );
  if(!name_dNdxCol.empty()) r_dNdxCol = new DataHandle<edm4hep::RecDqdxCollection> (name_dNdxCol, Gaudi::DataHandle::Reader, this);
  if(!name_tofCol.empty()) r_TofCol = new DataHandle<edm4hep::RecTofCollection>(name_tofCol, Gaudi::DataHandle::Reader, this);


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
  m_algorithmManager.RegisterAlgorithmFactory("ConeClusteringAlg",      new ConeClusteringAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("TrackClusterConnectingAlg",   new TrackClusterConnectingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("PFOReclusteringAlg",          new PFOReclusteringAlg::Factory);
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
      //if(type_AlgPars.value()[ialg].at(ipar)=="stringVec") m_settings.map_stringVecPars[name_AlgPars.value()[ialg].at(ipar)] = value_AlgPars.value()[ialg].at(ipar) ;
      if(type_AlgPars.value()[ialg].at(ipar)=="bool")   m_settings.map_boolPars[name_AlgPars.value()[ialg].at(ipar)] = (bool)std::stoi( (string)value_AlgPars.value()[ialg].at(ipar) );
    }

    m_algorithmManager.RegisterAlgorithm( name_Algs.value()[ialg], m_settings );
  }


  //Initialize services
  m_geosvc = service<IGeomSvc>("GeomSvc");
  if ( !m_geosvc )  throw "CyberPFAlg :Failed to find GeomSvc ...";

  m_dd4hep = m_geosvc->lcdd();
  if ( !m_dd4hep )  throw "CyberPFAlg :Failed to get dd4hep::Detector ...";
  
  //m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep);
  //m_volumeManager = m_dd4hep->volumeManager();

  dd4hep::rec::ECALSystemInfoData* EcalEndcapData = m_geosvc->getDD4HepGeo().child("EcalEndcap").extension<dd4hep::rec::ECALSystemInfoData>();
  
  for(int imodule=0; imodule<EcalEndcapData->ModuleInfos.size(); imodule++){
    dd4hep::rec::ECALModuleInfoStruct tmp_module = EcalEndcapData->ModuleInfos[imodule];
    for(int ilayer=0; ilayer<tmp_module.LayerInfos.size(); ilayer++){
      dd4hep::rec::ECALModuleInfoStruct::LayerInfo layer = tmp_module.LayerInfos[ilayer];
      std::tuple<int, int, int, int, int> tmp_key = std::make_tuple(tmp_module.moduleNumber, tmp_module.staveNumber, tmp_module.partNumber, layer.dlayerNumber, layer.slayerNumber);
      barNumberMapEndcapMap[tmp_key] = layer.barNumber;
    }
  }

  m_energycorsvc = service<ICrystalEcalSvc>("CrystalEcalEnergyCorrectionSvc");
  if ( !m_energycorsvc )  throw "CyberPFAlg :Failed to find CrystalEcalEnergyCorrectionSvc ...";
  //m_energycorsvc->initialize();

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
  std::cout<<"CyberPFAlg::initialize"<<std::endl;


  //Output collections
  w_ClusterCollection["EcalCluster"] = new DataHandle<edm4hep::ClusterCollection>("EcalCluster", Gaudi::DataHandle::Writer, this);
  w_ClusterCollection["EcalCore"]    = new DataHandle<edm4hep::ClusterCollection>("EcalCore", Gaudi::DataHandle::Writer, this);
  w_ClusterCollection["HcalCluster"] = new DataHandle<edm4hep::ClusterCollection>("HcalCluster", Gaudi::DataHandle::Writer, this);

  //Output ntuple for analysis.
  if(m_WriteAna){
    std::string s_outfile = m_filename;
    m_wfile = new TFile(s_outfile.c_str(), "recreate");
    t_MCParticle = new TTree("MCParticle", "MCParticle");
    t_SimBar = new TTree("SimBarHit", "SimBarHit");
    t_LocalMax = new TTree("LocalMax", "LocalMax");
    t_Layers = new TTree("RecLayers","RecLayers");
    t_Hough = new TTree("Hough", "Hough");
    t_Cone = new TTree("Cone", "Cone");
    t_TrackAxis = new TTree("TrackAxis", "TrackAxis");
    t_Axis = new TTree("Axis", "Axis");
    t_HalfCluster = new TTree("HalfCluster","HalfCluster");
    t_Tower = new TTree("Tower", "Tower");
    t_Cluster = new TTree("RecClusters", "RecClusters");
    t_Track = new TTree("Track", "Track");
    t_PFO = new TTree("PFO", "PFO");

    //MC particle 
    t_MCParticle->Branch("mcPdgid",     &m_mcPdgid);
    t_MCParticle->Branch("mcStatus",    &m_mcStatus);
    t_MCParticle->Branch("mcPx", &m_mcPx);
    t_MCParticle->Branch("mcPy", &m_mcPy);
    t_MCParticle->Branch("mcPz", &m_mcPz);
    t_MCParticle->Branch("mcEn", &m_mcEn);
    t_MCParticle->Branch("mcMass", &m_mcMass);
    t_MCParticle->Branch("mcCharge", &m_mcCharge);
    t_MCParticle->Branch("mcVTXx", &m_mcVTXx);
    t_MCParticle->Branch("mcVTXy", &m_mcVTXy);
    t_MCParticle->Branch("mcVTXz", &m_mcVTXz);
    t_MCParticle->Branch("mcEPx", &m_mcEPx);
    t_MCParticle->Branch("mcEPy", &m_mcEPy);
    t_MCParticle->Branch("mcEPz", &m_mcEPz);
    t_MCParticle->Branch("mcdepEn_ecal", &m_depEn_ecal);
    t_MCParticle->Branch("mcdepEn_hcal", &m_depEn_hcal);

    //Bar
    t_SimBar->Branch("totE_EcalSim", &m_totE_EcalSim);
    t_SimBar->Branch("simBar_x", &m_simBar_x);
    t_SimBar->Branch("simBar_y", &m_simBar_y);
    t_SimBar->Branch("simBar_z", &m_simBar_z);
    t_SimBar->Branch("simBar_length", &m_simBar_length);
    t_SimBar->Branch("simBar_nBarInLayer", &m_simBar_nBarInLayer);
    t_SimBar->Branch("simBar_T1", &m_simBar_T1);
    t_SimBar->Branch("simBar_T2", &m_simBar_T2);
    t_SimBar->Branch("simBar_Q1", &m_simBar_Q1);
    t_SimBar->Branch("simBar_Q2", &m_simBar_Q2);
    t_SimBar->Branch("simBar_module", &m_simBar_module);
    t_SimBar->Branch("simBar_dlayer", &m_simBar_dlayer);
    t_SimBar->Branch("simBar_stave", &m_simBar_stave);
    t_SimBar->Branch("simBar_slayer", &m_simBar_slayer);
    t_SimBar->Branch("simBar_bar", &m_simBar_bar);
    t_SimBar->Branch("simBar_system", &m_simBar_system);
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

    // Hough
    t_Hough->Branch("houghU_tag", &m_houghU_tag);
    t_Hough->Branch("houghU_type", &m_houghU_type);
    t_Hough->Branch("houghU_x", &m_houghU_x);
    t_Hough->Branch("houghU_y", &m_houghU_y);
    t_Hough->Branch("houghU_z", &m_houghU_z);
    t_Hough->Branch("houghU_E", &m_houghU_E);
    t_Hough->Branch("houghU_truth_tag", &m_houghU_truth_tag);
    t_Hough->Branch("houghU_truth_MC_px", &m_houghU_truth_MC_px);
    t_Hough->Branch("houghU_truth_MC_py", &m_houghU_truth_MC_py);
    t_Hough->Branch("houghU_truth_MC_pz", &m_houghU_truth_MC_pz);
    t_Hough->Branch("houghU_truth_MC_E", &m_houghU_truth_MC_E);
    t_Hough->Branch("houghU_truth_MC_weight", &m_houghU_truth_MC_weight);    
    t_Hough->Branch("houghU_hit_tag", &m_houghU_hit_tag);
    t_Hough->Branch("houghU_hit_x", &m_houghU_hit_x);
    t_Hough->Branch("houghU_hit_y", &m_houghU_hit_y);
    t_Hough->Branch("houghU_hit_z", &m_houghU_hit_z);
    t_Hough->Branch("houghU_hit_E", &m_houghU_hit_E);
    t_Hough->Branch("houghV_tag", &m_houghV_tag);
    t_Hough->Branch("houghV_type", &m_houghV_type);
    t_Hough->Branch("houghV_x", &m_houghV_x);
    t_Hough->Branch("houghV_y", &m_houghV_y);
    t_Hough->Branch("houghV_z", &m_houghV_z);
    t_Hough->Branch("houghV_E", &m_houghV_E);
    t_Hough->Branch("houghV_alpha", &m_houghV_alpha);
    t_Hough->Branch("houghV_rho", &m_houghV_rho);
    t_Hough->Branch("houghV_truth_tag", &m_houghV_truth_tag);
    t_Hough->Branch("houghV_truth_MC_px", &m_houghV_truth_MC_px);
    t_Hough->Branch("houghV_truth_MC_py", &m_houghV_truth_MC_py);
    t_Hough->Branch("houghV_truth_MC_pz", &m_houghV_truth_MC_pz);
    t_Hough->Branch("houghV_truth_MC_E", &m_houghV_truth_MC_E);
    t_Hough->Branch("houghV_truth_MC_weight", &m_houghV_truth_MC_weight);   
    t_Hough->Branch("houghV_hit_tag", &m_houghV_hit_tag);
    t_Hough->Branch("houghV_hit_x", &m_houghV_hit_x);
    t_Hough->Branch("houghV_hit_y", &m_houghV_hit_y);
    t_Hough->Branch("houghV_hit_z", &m_houghV_hit_z);
    t_Hough->Branch("houghV_hit_E", &m_houghV_hit_E);
    // Cone
    t_Cone->Branch("coneU_tag", &m_coneU_tag);
    t_Cone->Branch("coneU_type", &m_coneU_type);
    t_Cone->Branch("coneU_x", &m_coneU_x);
    t_Cone->Branch("coneU_y", &m_coneU_y);
    t_Cone->Branch("coneU_z", &m_coneU_z);
    t_Cone->Branch("coneU_E", &m_coneU_E);
    t_Cone->Branch("coneU_truth_tag", &m_coneU_truth_tag);
    t_Cone->Branch("coneU_truth_MC_px", &m_coneU_truth_MC_px);
    t_Cone->Branch("coneU_truth_MC_py", &m_coneU_truth_MC_py);
    t_Cone->Branch("coneU_truth_MC_pz", &m_coneU_truth_MC_pz);
    t_Cone->Branch("coneU_truth_MC_E", &m_coneU_truth_MC_E);
    t_Cone->Branch("coneU_truth_MC_weight", &m_coneU_truth_MC_weight);    
    t_Cone->Branch("coneU_hit_tag", &m_coneU_hit_tag);
    t_Cone->Branch("coneU_hit_x", &m_coneU_hit_x);
    t_Cone->Branch("coneU_hit_y", &m_coneU_hit_y);
    t_Cone->Branch("coneU_hit_z", &m_coneU_hit_z);
    t_Cone->Branch("coneU_hit_E", &m_coneU_hit_E);
    t_Cone->Branch("coneV_tag", &m_coneV_tag);
    t_Cone->Branch("coneV_type", &m_coneV_type);
    t_Cone->Branch("coneV_x", &m_coneV_x);
    t_Cone->Branch("coneV_y", &m_coneV_y);
    t_Cone->Branch("coneV_z", &m_coneV_z);
    t_Cone->Branch("coneV_E", &m_coneV_E);
    t_Cone->Branch("coneV_truth_tag", &m_coneV_truth_tag);
    t_Cone->Branch("coneV_truth_MC_px", &m_coneV_truth_MC_px);
    t_Cone->Branch("coneV_truth_MC_py", &m_coneV_truth_MC_py);
    t_Cone->Branch("coneV_truth_MC_pz", &m_coneV_truth_MC_pz);
    t_Cone->Branch("coneV_truth_MC_E", &m_coneV_truth_MC_E);
    t_Cone->Branch("coneV_truth_MC_weight", &m_coneV_truth_MC_weight);   
    t_Cone->Branch("coneV_hit_tag", &m_coneV_hit_tag);
    t_Cone->Branch("coneV_hit_x", &m_coneV_hit_x);
    t_Cone->Branch("coneV_hit_y", &m_coneV_hit_y);
    t_Cone->Branch("coneV_hit_z", &m_coneV_hit_z);
    t_Cone->Branch("coneV_hit_E", &m_coneV_hit_E);
    // Track Axis
    t_TrackAxis->Branch("trackU_tag", &m_trackU_tag);
    t_TrackAxis->Branch("trackU_type", &m_trackU_type);
    t_TrackAxis->Branch("trackU_x", &m_trackU_x);
    t_TrackAxis->Branch("trackU_y", &m_trackU_y);
    t_TrackAxis->Branch("trackU_z", &m_trackU_z);
    t_TrackAxis->Branch("trackU_E", &m_trackU_E);
    t_TrackAxis->Branch("trackU_truth_tag", &m_trackU_truth_tag);
    t_TrackAxis->Branch("trackU_truth_MC_px", &m_trackU_truth_MC_px);
    t_TrackAxis->Branch("trackU_truth_MC_py", &m_trackU_truth_MC_py);
    t_TrackAxis->Branch("trackU_truth_MC_pz", &m_trackU_truth_MC_pz);
    t_TrackAxis->Branch("trackU_truth_MC_E", &m_trackU_truth_MC_E);
    t_TrackAxis->Branch("trackU_truth_MC_weight", &m_trackU_truth_MC_weight);    
    t_TrackAxis->Branch("trackU_hit_tag", &m_trackU_hit_tag);
    t_TrackAxis->Branch("trackU_hit_x", &m_trackU_hit_x);
    t_TrackAxis->Branch("trackU_hit_y", &m_trackU_hit_y);
    t_TrackAxis->Branch("trackU_hit_z", &m_trackU_hit_z);
    t_TrackAxis->Branch("trackU_hit_E", &m_trackU_hit_E);
    t_TrackAxis->Branch("trackV_tag", &m_trackV_tag);
    t_TrackAxis->Branch("trackV_type", &m_trackV_type);
    t_TrackAxis->Branch("trackV_x", &m_trackV_x);
    t_TrackAxis->Branch("trackV_y", &m_trackV_y);
    t_TrackAxis->Branch("trackV_z", &m_trackV_z);
    t_TrackAxis->Branch("trackV_E", &m_trackV_E);
    t_TrackAxis->Branch("trackV_truth_tag", &m_trackV_truth_tag);
    t_TrackAxis->Branch("trackV_truth_MC_px", &m_trackV_truth_MC_px);
    t_TrackAxis->Branch("trackV_truth_MC_py", &m_trackV_truth_MC_py);
    t_TrackAxis->Branch("trackV_truth_MC_pz", &m_trackV_truth_MC_pz);
    t_TrackAxis->Branch("trackV_truth_MC_E", &m_trackV_truth_MC_E);
    t_TrackAxis->Branch("trackV_truth_MC_weight", &m_trackV_truth_MC_weight);   
    t_TrackAxis->Branch("trackV_hit_tag", &m_trackV_hit_tag);
    t_TrackAxis->Branch("trackV_hit_x", &m_trackV_hit_x);
    t_TrackAxis->Branch("trackV_hit_y", &m_trackV_hit_y);
    t_TrackAxis->Branch("trackV_hit_z", &m_trackV_hit_z);
    t_TrackAxis->Branch("trackV_hit_E", &m_trackV_hit_E);
    //Axis
    t_Axis->Branch("axisU_tag", &m_axisU_tag);
    t_Axis->Branch("axisU_type", &m_axisU_type);
    t_Axis->Branch("axisU_x", &m_axisU_x);
    t_Axis->Branch("axisU_y", &m_axisU_y);
    t_Axis->Branch("axisU_z", &m_axisU_z);
    t_Axis->Branch("axisU_E", &m_axisU_E);
    t_Axis->Branch("axisU_truth_tag", &m_axisU_truth_tag);
    t_Axis->Branch("axisU_truth_MC_px", &m_axisU_truth_MC_px);
    t_Axis->Branch("axisU_truth_MC_py", &m_axisU_truth_MC_py);
    t_Axis->Branch("axisU_truth_MC_pz", &m_axisU_truth_MC_pz);
    t_Axis->Branch("axisU_truth_MC_E", &m_axisU_truth_MC_E);
    t_Axis->Branch("axisU_truth_MC_weight", &m_axisU_truth_MC_weight);
    t_Axis->Branch("axisU_hit_tag", &m_axisU_hit_tag);
    t_Axis->Branch("axisU_hit_x", &m_axisU_hit_x);
    t_Axis->Branch("axisU_hit_y", &m_axisU_hit_y);
    t_Axis->Branch("axisU_hit_z", &m_axisU_hit_z);
    t_Axis->Branch("axisU_hit_E", &m_axisU_hit_E);
    t_Axis->Branch("axisV_tag", &m_axisV_tag);
    t_Axis->Branch("axisV_type", &m_axisV_type);
    t_Axis->Branch("axisV_x", &m_axisV_x);
    t_Axis->Branch("axisV_y", &m_axisV_y);
    t_Axis->Branch("axisV_z", &m_axisV_z);
    t_Axis->Branch("axisV_E", &m_axisV_E);
    t_Axis->Branch("axisV_truth_tag", &m_axisV_truth_tag);
    t_Axis->Branch("axisV_truth_MC_px", &m_axisV_truth_MC_px);
    t_Axis->Branch("axisV_truth_MC_py", &m_axisV_truth_MC_py);
    t_Axis->Branch("axisV_truth_MC_pz", &m_axisV_truth_MC_pz);
    t_Axis->Branch("axisV_truth_MC_E", &m_axisV_truth_MC_E);
    t_Axis->Branch("axisV_truth_MC_weight", &m_axisV_truth_MC_weight);
    t_Axis->Branch("axisV_hit_tag", &m_axisV_hit_tag);
    t_Axis->Branch("axisV_hit_x", &m_axisV_hit_x);
    t_Axis->Branch("axisV_hit_y", &m_axisV_hit_y);
    t_Axis->Branch("axisV_hit_z", &m_axisV_hit_z);
    t_Axis->Branch("axisV_hit_E", &m_axisV_hit_E);    

    t_Axis->Branch("emptyAxisU_tag", &m_emptyAxisU_tag);
    t_Axis->Branch("emptyAxisU_x", &m_emptyAxisU_x);
    t_Axis->Branch("emptyAxisU_y", &m_emptyAxisU_y);
    t_Axis->Branch("emptyAxisU_z", &m_emptyAxisU_z);
    t_Axis->Branch("emptyAxisU_E", &m_emptyAxisU_E);
    t_Axis->Branch("emptyAxisV_tag", &m_emptyAxisV_tag);
    t_Axis->Branch("emptyAxisV_x", &m_emptyAxisV_x);
    t_Axis->Branch("emptyAxisV_y", &m_emptyAxisV_y);
    t_Axis->Branch("emptyAxisV_z", &m_emptyAxisV_z);
    t_Axis->Branch("emptyAxisV_E", &m_emptyAxisV_E);

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
    t_Tower->Branch("Ntower", &m_Ntower);
    t_Tower->Branch("towerID_id1", &m_towerID_id1);
    t_Tower->Branch("towerID_id2", &m_towerID_id2);
    t_Tower->Branch("NclusU", &m_NclusU);
    t_Tower->Branch("NclusV", &m_NclusV);
    t_Tower->Branch("totEn", &m_totEn);
    t_Tower->Branch("totEn_U", &m_totEn_U);
    t_Tower->Branch("totEn_V", &m_totEn_V);
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
    t_Cluster->Branch("EcalClus_Escale", &m_EcalClus_Escale);
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
    t_Cluster->Branch("EcalClus_trk_location", &m_EcalClus_trk_location);
    t_Cluster->Branch("EcalClus_trk_tag", &m_EcalClus_trk_tag);
    t_Cluster->Branch("EcalClus_trk_d0", &m_EcalClus_trk_d0);
    t_Cluster->Branch("EcalClus_trk_z0", &m_EcalClus_trk_z0);
    t_Cluster->Branch("EcalClus_trk_phi", &m_EcalClus_trk_phi);
    t_Cluster->Branch("EcalClus_trk_tanL", &m_EcalClus_trk_tanL);
    t_Cluster->Branch("EcalClus_trk_omega", &m_EcalClus_trk_omega);
    t_Cluster->Branch("EcalClus_trk_kappa", &m_EcalClus_trk_kappa);
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
    t_Cluster->Branch("HcalClus_nHit", &m_HcalClus_nHit);
    t_Cluster->Branch("HcalClus_nTrk", &m_HcalClus_nTrk);
    t_Cluster->Branch("HcalClus_ptrk", &m_HcalClus_pTrk);
    t_Cluster->Branch("HcalClus_hit_x", &m_HcalClus_hit_x);
    t_Cluster->Branch("HcalClus_hit_y", &m_HcalClus_hit_y);
    t_Cluster->Branch("HcalClus_hit_z", &m_HcalClus_hit_z);
    t_Cluster->Branch("HcalClus_hit_E", &m_HcalClus_hit_E);
    t_Cluster->Branch("HcalClus_hit_tag", &m_HcalClus_hit_tag);
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
    t_Cluster->Branch("SimpleHcalClus_nHit", &m_SimpleHcalClus_nHit);
    t_Cluster->Branch("SimpleHcalClus_nTrk", &m_SimpleHcalClus_nTrk);
    t_Cluster->Branch("SimpleHcalClus_ptrk", &m_SimpleHcalClus_pTrk);
    t_Cluster->Branch("SimpleHcalClus_hit_x", &m_SimpleHcalClus_hit_x);
    t_Cluster->Branch("SimpleHcalClus_hit_y", &m_SimpleHcalClus_hit_y);
    t_Cluster->Branch("SimpleHcalClus_hit_z", &m_SimpleHcalClus_hit_z);
    t_Cluster->Branch("SimpleHcalClus_hit_E", &m_SimpleHcalClus_hit_E);
    t_Cluster->Branch("SimpleHcalClus_hit_tag", &m_SimpleHcalClus_hit_tag);
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

    // Tracks
    t_Track->Branch("m_Ntrk", &m_Ntrk);
    t_Track->Branch("m_type", &m_type);
    t_Track->Branch("m_Nhit", &m_Nhit);
    t_Track->Branch("m_pid", &m_pid);
    t_Track->Branch("m_pid_truth", &m_pid_truth);
    t_Track->Branch("m_trk_px", &m_trk_px);
    t_Track->Branch("m_trk_py", &m_trk_py);
    t_Track->Branch("m_trk_pz", &m_trk_pz);
    t_Track->Branch("m_trk_p", &m_trk_p);
    t_Track->Branch("m_trk_truthweight", &m_trk_truthweight);
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
    t_Track->Branch("m_trkstate_x_ECAL", &m_trkstate_x_ECAL);
    t_Track->Branch("m_trkstate_y_ECAL", &m_trkstate_y_ECAL);
    t_Track->Branch("m_trkstate_z_ECAL", &m_trkstate_z_ECAL);
    t_Track->Branch("m_trkstate_tag_ECAL", &m_trkstate_tag_ECAL);
    t_Track->Branch("m_trkstate_x_HCAL", &m_trkstate_x_HCAL);
    t_Track->Branch("m_trkstate_y_HCAL", &m_trkstate_y_HCAL);
    t_Track->Branch("m_trkstate_z_HCAL", &m_trkstate_z_HCAL);
    t_Track->Branch("m_trkstate_tag_HCAL", &m_trkstate_tag_HCAL);

    //PFOs
    t_PFO->Branch("pfo_tag", &pfo_tag);
    t_PFO->Branch("pfo_n_track", &pfo_n_track);
    t_PFO->Branch("pfo_n_ecal_clus", &pfo_n_ecal_clus);
    t_PFO->Branch("pfo_n_hcal_clus", &pfo_n_hcal_clus);
    t_PFO->Branch("pfo_ecal_tag", &pfo_ecal_tag);
    t_PFO->Branch("pfo_hcal_tag", &pfo_hcal_tag);
    t_PFO->Branch("pfo_ecal_clus_x", &pfo_ecal_clus_x);
    t_PFO->Branch("pfo_ecal_clus_y", &pfo_ecal_clus_y);
    t_PFO->Branch("pfo_ecal_clus_z", &pfo_ecal_clus_z);
    t_PFO->Branch("pfo_ecal_clus_E", &pfo_ecal_clus_E);
    t_PFO->Branch("pfo_ecal_clus_Escale", &pfo_ecal_clus_Escale);
    t_PFO->Branch("pfo_hcal_clus_x", &pfo_hcal_clus_x);
    t_PFO->Branch("pfo_hcal_clus_y", &pfo_hcal_clus_y);
    t_PFO->Branch("pfo_hcal_clus_z", &pfo_hcal_clus_z);
    t_PFO->Branch("pfo_hcal_clus_E", &pfo_hcal_clus_E);
    t_PFO->Branch("pfo_trk_tag", &pfo_trk_tag);
    t_PFO->Branch("pfo_trk_location", &pfo_trk_location);
    t_PFO->Branch("pfo_trk_d0", &pfo_trk_d0);
    t_PFO->Branch("pfo_trk_z0", &pfo_trk_z0);
    t_PFO->Branch("pfo_trk_phi", &pfo_trk_phi);
    t_PFO->Branch("pfo_trk_tanL", &pfo_trk_tanL);
    t_PFO->Branch("pfo_trk_omega", &pfo_trk_omega);
    t_PFO->Branch("pfo_trk_kappa", &pfo_trk_kappa);    

  }

  return GaudiAlgorithm::initialize();
}

StatusCode CyberPFAlg::execute()
{
// clock_t yyy_start, yyy_endrec, yyy_endfill;
// yyy_start = clock(); // 记录开始时间

  if(_nEvt==0) std::cout<<"CyberPFAlg::execute Start"<<std::endl;
  std::cout<<"Processing event: "<<_nEvt<<std::endl;

  if(_nEvt<m_Nskip){ _nEvt++;  return StatusCode::SUCCESS; }

  //InitializeForNewEvent(); 
  CyberDataCol     m_DataCol;
  m_DataCol.Clear();
  m_DataCol.EnergyCorrSvc = m_energycorsvc; 
  m_DataCol.tofCol = const_cast<edm4hep::RecTofCollection*>(r_TofCol->get());
  m_DataCol.dNdxCol = const_cast<edm4hep::RecDqdxCollection*>(r_dNdxCol->get()); 


  //Readin collections 
  m_pMCParticleCreator->CreateMCParticle( m_DataCol, *r_MCParticleCol );
  if(m_useMCPTrk) m_pTrackCreator->CreateTracksFromMCParticle(m_DataCol, *r_MCParticleCol);
  else m_pTrackCreator->CreateTracks( m_DataCol, r_TrackCols, r_MCPTrkAssoCol );
  m_pCaloHitsCreator->CreateCaloHits( m_DataCol, r_CaloHitCols, map_readout_decoder, map_CaloMCPAssoCols, m_geosvc, barNumberMapEndcapMap);

  //Perform PFA algorithm
  m_algorithmManager.RunAlgorithm( m_DataCol );

  m_pOutputCreator->CreateOutputCollections( m_DataCol, 
                                             w_RecEcalCol, 
                                             w_RecCoreCol, 
                                             w_RecHcalCol,
                                             w_RecTrkCol, 
                                             w_ClusterCollection, 
                                             w_ReconstructedParticleCollection);

// yyy_endrec = clock();  // 重建结束的时间

  if(m_WriteAna){
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
      m_mcVTXx.push_back( m_MCPCol[imc].getVertex()[0] );
      m_mcVTXy.push_back( m_MCPCol[imc].getVertex()[1] );
      m_mcVTXz.push_back( m_MCPCol[imc].getVertex()[2] );
      m_mcEPx.push_back( m_MCPCol[imc].getEndpoint()[0] );
      m_mcEPy.push_back( m_MCPCol[imc].getEndpoint()[1] );
      m_mcEPz.push_back( m_MCPCol[imc].getEndpoint()[2] );
      //double tmp_phi = std::atan2(m_MCPCol[imc].getMomentum()[1], m_MCPCol[imc].getMomentum()[0])* 180.0 / M_PI;
      //if (tmp_phi < 0) tmp_phi += 360.0;
      //double tmp_theta = std::atan2(m_MCPCol[imc].getMomentum()[2], sqrt(m_MCPCol[imc].getMomentum()[1]*m_MCPCol[imc].getMomentum()[1]+m_MCPCol[imc].getMomentum()[0]*m_MCPCol[imc].getMomentum()[0]))* 180.0 / M_PI + 90; 
      //cout<<"MCParticle: "<<imc<<" PDG: "<<m_MCPCol[imc].getPDG()<<" Theta: "<<tmp_theta<<" Phi: "<<tmp_phi<<endl;
   
      double EnDep_ecal = GetParticleDepEnergy(m_MCPCol[imc], m_DataCol.map_BarCol["BarCol"]);
      double EnDep_hcal = GetParticleDepEnergy(m_MCPCol[imc], m_DataCol.map_CaloHit["HCALBarrel"]);
      m_depEn_ecal.push_back(EnDep_ecal);
      m_depEn_hcal.push_back(EnDep_hcal);
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
      m_simBar_length.push_back(p_hitbar->getBarLength());
      m_simBar_nBarInLayer.push_back(p_hitbar->getNBarInLayer());
      m_simBar_Q1.push_back(p_hitbar->getQ1());
      m_simBar_Q2.push_back(p_hitbar->getQ2());
      m_simBar_T1.push_back(p_hitbar->getT1());
      m_simBar_T2.push_back(p_hitbar->getT2());
      m_simBar_module.push_back(p_hitbar->getModule());
      m_simBar_dlayer.push_back(p_hitbar->getDlayer());
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
   
    std::vector<Cyber::CaloHit*> m_hcalHitsCol; m_hcalHitsCol.clear();
    for(int ih=0; ih<m_DataCol.map_CaloHit["HCALBarrel"].size(); ih++)
      m_hcalHitsCol.push_back( m_DataCol.map_CaloHit["HCALBarrel"][ih].get() );
    for(int ih=0; ih<m_DataCol.map_CaloHit["HCALEndcaps"].size(); ih++)
      m_hcalHitsCol.push_back( m_DataCol.map_CaloHit["HCALEndcaps"][ih].get() );   

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
    std::vector<Cyber::CaloHalfCluster*> m_halfclusters; m_halfclusters.clear();
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
   
    std::vector<const Cyber::CaloHalfCluster*> m_halfclusterV; m_halfclusterV.clear();
    std::vector<const Cyber::CaloHalfCluster*> m_halfclusterU; m_halfclusterU.clear();
    for(int i=0; i<m_DataCol.map_HalfCluster["HalfClusterColU"].size(); i++){
      m_halfclusterU.push_back( m_DataCol.map_HalfCluster["HalfClusterColU"][i].get() );
    }
    for(int i=0; i<m_DataCol.map_HalfCluster["HalfClusterColV"].size(); i++){
      m_halfclusterV.push_back( m_DataCol.map_HalfCluster["HalfClusterColV"][i].get() );
    }
    // Hough
    ClearHough();
    int houghU_index=0;
    int houghV_index=0;
    for(int i=0; i<m_halfclusterU.size(); i++){  // loop half cluster U
      std::vector<const Cyber::CaloHalfCluster*> m_mergedaxisU = m_halfclusterU[i]->getHalfClusterCol("HoughAxis");
      for(int ita=0; ita<m_mergedaxisU.size(); ita++){ // loop  axis U
        // General information of the axis
        m_houghU_tag.push_back(houghU_index);
        m_houghU_type.push_back(m_mergedaxisU[ita]->getType());
        m_houghU_x.push_back(m_mergedaxisU[ita]->getPos().x());
        m_houghU_y.push_back(m_mergedaxisU[ita]->getPos().y());
        m_houghU_z.push_back(m_mergedaxisU[ita]->getPos().z());
        m_houghU_E.push_back(m_mergedaxisU[ita]->getEnergy());
   
        // MC truth information of the Axis
        auto truthMap = m_mergedaxisU[ita]->getLinkedMCP();
        for(auto iter: truthMap){
          m_houghU_truth_tag.push_back(houghU_index);
          m_houghU_truth_MC_px.push_back(iter.first.getMomentum().x);
          m_houghU_truth_MC_py.push_back(iter.first.getMomentum().y);
          m_houghU_truth_MC_pz.push_back(iter.first.getMomentum().z);
          m_houghU_truth_MC_E.push_back(iter.first.getEnergy());
          m_houghU_truth_MC_weight.push_back(iter.second);
        }
   
        // Hits on axis
        for(int ilm=0; ilm<m_mergedaxisU[ita]->getCluster().size(); ilm++){ // loop local max
          m_houghU_hit_tag.push_back(houghU_index);
          m_houghU_hit_x.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getPos().x() );
          m_houghU_hit_y.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getPos().y() );
          m_houghU_hit_z.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getPos().z() );
          m_houghU_hit_E.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getEnergy()  );
        }
   
        houghU_index++;
      }
    }
    for(int i=0; i<m_halfclusterV.size(); i++){  // loop half cluster V
      std::vector<const Cyber::CaloHalfCluster*> m_mergedaxisV = m_halfclusterV[i]->getHalfClusterCol("HoughAxis");
      for(int ita=0; ita<m_mergedaxisV.size(); ita++){ // loop  axis V
        // General information of the axis
        m_houghV_tag.push_back(houghV_index);
        m_houghV_type.push_back(m_mergedaxisV[ita]->getType());
        m_houghV_x.push_back(m_mergedaxisV[ita]->getPos().x());
        m_houghV_y.push_back(m_mergedaxisV[ita]->getPos().y());
        m_houghV_z.push_back(m_mergedaxisV[ita]->getPos().z());
        m_houghV_E.push_back(m_mergedaxisV[ita]->getEnergy());
        m_houghV_alpha.push_back(m_mergedaxisV[ita]->getHoughAlpha());
        m_houghV_rho.push_back(m_mergedaxisV[ita]->getHoughRho());
   
        // MC truth information of the Axis
        auto truthMap = m_mergedaxisV[ita]->getLinkedMCP();
        for(auto iter: truthMap){
          m_houghV_truth_tag.push_back(houghV_index);
          m_houghV_truth_MC_px.push_back(iter.first.getMomentum().x);
          m_houghV_truth_MC_py.push_back(iter.first.getMomentum().y);
          m_houghV_truth_MC_pz.push_back(iter.first.getMomentum().z);
          m_houghV_truth_MC_E.push_back(iter.first.getEnergy());
          m_houghV_truth_MC_weight.push_back(iter.second);
        }
   
        // Hits on axis
        for(int ilm=0; ilm<m_mergedaxisV[ita]->getCluster().size(); ilm++){
          m_houghV_hit_tag.push_back(houghV_index);
          m_houghV_hit_x.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getPos().x());
          m_houghV_hit_y.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getPos().y());
          m_houghV_hit_z.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getPos().z());
          m_houghV_hit_E.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getEnergy());
        }
   
        houghV_index++;
      }
    }
    t_Hough->Fill();
    // Cone
    ClearCone();
    int coneU_index=0;
    int coneV_index=0;
    for(int i=0; i<m_halfclusterU.size(); i++){  // loop half cluster U
      std::vector<const Cyber::CaloHalfCluster*> m_mergedaxisU = m_halfclusterU[i]->getHalfClusterCol("ConeAxis");
      for(int ita=0; ita<m_mergedaxisU.size(); ita++){ // loop  axis U
        // General information of the axis
        m_coneU_tag.push_back(coneU_index);
        m_coneU_type.push_back(m_mergedaxisU[ita]->getType());
        m_coneU_x.push_back(m_mergedaxisU[ita]->getPos().x());
        m_coneU_y.push_back(m_mergedaxisU[ita]->getPos().y());
        m_coneU_z.push_back(m_mergedaxisU[ita]->getPos().z());
        m_coneU_E.push_back(m_mergedaxisU[ita]->getEnergy());
   
        // MC truth information of the Axis
        auto truthMap = m_mergedaxisU[ita]->getLinkedMCP();
        for(auto iter: truthMap){
          m_coneU_truth_tag.push_back(coneU_index);
          m_coneU_truth_MC_px.push_back(iter.first.getMomentum().x);
          m_coneU_truth_MC_py.push_back(iter.first.getMomentum().y);
          m_coneU_truth_MC_pz.push_back(iter.first.getMomentum().z);
          m_coneU_truth_MC_E.push_back(iter.first.getEnergy());
          m_coneU_truth_MC_weight.push_back(iter.second);
        }
   
        // Hits on axis
        for(int ilm=0; ilm<m_mergedaxisU[ita]->getCluster().size(); ilm++){ // loop local max
          m_coneU_hit_tag.push_back(coneU_index);
          m_coneU_hit_x.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getPos().x() );
          m_coneU_hit_y.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getPos().y() );
          m_coneU_hit_z.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getPos().z() );
          m_coneU_hit_E.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getEnergy()  );
        }
   
        coneU_index++;
      }
    }
    for(int i=0; i<m_halfclusterV.size(); i++){  // loop half cluster V
      std::vector<const Cyber::CaloHalfCluster*> m_mergedaxisV = m_halfclusterV[i]->getHalfClusterCol("ConeAxis");
      for(int ita=0; ita<m_mergedaxisV.size(); ita++){ // loop  axis V
        // General information of the axis
        m_coneV_tag.push_back(coneV_index);
        m_coneV_type.push_back(m_mergedaxisV[ita]->getType());
        m_coneV_x.push_back(m_mergedaxisV[ita]->getPos().x());
        m_coneV_y.push_back(m_mergedaxisV[ita]->getPos().y());
        m_coneV_z.push_back(m_mergedaxisV[ita]->getPos().z());
        m_coneV_E.push_back(m_mergedaxisV[ita]->getEnergy());
   
        // MC truth information of the Axis
        auto truthMap = m_mergedaxisV[ita]->getLinkedMCP();
        for(auto iter: truthMap){
          m_coneV_truth_tag.push_back(coneV_index);
          m_coneV_truth_MC_px.push_back(iter.first.getMomentum().x);
          m_coneV_truth_MC_py.push_back(iter.first.getMomentum().y);
          m_coneV_truth_MC_pz.push_back(iter.first.getMomentum().z);
          m_coneV_truth_MC_E.push_back(iter.first.getEnergy());
          m_coneV_truth_MC_weight.push_back(iter.second);
        }
   
        // Hits on axis
        for(int ilm=0; ilm<m_mergedaxisV[ita]->getCluster().size(); ilm++){
          m_coneV_hit_tag.push_back(coneV_index);
          m_coneV_hit_x.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getPos().x());
          m_coneV_hit_y.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getPos().y());
          m_coneV_hit_z.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getPos().z());
          m_coneV_hit_E.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getEnergy());
        }
   
        coneV_index++;
      }
    }
    t_Cone->Fill();
    // Track axis
    ClearTrackAxis();
    int trackU_index=0;
    int trackV_index=0;
    for(int i=0; i<m_halfclusterU.size(); i++){  // loop half cluster U
      std::vector<const Cyber::CaloHalfCluster*> m_mergedaxisU = m_halfclusterU[i]->getHalfClusterCol("TrackAxis");
      for(int ita=0; ita<m_mergedaxisU.size(); ita++){ // loop  axis U
        // General information of the axis
        m_trackU_tag.push_back(trackU_index);
        m_trackU_type.push_back(m_mergedaxisU[ita]->getType());
        m_trackU_x.push_back(m_mergedaxisU[ita]->getPos().x());
        m_trackU_y.push_back(m_mergedaxisU[ita]->getPos().y());
        m_trackU_z.push_back(m_mergedaxisU[ita]->getPos().z());
        m_trackU_E.push_back(m_mergedaxisU[ita]->getEnergy());
   
        // MC truth information of the Axis
        auto truthMap = m_mergedaxisU[ita]->getLinkedMCP();
        for(auto iter: truthMap){
          m_trackU_truth_tag.push_back(trackU_index);
          m_trackU_truth_MC_px.push_back(iter.first.getMomentum().x);
          m_trackU_truth_MC_py.push_back(iter.first.getMomentum().y);
          m_trackU_truth_MC_pz.push_back(iter.first.getMomentum().z);
          m_trackU_truth_MC_E.push_back(iter.first.getEnergy());
          m_trackU_truth_MC_weight.push_back(iter.second);
        }
   
        // Hits on axis
        for(int ilm=0; ilm<m_mergedaxisU[ita]->getCluster().size(); ilm++){ // loop local max
          m_trackU_hit_tag.push_back(trackU_index);
          m_trackU_hit_x.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getPos().x() );
          m_trackU_hit_y.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getPos().y() );
          m_trackU_hit_z.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getPos().z() );
          m_trackU_hit_E.push_back( m_mergedaxisU[ita]->getCluster()[ilm]->getEnergy()  );
        }
   
        trackU_index++;
      }
    }
    for(int i=0; i<m_halfclusterV.size(); i++){  // loop half cluster V
      std::vector<const Cyber::CaloHalfCluster*> m_mergedaxisV = m_halfclusterV[i]->getHalfClusterCol("TrackAxis");
      for(int ita=0; ita<m_mergedaxisV.size(); ita++){ // loop  axis V
        // General information of the axis
        m_trackV_tag.push_back(trackV_index);
        m_trackV_type.push_back(m_mergedaxisV[ita]->getType());
        m_trackV_x.push_back(m_mergedaxisV[ita]->getPos().x());
        m_trackV_y.push_back(m_mergedaxisV[ita]->getPos().y());
        m_trackV_z.push_back(m_mergedaxisV[ita]->getPos().z());
        m_trackV_E.push_back(m_mergedaxisV[ita]->getEnergy());
   
        // MC truth information of the Axis
        auto truthMap = m_mergedaxisV[ita]->getLinkedMCP();
        for(auto iter: truthMap){
          m_trackV_truth_tag.push_back(trackV_index);
          m_trackV_truth_MC_px.push_back(iter.first.getMomentum().x);
          m_trackV_truth_MC_py.push_back(iter.first.getMomentum().y);
          m_trackV_truth_MC_pz.push_back(iter.first.getMomentum().z);
          m_trackV_truth_MC_E.push_back(iter.first.getEnergy());
          m_trackV_truth_MC_weight.push_back(iter.second);
        }
   
        // Hits on axis
        for(int ilm=0; ilm<m_mergedaxisV[ita]->getCluster().size(); ilm++){
          m_trackV_hit_tag.push_back(trackV_index);
          m_trackV_hit_x.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getPos().x());
          m_trackV_hit_y.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getPos().y());
          m_trackV_hit_z.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getPos().z());
          m_trackV_hit_E.push_back(m_mergedaxisV[ita]->getCluster()[ilm]->getEnergy());
        }
   
        trackV_index++;
      }
    }
    t_TrackAxis->Fill();
    //Axis
    ClearAxis();
    int axisU_index=0;
    int axisV_index=0;
    for(int i=0; i<m_halfclusterU.size(); i++){  // loop half cluster U
      std::vector<const Cyber::CaloHalfCluster*> m_mergedaxisU = m_halfclusterU[i]->getHalfClusterCol("MergedAxis");
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
      std::vector<const Cyber::CaloHalfCluster*> m_mergedaxisV = m_halfclusterV[i]->getHalfClusterCol("MergedAxis");
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
    
    t_Axis->Fill();
   
   
    //Half cluster
    ClearHalfCluster();
    m_halfclusterV.clear();
    m_halfclusterU.clear();
    m_totE_HFClusV = 0;
    m_totE_HFClusU = 0;
    //for(int i=0; i<m_DataCol.map_HalfCluster["ESHalfClusterU"].size(); i++){
    //  m_halfclusterU.push_back( m_DataCol.map_HalfCluster["ESHalfClusterU"][i]->getHalfClusterCol("MergedAxis")[0] );
    //}
    //for(int i=0; i<m_DataCol.map_HalfCluster["ESHalfClusterV"].size(); i++){
    //  m_halfclusterV.push_back( m_DataCol.map_HalfCluster["ESHalfClusterV"][i]->getHalfClusterCol("MergedAxis")[0] );
    //}
    for(int i=0; i<m_DataCol.map_HalfCluster["HalfClusterColU"].size(); i++){
      m_halfclusterU.push_back( m_DataCol.map_HalfCluster["HalfClusterColU"][i].get() );
    }
    for(int i=0; i<m_DataCol.map_HalfCluster["HalfClusterColV"].size(); i++){
      m_halfclusterV.push_back( m_DataCol.map_HalfCluster["HalfClusterColV"][i].get() );
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
   
        // Bars (hits)
        for(int ilm=0; ilm<m_halfclusterV[i]->getBars().size(); ilm++){ 
          m_HalfClusterV_hit_tag.push_back(i);
          m_HalfClusterV_hit_x.push_back( m_halfclusterV[i]->getBars()[ilm]->getPosition().x() );
          m_HalfClusterV_hit_y.push_back( m_halfclusterV[i]->getBars()[ilm]->getPosition().y() );
          m_HalfClusterV_hit_z.push_back( m_halfclusterV[i]->getBars()[ilm]->getPosition().z() );
          m_HalfClusterV_hit_E.push_back( m_halfclusterV[i]->getBars()[ilm]->getEnergy()  );
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
   
        // Bars (hits)
        for(int ilm=0; ilm<m_halfclusterU[i]->getBars().size(); ilm++){
          m_HalfClusterU_hit_tag.push_back(i);
          m_HalfClusterU_hit_x.push_back( m_halfclusterU[i]->getBars()[ilm]->getPosition().x() );
          m_HalfClusterU_hit_y.push_back( m_halfclusterU[i]->getBars()[ilm]->getPosition().y() );
          m_HalfClusterU_hit_z.push_back( m_halfclusterU[i]->getBars()[ilm]->getPosition().z() );
          m_HalfClusterU_hit_E.push_back( m_halfclusterU[i]->getBars()[ilm]->getEnergy()  );
        }
    }
    t_HalfCluster->Fill();
   
   
    //Tower
    ClearTower();
    std::vector<std::shared_ptr<Cyber::Calo3DCluster>> m_tower = m_DataCol.map_CaloCluster["ESTower"];
    m_Ntower = m_tower.size();
    for(int it=0; it<m_tower.size(); it++){
      //ClearTower();
      m_towerID_id1.push_back(m_tower[it]->getTowerID()[0][1]);
      m_towerID_id2.push_back(m_tower[it]->getTowerID()[0][2]);

      std::vector<const CaloHalfCluster*> m_HFClusU = m_tower[it]->getHalfClusterUCol("ESHalfClusterU");
      std::vector<const CaloHalfCluster*> m_HFClusV = m_tower[it]->getHalfClusterVCol("ESHalfClusterV");
   
      m_NclusU.push_back(m_HFClusU.size());
      m_NclusV.push_back(m_HFClusV.size());
      m_totEn.push_back(m_tower[it]->getEnergy());
      float tmp_totEn_U = 0.;
      float tmp_totEn_V = 0.;
      for(int ic=0; ic<m_HFClusU.size(); ic++){
        m_HalfClusterU_tag.push_back(it);
        m_HalfClusterU_x.push_back(m_HFClusU[ic]->getPos().x());
        m_HalfClusterU_y.push_back(m_HFClusU[ic]->getPos().y());
        m_HalfClusterU_z.push_back(m_HFClusU[ic]->getPos().z());
        m_HalfClusterU_E.push_back(m_HFClusU[ic]->getEnergy());
        m_HalfClusterU_type.push_back(m_HFClusU[ic]->getType());
        m_HalfClusterU_nTrk.push_back(m_HFClusU[ic]->getAssociatedTracks().size());
        tmp_totEn_U += m_HFClusU[ic]->getEnergy();
      }
   
      for(int ic=0; ic<m_HFClusV.size(); ic++){
        m_HalfClusterV_tag.push_back(it);
        m_HalfClusterV_x.push_back(m_HFClusV[ic]->getPos().x());
        m_HalfClusterV_y.push_back(m_HFClusV[ic]->getPos().y());
        m_HalfClusterV_z.push_back(m_HFClusV[ic]->getPos().z());
        m_HalfClusterV_E.push_back(m_HFClusV[ic]->getEnergy());
        m_HalfClusterV_type.push_back(m_HFClusV[ic]->getType());
        m_HalfClusterV_nTrk.push_back(m_HFClusV[ic]->getAssociatedTracks().size());
        tmp_totEn_V += m_HFClusV[ic]->getEnergy();
      }
      m_totEn_U.push_back(tmp_totEn_U);
      m_totEn_V.push_back(tmp_totEn_V);
    }
    t_Tower->Fill();

    cout<<"  Write 3D cluster"<<endl;
    //3D cluster
    ClearCluster();
    std::vector<std::shared_ptr<Cyber::Calo3DCluster>> m_EcalClusterCol = m_DataCol.map_CaloCluster["TrkMergedECAL"];
    std::vector<std::shared_ptr<Cyber::Calo3DCluster>> m_HcalClusterCol = m_DataCol.map_CaloCluster["HCALCluster"];
    std::vector<std::shared_ptr<Cyber::Calo3DCluster>> m_SimpleHcalClusterCol = m_DataCol.map_CaloCluster["SimpleHCALCluster"];
    m_totE_Ecal = 0.;
    m_totE_Hcal = 0.;
    m_Nclus_Ecal = m_EcalClusterCol.size();
    m_Nclus_Hcal = m_SimpleHcalClusterCol.size();
    for(int icl=0; icl<m_EcalClusterCol.size(); icl++){
      m_EcalClus_x.push_back(m_EcalClusterCol[icl]->getShowerCenter().x());
      m_EcalClus_y.push_back(m_EcalClusterCol[icl]->getShowerCenter().y());
      m_EcalClus_z.push_back(m_EcalClusterCol[icl]->getShowerCenter().z());
      m_EcalClus_E.push_back(m_EcalClusterCol[icl]->getLongiE());
      m_EcalClus_nTrk.push_back(m_EcalClusterCol[icl]->getAssociatedTracks().size());
   
      //double tmp_phi = std::atan2(m_EcalClusterCol[icl]->getShowerCenter().y(), m_EcalClusterCol[icl]->getShowerCenter().x())* 180.0 / M_PI;
      //if (tmp_phi < 0) tmp_phi += 360.0;
      //double tmp_theta = std::atan2(m_EcalClusterCol[icl]->getShowerCenter().z(), m_EcalClusterCol[icl]->getShowerCenter().Perp())* 180.0 / M_PI + 90; 
      //cout<<" Theta: "<<tmp_theta<<" Phi: "<<tmp_phi<<endl;
      //m_EcalClus_Escale.push_back(m_energycorsvc->energyCorrection(m_EcalClusterCol[icl]->getLongiE(), tmp_phi, tmp_theta));
      m_EcalClus_Escale.push_back(m_EcalClusterCol[icl]->getLongiE());
   
   
      if(m_EcalClusterCol[icl]->getAssociatedTracks().size()==1){
        const Track* trk = m_EcalClusterCol[icl]->getAssociatedTracks()[0];
        m_EcalClus_pTrk.push_back(trk->getMomentum());
   
        std::vector<TrackState> AllTrackStates = trk->getAllTrackStates();
        for(int istate=0; istate<AllTrackStates.size(); istate++){
          m_EcalClus_trk_tag.push_back(icl);
          m_EcalClus_trk_d0.push_back(AllTrackStates[istate].D0);
          m_EcalClus_trk_z0.push_back(AllTrackStates[istate].Z0);
          m_EcalClus_trk_phi.push_back(AllTrackStates[istate].phi0);
          m_EcalClus_trk_tanL.push_back( AllTrackStates[istate].tanLambda );
          m_EcalClus_trk_kappa.push_back( AllTrackStates[istate].Kappa);
          m_EcalClus_trk_omega.push_back( AllTrackStates[istate].Omega );
          m_EcalClus_trk_location.push_back( AllTrackStates[istate].location );
        }
   
      }
      else
        m_EcalClus_pTrk.push_back(-99);
   
      m_EcalClus_typeU.push_back(m_EcalClusterCol[icl]->getHalfClusterUCol("LinkedLongiCluster")[0]->getType());
      m_EcalClus_typeV.push_back(m_EcalClusterCol[icl]->getHalfClusterVCol("LinkedLongiCluster")[0]->getType());
      for(int ii=0; ii<m_EcalClusterCol[icl]->getHalfClusterUCol("LinkedLongiCluster").size(); ii++){
        for(int ihit=0; ihit<m_EcalClusterCol[icl]->getHalfClusterUCol("LinkedLongiCluster")[ii]->getBars().size(); ihit++){
          auto shower = m_EcalClusterCol[icl]->getHalfClusterUCol("LinkedLongiCluster")[ii]->getBars()[ihit];
          m_EcalClus_hitU_tag.push_back(icl);
          m_EcalClus_hitU_x.push_back(shower->getPosition().x());
          m_EcalClus_hitU_y.push_back(shower->getPosition().y());
          m_EcalClus_hitU_z.push_back(shower->getPosition().z());
          m_EcalClus_hitU_E.push_back(shower->getEnergy());
        }
      }
      for(int ii=0; ii<m_EcalClusterCol[icl]->getHalfClusterVCol("LinkedLongiCluster").size(); ii++){
        for(int ihit=0; ihit<m_EcalClusterCol[icl]->getHalfClusterVCol("LinkedLongiCluster")[ii]->getBars().size(); ihit++){
          auto shower = m_EcalClusterCol[icl]->getHalfClusterVCol("LinkedLongiCluster")[ii]->getBars()[ihit];
          m_EcalClus_hitV_tag.push_back(icl);
          m_EcalClus_hitV_x.push_back(shower->getPosition().x());
          m_EcalClus_hitV_y.push_back(shower->getPosition().y());
          m_EcalClus_hitV_z.push_back(shower->getPosition().z());
          m_EcalClus_hitV_E.push_back(shower->getEnergy());
        }
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
      m_HcalClus_nHit.push_back(m_HcalClusterCol[icl]->getCaloHits().size());
   
      for(int ih=0; ih<m_HcalClusterCol[icl]->getCaloHits().size(); ih++){
        m_HcalClus_hit_tag.push_back(icl);
        m_HcalClus_hit_x.push_back(m_HcalClusterCol[icl]->getCaloHits()[ih]->getPosition().x());
        m_HcalClus_hit_y.push_back(m_HcalClusterCol[icl]->getCaloHits()[ih]->getPosition().y());
        m_HcalClus_hit_z.push_back(m_HcalClusterCol[icl]->getCaloHits()[ih]->getPosition().z());
        m_HcalClus_hit_E.push_back(m_HcalClusterCol[icl]->getCaloHits()[ih]->getEnergy());
      }
   
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
      m_SimpleHcalClus_nHit.push_back(m_SimpleHcalClusterCol[icl]->getCaloHits().size());
   
      for(int ih=0; ih<m_SimpleHcalClusterCol[icl]->getCaloHits().size(); ih++){
        m_SimpleHcalClus_hit_tag.push_back(icl);
        m_SimpleHcalClus_hit_x.push_back(m_SimpleHcalClusterCol[icl]->getCaloHits()[ih]->getPosition().x());
        m_SimpleHcalClus_hit_y.push_back(m_SimpleHcalClusterCol[icl]->getCaloHits()[ih]->getPosition().y());
        m_SimpleHcalClus_hit_z.push_back(m_SimpleHcalClusterCol[icl]->getCaloHits()[ih]->getPosition().z());
        m_SimpleHcalClus_hit_E.push_back(m_SimpleHcalClusterCol[icl]->getCaloHits()[ih]->getEnergy());
      }
   
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
   
    // Save Track info
    ClearTrack();
    std::vector<Cyber::Track*> m_trkCol; 
    for(int it=0; it<m_DataCol.TrackCol.size(); it++)
      m_trkCol.push_back( m_DataCol.TrackCol[it].get() );
   
    m_Ntrk = m_trkCol.size();
    for(int itrk=0; itrk<m_Ntrk; itrk++){
      m_type.push_back(m_trkCol[itrk]->getType());
      m_Nhit.push_back(m_trkCol[itrk]->getTrackerHits());
      m_pid.push_back(m_trkCol[itrk]->getPID());
      m_pid_truth.push_back(m_trkCol[itrk]->getLeadingMCP().getPDG());
      m_trk_truthweight.push_back( m_trkCol[itrk]->getLeadingMCPweight() );
      m_trk_px.push_back( m_trkCol[itrk]->getP3().x() );
      m_trk_py.push_back( m_trkCol[itrk]->getP3().y() );
      m_trk_pz.push_back( m_trkCol[itrk]->getP3().z() );
      m_trk_p.push_back( m_trkCol[itrk]->getMomentum() );
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
      }
      std::vector<TrackState> EcalTrackStates = m_trkCol[itrk]->getTrackStates("Ecal");
      for(int istate=0; istate<EcalTrackStates.size(); istate++){
        m_trkstate_x_ECAL.push_back(EcalTrackStates[istate].referencePoint.X());
        m_trkstate_y_ECAL.push_back(EcalTrackStates[istate].referencePoint.Y());
        m_trkstate_z_ECAL.push_back(EcalTrackStates[istate].referencePoint.Z());
        m_trkstate_tag_ECAL.push_back(itrk);
      }
      std::vector<TrackState> HcalTrackStates = m_trkCol[itrk]->getTrackStates("Hcal");
      for(int istate=0; istate<HcalTrackStates.size(); istate++){
        m_trkstate_x_ECAL.push_back(HcalTrackStates[istate].referencePoint.X());
        m_trkstate_y_ECAL.push_back(HcalTrackStates[istate].referencePoint.Y());
        m_trkstate_z_ECAL.push_back(HcalTrackStates[istate].referencePoint.Z());
        m_trkstate_tag_HCAL.push_back(itrk);
      }
    }
    t_Track->Fill();
   
    // yyy: pfo
    ClearPFO();
    std::vector<Cyber::PFObject*> m_pfobjects; m_pfobjects.clear();
    for(int ip=0; ip<m_DataCol.map_PFObjects["outputPFO"].size(); ip++)
      m_pfobjects.push_back(m_DataCol.map_PFObjects["outputPFO"][ip].get());
   
   
    for(int ip=0; ip<m_pfobjects.size(); ip++){
      std::vector<const Track*> t_tracks = m_pfobjects[ip]->getTracks();
      std::vector<const Calo3DCluster*> t_ecal_clusters = m_pfobjects[ip]->getECALClusters();
      std::vector<const Calo3DCluster*> t_hcal_clusters =  m_pfobjects[ip]->getHCALClusters();
   
      pfo_tag.push_back(ip);
      pfo_n_track.push_back(t_tracks.size());
      pfo_n_ecal_clus.push_back(t_ecal_clusters.size());
      pfo_n_hcal_clus.push_back(t_hcal_clusters.size());
   
      for(int it=0; it<t_tracks.size(); it++){
        std::vector<TrackState> AllTrackStates = t_tracks[it]->getAllTrackStates();
        for(int istate=0; istate<AllTrackStates.size(); istate++){
          pfo_trk_tag.push_back(ip);
          pfo_trk_d0.push_back( AllTrackStates[istate].D0 );
          pfo_trk_z0.push_back( AllTrackStates[istate].Z0 );
          pfo_trk_phi.push_back( AllTrackStates[istate].phi0 );
          pfo_trk_tanL.push_back( AllTrackStates[istate].tanLambda );
          pfo_trk_kappa.push_back( AllTrackStates[istate].Kappa);
          pfo_trk_omega.push_back( AllTrackStates[istate].Omega );
          pfo_trk_location.push_back( AllTrackStates[istate].location );
        }
      }
   
      for(int ie=0; ie<t_ecal_clusters.size(); ie++){
        pfo_ecal_tag.push_back(ip);
        pfo_ecal_clus_x.push_back(t_ecal_clusters[ie]->getShowerCenter().x());
        pfo_ecal_clus_y.push_back(t_ecal_clusters[ie]->getShowerCenter().y());
        pfo_ecal_clus_z.push_back(t_ecal_clusters[ie]->getShowerCenter().z());
        pfo_ecal_clus_E.push_back(t_ecal_clusters[ie]->getLongiE());
   
        //double tmp_phi = std::atan2(t_ecal_clusters[ie]->getShowerCenter().y(), t_ecal_clusters[ie]->getShowerCenter().x())* 180.0 / M_PI;
        //if (tmp_phi < 0) tmp_phi += 360.0;
        //double tmp_theta = std::atan2(t_ecal_clusters[ie]->getShowerCenter().z(), t_ecal_clusters[ie]->getShowerCenter().Perp())* 180.0 / M_PI + 90; 
        pfo_ecal_clus_Escale.push_back(t_ecal_clusters[ie]->getLongiE());
   
      }
      for(int ih=0; ih<t_hcal_clusters.size(); ih++){
        pfo_hcal_tag.push_back(ip);
        pfo_hcal_clus_x.push_back(t_hcal_clusters[ih]->getHitCenter().x());
        pfo_hcal_clus_y.push_back(t_hcal_clusters[ih]->getHitCenter().y());
        pfo_hcal_clus_z.push_back(t_hcal_clusters[ih]->getHitCenter().z());
        pfo_hcal_clus_E.push_back(t_hcal_clusters[ih]->getHitsE());
      }
    }
    t_PFO->Fill();

  }

  //Clean Events
  //system("/cefs/higgs/songwz/winter22/CEPCSW/workarea/memory/memory_test.sh before_clean");
  m_DataCol.Clear();
  //system("/cefs/higgs/songwz/winter22/CEPCSW/workarea/memory/memory_test.sh event_end");


  std::cout<<"Event: "<<_nEvt<<" is done"<<std::endl;
  _nEvt ++ ;
  return StatusCode::SUCCESS;
}

StatusCode CyberPFAlg::finalize()
{
  if(m_WriteAna){
    m_wfile->cd();
    t_MCParticle->Write();
    t_SimBar->Write();
    t_LocalMax->Write();
    t_Layers->Write();
    t_Hough->Write();
    t_Cone->Write();
    t_TrackAxis->Write();
    t_Axis->Write();
    t_HalfCluster->Write();
    t_Tower->Write();
    t_Cluster->Write();
    t_Track->Write();
    t_PFO->Write();
    m_wfile->Close();
    delete m_wfile, t_MCParticle, t_SimBar, t_LocalMax, t_Layers, t_Hough, t_Cone, t_TrackAxis, t_Axis, t_HalfCluster, t_Tower, t_Cluster, t_Track, t_PFO; 
  }

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
  for(auto iter : w_ClusterCollection) delete iter.second;
  w_ClusterCollection.clear(); 
  //r_ECalHitCols.clear();
  //r_HCalHitCols.clear();
  r_CaloHitCols.clear();
  //m_energycorsvc->finalize();
  //delete m_cellIDConverter, m_geosvc;
  info() << "Processed " << _nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}

void CyberPFAlg::ClearMCParticle(){
  m_mcPdgid.clear();
  m_mcStatus.clear();
  m_mcPx.clear();
  m_mcPy.clear();
  m_mcPz.clear();
  m_mcEn.clear();
  m_mcMass.clear();
  m_mcCharge.clear();
  m_mcVTXx.clear();
  m_mcVTXy.clear();
  m_mcVTXz.clear();
  m_mcEPx.clear();
  m_mcEPy.clear();
  m_mcEPz.clear();
  m_depEn_ecal.clear();
  m_depEn_hcal.clear();
}

void CyberPFAlg::ClearBar(){
  m_totE_EcalSim = -99;
  m_totE_HcalSim = -99;
  m_simBar_x.clear();
  m_simBar_y.clear();
  m_simBar_z.clear();
  m_simBar_length.clear();
  m_simBar_nBarInLayer.clear();
  m_simBar_T1.clear();
  m_simBar_T2.clear();
  m_simBar_Q1.clear();
  m_simBar_Q2.clear();
  m_simBar_module.clear();
  m_simBar_dlayer.clear();
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

void CyberPFAlg::ClearLocalMax(){
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

void CyberPFAlg::ClearLayer(){
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


void CyberPFAlg::ClearHough(){
  m_houghU_tag.clear();
  m_houghU_type.clear();
  m_houghU_x.clear();
  m_houghU_y.clear();
  m_houghU_z.clear();
  m_houghU_E.clear();
  m_houghU_truth_tag.clear();
  m_houghU_truth_MC_px.clear();
  m_houghU_truth_MC_py.clear();
  m_houghU_truth_MC_pz.clear();
  m_houghU_truth_MC_E.clear();
  m_houghU_truth_MC_weight.clear();
  m_houghU_hit_tag.clear();
  m_houghU_hit_x.clear();
  m_houghU_hit_y.clear();
  m_houghU_hit_z.clear();
  m_houghU_hit_E.clear();
  m_houghV_tag.clear();
  m_houghV_type.clear();
  m_houghV_x.clear();
  m_houghV_y.clear();
  m_houghV_z.clear();
  m_houghV_E.clear();
  m_houghV_alpha.clear();
  m_houghV_rho.clear();
  m_houghV_truth_tag.clear();
  m_houghV_truth_MC_px.clear();
  m_houghV_truth_MC_py.clear();
  m_houghV_truth_MC_pz.clear();
  m_houghV_truth_MC_E.clear();
  m_houghV_truth_MC_weight.clear();
  m_houghV_hit_tag.clear();
  m_houghV_hit_x.clear();
  m_houghV_hit_y.clear();
  m_houghV_hit_z.clear();
  m_houghV_hit_E.clear();
}


void CyberPFAlg::ClearCone(){
  m_coneU_tag.clear();
  m_coneU_type.clear();
  m_coneU_x.clear();
  m_coneU_y.clear();
  m_coneU_z.clear();
  m_coneU_E.clear();
  m_coneU_truth_tag.clear();
  m_coneU_truth_MC_px.clear();
  m_coneU_truth_MC_py.clear();
  m_coneU_truth_MC_pz.clear();
  m_coneU_truth_MC_E.clear();
  m_coneU_truth_MC_weight.clear();
  m_coneU_hit_tag.clear();
  m_coneU_hit_x.clear();
  m_coneU_hit_y.clear();
  m_coneU_hit_z.clear();
  m_coneU_hit_E.clear();
  m_coneV_tag.clear();
  m_coneV_type.clear();
  m_coneV_x.clear();
  m_coneV_y.clear();
  m_coneV_z.clear();
  m_coneV_E.clear();
  m_coneV_truth_tag.clear();
  m_coneV_truth_MC_px.clear();
  m_coneV_truth_MC_py.clear();
  m_coneV_truth_MC_pz.clear();
  m_coneV_truth_MC_E.clear();
  m_coneV_truth_MC_weight.clear();
  m_coneV_hit_tag.clear();
  m_coneV_hit_x.clear();
  m_coneV_hit_y.clear();
  m_coneV_hit_z.clear();
  m_coneV_hit_E.clear();
}


void CyberPFAlg::ClearTrackAxis(){
  m_trackU_tag.clear();
  m_trackU_type.clear();
  m_trackU_x.clear();
  m_trackU_y.clear();
  m_trackU_z.clear();
  m_trackU_E.clear();
  m_trackU_truth_tag.clear();
  m_trackU_truth_MC_px.clear();
  m_trackU_truth_MC_py.clear();
  m_trackU_truth_MC_pz.clear();
  m_trackU_truth_MC_E.clear();
  m_trackU_truth_MC_weight.clear();
  m_trackU_hit_tag.clear();
  m_trackU_hit_x.clear();
  m_trackU_hit_y.clear();
  m_trackU_hit_z.clear();
  m_trackU_hit_E.clear();
  m_trackV_tag.clear();
  m_trackV_type.clear();
  m_trackV_x.clear();
  m_trackV_y.clear();
  m_trackV_z.clear();
  m_trackV_E.clear();
  m_trackV_truth_tag.clear();
  m_trackV_truth_MC_px.clear();
  m_trackV_truth_MC_py.clear();
  m_trackV_truth_MC_pz.clear();
  m_trackV_truth_MC_E.clear();
  m_trackV_truth_MC_weight.clear();
  m_trackV_hit_tag.clear();
  m_trackV_hit_x.clear();
  m_trackV_hit_y.clear();
  m_trackV_hit_z.clear();
  m_trackV_hit_E.clear();
}


void CyberPFAlg::ClearAxis(){
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

void CyberPFAlg::ClearHalfCluster(){
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

void CyberPFAlg::ClearTower(){
  m_Ntower = 0;
  m_towerID_id1.clear();
  m_towerID_id2.clear();
  m_NclusU.clear();
  m_NclusV.clear();
  m_totEn.clear();
  m_totEn_U.clear();
  m_totEn_V.clear();
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

void CyberPFAlg::ClearCluster(){
  m_totE_Ecal = -99.;
  m_totE_Hcal = -99.;
  m_Nclus_Ecal = -99;
  m_Nclus_Hcal = -99;
  m_EcalClus_x.clear();
  m_EcalClus_y.clear();
  m_EcalClus_z.clear();
  m_EcalClus_E.clear();
  m_EcalClus_Escale.clear();
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
  m_EcalClus_trk_location.clear();
  m_EcalClus_trk_tag.clear();
  m_EcalClus_trk_d0.clear();
  m_EcalClus_trk_z0.clear();
  m_EcalClus_trk_phi.clear();
  m_EcalClus_trk_tanL.clear();
  m_EcalClus_trk_kappa.clear();
  m_EcalClus_trk_omega.clear();
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
  m_HcalClus_nHit.clear();
  m_HcalClus_nTrk.clear();
  m_HcalClus_pTrk.clear();
  m_HcalClus_hit_x.clear();
  m_HcalClus_hit_y.clear();
  m_HcalClus_hit_z.clear();
  m_HcalClus_hit_E.clear();
  m_HcalClus_hit_tag.clear();
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
  m_SimpleHcalClus_nHit.clear();
  m_SimpleHcalClus_nTrk.clear();
  m_SimpleHcalClus_pTrk.clear();
  m_SimpleHcalClus_hit_x.clear();
  m_SimpleHcalClus_hit_y.clear();
  m_SimpleHcalClus_hit_z.clear();
  m_SimpleHcalClus_hit_E.clear();
  m_SimpleHcalClus_hit_tag.clear();
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

void CyberPFAlg::ClearTrack(){
  m_type.clear();
  m_Nhit.clear();
  m_pid.clear();
  m_pid_truth.clear();
  m_trk_px.clear();
  m_trk_py.clear();
  m_trk_pz.clear();
  m_trk_p.clear();
  m_trk_truthweight.clear();
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
  m_trkstate_x_ECAL.clear();
  m_trkstate_y_ECAL.clear();
  m_trkstate_z_ECAL.clear();
  m_trkstate_x_HCAL.clear();
  m_trkstate_y_HCAL.clear();
  m_trkstate_z_HCAL.clear();
  m_trkstate_tag_ECAL.clear();
  m_trkstate_tag_HCAL.clear();

}

void CyberPFAlg::ClearPFO(){
  pfo_tag.clear();
  pfo_n_track.clear();
  pfo_n_ecal_clus.clear();
  pfo_n_hcal_clus.clear();
  pfo_trk_tag.clear();
  pfo_trk_d0.clear();
  pfo_trk_z0.clear();
  pfo_trk_phi.clear();
  pfo_trk_tanL.clear();
  pfo_trk_kappa.clear();
  pfo_trk_omega.clear();
  pfo_trk_location.clear();
  pfo_ecal_tag.clear();
  pfo_ecal_clus_x.clear();
  pfo_ecal_clus_y.clear();
  pfo_ecal_clus_z.clear();
  pfo_ecal_clus_E.clear();
  pfo_ecal_clus_Escale.clear();
  pfo_hcal_tag.clear();
  pfo_hcal_clus_x.clear();
  pfo_hcal_clus_y.clear();
  pfo_hcal_clus_z.clear();
  pfo_hcal_clus_E.clear();
}


double CyberPFAlg::GetParticleDepEnergy(edm4hep::MCParticle& mcp, std::vector<std::shared_ptr<Cyber::CaloUnit>>& barcol){

  double EnDep = 0.;
  for(int i=0; i<barcol.size(); i++){
    std::vector< std::pair<edm4hep::MCParticle, float> > mcp_map = barcol[i]->getLinkedMCP();
    for(int ip=0; ip<mcp_map.size(); ip++){
      if(mcp_map[ip].first==mcp){
        EnDep += barcol[i]->getEnergy() * mcp_map[ip].second;
      }
    }
  }
  return EnDep;

}

double CyberPFAlg::GetParticleDepEnergy(edm4hep::MCParticle& mcp, std::vector<std::shared_ptr<Cyber::CaloHit>>& hitcol){

  double EnDep = 0.;
  for(int i=0; i<hitcol.size(); i++){
    std::vector< std::pair<edm4hep::MCParticle, float> > mcp_map = hitcol[i]->getLinkedMCP();
    for(int ip=0; ip<mcp_map.size(); ip++){
      if(mcp_map[ip].first==mcp){
        EnDep += hitcol[i]->getEnergy() * mcp_map[ip].second;
      }
    }
  }
  return EnDep;

}

#endif
