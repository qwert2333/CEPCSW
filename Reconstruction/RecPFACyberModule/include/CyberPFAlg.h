//=============================================================
// CyberPFA: a PFA developed for CEPC referenece detector
// Ver. CyberPFA-5.0.1(2025.01.09)
//-------------------------------------------------------------
//  Author: Fangyi Guo, Yang Zhang, Weizheng Song, Shengsen Sun
//          (IHEP, CAS)
//  Contact: guofangyi@ihep.ac.cn,
//           sunss@ihep.ac.cn
//=============================================================
#ifndef PANDORAPLUS_ALG_H
#define PANDORAPLUS_ALG_H

#include <string>
#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/Segmentations.h>
#include "DetInterface/IGeomSvc.h"
#include "DetIdentifier/CEPCDetectorData.h"
#include <CrystalEcalSvc/ICrystalEcalSvc.h>

#include "k4FWCore/PodioDataSvc.h"
#include "podio/CollectionBase.h"
#include "podio/ROOTFrameWriter.h"

#include "CyberDataCol.h"
#include "Tools/MCParticleCreator.h"
#include "Tools/TrackCreator.h"
#include "Tools/CaloHitsCreator.h"
#include "Tools/OutputCreator.h"
#include "Tools/AlgorithmManager.h"
#include "Algorithm/ExampleAlg.h"
#include "Algorithm/GlobalClusteringAlg.h"
#include "Algorithm/HcalClusteringAlg.h"
#include "Algorithm/LocalMaxFindingAlg.h"
#include "Algorithm/TrackMatchingAlg.h"
#include "Algorithm/HoughClusteringAlg.h"
#include "Algorithm/ConeClustering2DAlg.h"
#include "Algorithm/AxisMergingAlg.h"
#include "Algorithm/EnergySplittingAlg.h"
#include "Algorithm/EnergyTimeMatchingAlg.h"
#include "Algorithm/ConeClusteringAlg.h"
#include "Algorithm/TrackExtrapolatingAlg.h"
#include "Algorithm/PFOCreatingAlg.h"
#include "Algorithm/TrackClusterConnectingAlg.h"
#include "Algorithm/PFOReclusteringAlg.h"

#include "Algorithm/TruthClusteringAlg.h"
#include "Algorithm/TruthTrackMatchingAlg.h"
#include "Algorithm/TruthPatternRecAlg.h"
#include "Algorithm/TruthEnergySplittingAlg.h"
#include "Algorithm/TruthMatchingAlg.h"
#include "Algorithm/TruthClusterMergingAlg.h"

#include "TVector3.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include <cstdlib>

using namespace Cyber;
using namespace std;

class CyberPFAlg : public GaudiAlgorithm
{
 
public:
 
  CyberPFAlg(const std::string& name, ISvcLocator* svcLoc);
 
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual StatusCode initialize() ;
 
  /** Called for every event - the working horse.
   */
  virtual StatusCode execute() ; 
 
  /** Called after data processing for clean up.
   */
  virtual StatusCode finalize() ;

protected:

  int _nEvt ;
  TRandom3 rndm;

  SmartIF<IGeomSvc> m_geosvc;
  SmartIF<ICrystalEcalSvc> m_energycorsvc;
  std::map<std::string, dd4hep::DDSegmentation::BitFieldCoder*> map_readout_decoder;
  dd4hep::Detector* m_dd4hep;
  //dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
  //dd4hep::VolumeManager m_volumeManager;
  std::map<std::tuple<int, int, int, int, int>, int> barNumberMapEndcapMap;
  //DataCollection: moved into execute() to ensure everything can be cleand after one event. 
  //CyberDataCol     m_DataCol; 


  //Creators and their setting
  MCParticleCreator       *m_pMCParticleCreator;
  TrackCreator            *m_pTrackCreator;
  CaloHitsCreator         *m_pCaloHitsCreator;
  OutputCreator           *m_pOutputCreator;

  Settings   m_pMCParticleCreatorSettings;
  Settings   m_pTrackCreatorSettings;
  Settings   m_CaloHitsCreatorSettings;
  Settings   m_OutputCreatorSettings;


  //Parameters for PFA algorithm
  Settings m_GlobalSettings; 


  //Algorithm for PFA
  Cyber::AlgorithmManager m_algorithmManager; 


  //Readin collection names
  Gaudi::Property< std::string > name_MCParticleCol{ this, "MCParticleCollection", "MCParticle" };
  Gaudi::Property< std::string > name_MCPTrkAssoCol{this, "MCRecoTrackParticleAssociationCollection", "MarlinTrkAssociation"};
  Gaudi::Property< std::vector<std::string> > name_TrackCol{ this, "TrackCollections", {"MarlinTrkTracks"} };
  Gaudi::Property< std::string > name_dNdxCol{this, "DndxCollection", "DndxTracks"};
  Gaudi::Property< std::string > name_tofCol{this, "TOFCollection", "RecTofCollection"};
  Gaudi::Property< std::vector<std::string> > name_EcalHits{ this, "ECalCaloHitCollections", {"ECALBarrel"} };
  Gaudi::Property< std::vector<std::string> > name_EcalReadout{ this, "ECalReadOutNames", {"EcalBarrelCollection"} }; 
  Gaudi::Property< std::vector<std::string> > name_EcalMCPAssociation{ this, "ECalMCPAssociationName", {"ECALBarrelParticleAssoCol"} };
  Gaudi::Property< std::vector<std::string> > name_HcalHits{ this, "HCalCaloHitCollections", {"HCALBarrel"} };
  Gaudi::Property< std::vector<std::string> > name_HcalReadout{ this, "HCalReadOutNames", {"HcalBarrelCollection"} }; 
  Gaudi::Property< std::vector<std::string> > name_HcalMCPAssociation{ this, "HCalMCPAssociationName", {"HCALBarrelParticleAssoCol"} };
  

  //---Readin collections
  typedef DataHandle<edm4hep::TrackCollection>                          TrackType; 
  typedef DataHandle<edm4hep::CalorimeterHitCollection>                 CaloType; 
  typedef DataHandle<edm4hep::MCRecoCaloParticleAssociationCollection>  CaloParticleAssoType; 
  DataHandle<edm4hep::MCParticleCollection>* r_MCParticleCol; 
  DataHandle<edm4hep::MCRecoTrackParticleAssociationCollection>* r_MCPTrkAssoCol;  
  std::vector<TrackType*> r_TrackCols; 
  //std::vector<CaloType*>  r_ECalHitCols; 
  //std::vector<CaloType*>  r_HCalHitCols; 
  std::vector<CaloType*>  r_CaloHitCols; 
  std::map<std::string, CaloParticleAssoType*> map_CaloMCPAssoCols;
  DataHandle<edm4hep::RecTofCollection>* r_TofCol;
  DataHandle<edm4hep::RecDqdxCollection>* r_dNdxCol;

  //Global parameters.
  Gaudi::Property<float>  m_BField{this,  "BField", 3., "Magnetic field"};
  Gaudi::Property<float>  m_seed{this,    "Seed", 2131, "Random Seed"};
  Gaudi::Property<int>    m_Debug{this,   "Debug", 0, "Debug level"};
  Gaudi::Property<int>    m_Nskip{this,   "SkipEvt", 0, "Skip event"};
  Gaudi::Property<std::string>   m_EcalType{this, "EcalType", "BarEcal", "ECAL type"};
  Gaudi::Property<bool>   m_useMCPTrk{this,  "UseMCPTrack", 1, "Use track from MCParticle extrapolation"};
  Gaudi::Property<bool>   m_useTruthMatchTrk{this,  "UseTruthMatchTrack", 1, "Use track from MCParticle extrapolation"};
  Gaudi::Property<bool>   m_doCleanTrack{this,  "DoCleanTrack", 1, "Do clean tracks"};
  Gaudi::Property<std::string>   m_trackIDFile{this,  "TrackIDFile", "", "BDT weight file for track ID"};
  Gaudi::Property<std::string>   m_trackIDMethod{this,  "TrackIDMethod", "BDTG", "BDT weight file for track ID"};
  Gaudi::Property<float>  m_EcalChargedCalib{this,  "EcalChargedCalib", 1.26, "ECAL global calibration"};
  Gaudi::Property<float>  m_HcalChargedCalib{this,  "HcalChargedCalib", 4.,  "HCAL global calibration"};
  Gaudi::Property<float>  m_EcalNeutralCalib{this,  "EcalNeutralCalib", 1., "ECAL global calibration"};
  Gaudi::Property<float>  m_HcalNeutralCalib{this,  "HcalNeutralCalib", 4.,  "HCAL global calibration"};
  
  //Algorithms: 
  typedef std::vector<std::string> StringVector;
  Gaudi::Property< StringVector > name_Algs{ this, "AlgList", {} };
  Gaudi::Property< std::vector<StringVector> > name_AlgPars{ this, "AlgParNames", {} };
  Gaudi::Property< std::vector<StringVector> > type_AlgPars{ this, "AlgParTypes", {} };
  Gaudi::Property< std::vector<StringVector> > value_AlgPars{this, "AlgParValues", {} };


  // Output collections
  DataHandle<edm4hep::CalorimeterHitCollection>             w_RecEcalCol{"RecECALBarrel", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::CalorimeterHitCollection>             w_RecCoreCol{"RecECALBarrelCore", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::CalorimeterHitCollection>             w_RecHcalCol{"RecHCALBarrel", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::TrackCollection>                      w_RecTrkCol{"RecTracks", Gaudi::DataHandle::Writer, this};
  //Gaudi::Property< std::string > name_EcalCluster{this, "OutputEcalCluster", "TrkMergedECAL"};
  //Gaudi::Property< std::string > name_EcalCore   {this, "OutputEcalCore", "EcalCore"};
  //Gaudi::Property< std::string > name_HcalCluster{this, "OutputHcalCluster", "HcalCluster"};
  Gaudi::Property< std::string > name_PFObject   {this, "OutputPFO", "outputPFO"};
  DataHandle<edm4hep::ReconstructedParticleCollection>      w_ReconstructedParticleCollection {"PandoraPFOs"    ,Gaudi::DataHandle::Writer, this};
  std::map<std::string, DataHandle<edm4hep::ClusterCollection>* > w_ClusterCollection;


  //For Ana
  Gaudi::Property<bool>  m_WriteAna {this, "WriteAna", false, "Write Ntuples for analysis"};
  Gaudi::Property<std::string> m_filename{this, "AnaFileName", "testout.root", "Output file name"};

  typedef std::vector<float> FloatVec;
  typedef std::vector<int>   IntVec;

  TFile* m_wfile;

  // MC particle
  TTree *t_MCParticle;
  int m_Nmc;
  IntVec m_mcPdgid, m_mcStatus;
  FloatVec m_mcPx, m_mcPy, m_mcPz, m_mcEn, m_mcMass, m_mcCharge;
  FloatVec m_mcVTXx, m_mcVTXy, m_mcVTXz, m_mcEPx, m_mcEPy, m_mcEPz, m_depEn_ecal, m_depEn_hcal;

  //Raw bars and hits
  TTree* t_SimBar;
  float m_totE_EcalSim, m_totE_HcalSim;
  FloatVec m_simBar_x, m_simBar_y, m_simBar_z, m_simBar_length, m_simBar_nBarInLayer, m_simBar_T1, m_simBar_T2, m_simBar_Q1, m_simBar_Q2; 
  FloatVec m_simBar_truthMC_tag, m_simBar_truthMC_pid, m_simBar_truthMC_px, m_simBar_truthMC_py, m_simBar_truthMC_pz, m_simBar_truthMC_E, 
           m_simBar_truthMC_EPx, m_simBar_truthMC_EPy, m_simBar_truthMC_EPz, m_simBar_truthMC_weight;   
  IntVec m_simBar_dlayer, m_simBar_stave, m_simBar_slayer, m_simBar_module, m_simBar_bar, m_simBar_system;
  FloatVec m_HcalHit_x, m_HcalHit_y, m_HcalHit_z, m_HcalHit_E, 
           m_HcalHit_truthMC_tag, m_HcalHit_truthMC_pid, m_HcalHit_truthMC_px, m_HcalHit_truthMC_py, m_HcalHit_truthMC_pz, m_HcalHit_truthMC_E,
           m_HcalHit_truthMC_EPx, m_HcalHit_truthMC_EPy, m_HcalHit_truthMC_EPz, m_HcalHit_truthMC_weight;
  IntVec   m_HcalHit_layer;

  //localMax
  TTree *t_LocalMax;
  int m_NlmU, m_NlmV;
  IntVec m_localMaxU_mc_pdg, m_localMaxV_mc_pdg, m_localMaxU_mc_tag, m_localMaxV_mc_tag;
  FloatVec m_localMaxU_tag, m_localMaxU_x, m_localMaxU_y, m_localMaxU_z, m_localMaxU_E,
            m_localMaxU_mc_px, m_localMaxU_mc_py, m_localMaxU_mc_pz, m_localMaxU_mc_weight;
  FloatVec m_localMaxV_tag, m_localMaxV_x, m_localMaxV_y, m_localMaxV_z, m_localMaxV_E,
            m_localMaxV_mc_px, m_localMaxV_mc_py, m_localMaxV_mc_pz, m_localMaxV_mc_weight;


  //1D showers
  TTree *t_Layers;
  int m_NshowerU, m_NshowerV;
  IntVec m_barShowerU_mc_pdg, m_barShowerV_mc_pdg, m_barShowerU_mc_tag, m_barShowerV_mc_tag;
  FloatVec m_barShowerU_tag, m_barShowerU_x, m_barShowerU_y, m_barShowerU_z, m_barShowerU_E,
            m_barShowerU_mc_px, m_barShowerU_mc_py, m_barShowerU_mc_pz, m_barShowerU_mc_weight;
  FloatVec m_barShowerV_tag, m_barShowerV_x, m_barShowerV_y, m_barShowerV_z, m_barShowerV_E,
            m_barShowerV_mc_px, m_barShowerV_mc_py, m_barShowerV_mc_pz, m_barShowerV_mc_weight;
  // Hough axis
  TTree * t_Hough;
  IntVec m_houghU_type, m_houghV_type;
  FloatVec  m_houghU_tag, m_houghU_x, m_houghU_y, m_houghU_z, m_houghU_E,
            m_houghU_truth_tag, m_houghU_truth_MC_px, m_houghU_truth_MC_py, m_houghU_truth_MC_pz, m_houghU_truth_MC_E, m_houghU_truth_MC_weight,
            m_houghU_hit_tag, m_houghU_hit_x, m_houghU_hit_y, m_houghU_hit_z, m_houghU_hit_E;
  FloatVec  m_houghV_tag, m_houghV_x, m_houghV_y, m_houghV_z, m_houghV_E, m_houghV_alpha, m_houghV_rho,
            m_houghV_truth_tag, m_houghV_truth_MC_px, m_houghV_truth_MC_py, m_houghV_truth_MC_pz, m_houghV_truth_MC_E, m_houghV_truth_MC_weight,
            m_houghV_hit_tag, m_houghV_hit_x, m_houghV_hit_y, m_houghV_hit_z, m_houghV_hit_E;
  // Cone axis
  TTree * t_Cone;
  IntVec m_coneU_type, m_coneV_type;
  FloatVec  m_coneU_tag, m_coneU_x, m_coneU_y, m_coneU_z, m_coneU_E,
            m_coneU_truth_tag, m_coneU_truth_MC_px, m_coneU_truth_MC_py, m_coneU_truth_MC_pz, m_coneU_truth_MC_E, m_coneU_truth_MC_weight,
            m_coneU_hit_tag, m_coneU_hit_x, m_coneU_hit_y, m_coneU_hit_z, m_coneU_hit_E;
  FloatVec  m_coneV_tag, m_coneV_x, m_coneV_y, m_coneV_z, m_coneV_E,
            m_coneV_truth_tag, m_coneV_truth_MC_px, m_coneV_truth_MC_py, m_coneV_truth_MC_pz, m_coneV_truth_MC_E, m_coneV_truth_MC_weight,
            m_coneV_hit_tag, m_coneV_hit_x, m_coneV_hit_y, m_coneV_hit_z, m_coneV_hit_E;
  // Track axis
  TTree * t_TrackAxis;
  IntVec m_trackU_type, m_trackV_type;
  FloatVec  m_trackU_tag, m_trackU_x, m_trackU_y, m_trackU_z, m_trackU_E,
            m_trackU_truth_tag, m_trackU_truth_MC_px, m_trackU_truth_MC_py, m_trackU_truth_MC_pz, m_trackU_truth_MC_E, m_trackU_truth_MC_weight,
            m_trackU_hit_tag, m_trackU_hit_x, m_trackU_hit_y, m_trackU_hit_z, m_trackU_hit_E;
  FloatVec  m_trackV_tag, m_trackV_x, m_trackV_y, m_trackV_z, m_trackV_E,
            m_trackV_truth_tag, m_trackV_truth_MC_px, m_trackV_truth_MC_py, m_trackV_truth_MC_pz, m_trackV_truth_MC_E, m_trackV_truth_MC_weight,
            m_trackV_hit_tag, m_trackV_hit_x, m_trackV_hit_y, m_trackV_hit_z, m_trackV_hit_E;
  // axis
  TTree *t_Axis;
  IntVec m_axisU_type, m_axisV_type;
  FloatVec  m_axisU_tag, m_axisU_x, m_axisU_y, m_axisU_z, m_axisU_E,
            m_axisU_truth_tag, m_axisU_truth_MC_px, m_axisU_truth_MC_py, m_axisU_truth_MC_pz, m_axisU_truth_MC_E, m_axisU_truth_MC_weight,
            m_axisU_hit_tag, m_axisU_hit_x, m_axisU_hit_y, m_axisU_hit_z, m_axisU_hit_E;
  FloatVec  m_axisV_tag, m_axisV_x, m_axisV_y, m_axisV_z, m_axisV_E,
            m_axisV_truth_tag, m_axisV_truth_MC_px, m_axisV_truth_MC_py, m_axisV_truth_MC_pz, m_axisV_truth_MC_E, m_axisV_truth_MC_weight,
            m_axisV_hit_tag, m_axisV_hit_x, m_axisV_hit_y, m_axisV_hit_z, m_axisV_hit_E;
  FloatVec  m_emptyAxisU_tag, m_emptyAxisU_x, m_emptyAxisU_y, m_emptyAxisU_z, m_emptyAxisU_E; 
  FloatVec  m_emptyAxisV_tag, m_emptyAxisV_x, m_emptyAxisV_y, m_emptyAxisV_z, m_emptyAxisV_E; 


  //HalfCluster after energy splitting
  TTree *t_HalfCluster;
  float m_totE_HFClusU, m_totE_HFClusV;
  FloatVec m_HalfClusterV_x, m_HalfClusterV_y, m_HalfClusterV_z, m_HalfClusterV_E, m_HalfClusterV_tag, m_HalfClusterV_type, m_HalfClusterV_nTrk;
  FloatVec m_HalfClusterV_hit_x, m_HalfClusterV_hit_y, m_HalfClusterV_hit_z, m_HalfClusterV_hit_E, m_HalfClusterV_hit_tag;
  FloatVec m_HalfClusterV_truth_tag, m_HalfClusterV_truthMC_px, m_HalfClusterV_truthMC_py, m_HalfClusterV_truthMC_pz, m_HalfClusterV_truthMC_E, m_HalfClusterV_truthMC_weight;
  FloatVec m_HalfClusterU_x, m_HalfClusterU_y, m_HalfClusterU_z, m_HalfClusterU_E, m_HalfClusterU_tag, m_HalfClusterU_type, m_HalfClusterU_nTrk;
  FloatVec m_HalfClusterU_hit_x, m_HalfClusterU_hit_y, m_HalfClusterU_hit_z, m_HalfClusterU_hit_E, m_HalfClusterU_hit_tag;
  FloatVec m_HalfClusterU_truth_tag, m_HalfClusterU_truthMC_px, m_HalfClusterU_truthMC_py, m_HalfClusterU_truthMC_pz, m_HalfClusterU_truthMC_E, m_HalfClusterU_truthMC_weight;


  //Tower
  TTree *t_Tower; 
  int m_Ntower;
  IntVec m_towerID_id1; 
  IntVec m_towerID_id2; 
  IntVec m_NclusU, m_NclusV;
  FloatVec m_totEn, m_totEn_U, m_totEn_V;
  

  //3D clusters
  TTree *t_Cluster;
  float m_totE_Ecal, m_totE_Hcal;
  int m_Nclus_Ecal, m_Nclus_Hcal;
  IntVec m_EcalClus_trk_location, m_EcalClus_trk_tag;
  FloatVec m_EcalClus_x, m_EcalClus_y, m_EcalClus_z, m_EcalClus_E, m_EcalClus_Escale, m_EcalClus_nTrk, m_EcalClus_pTrk, m_EcalClus_typeU, m_EcalClus_typeV,
          m_EcalClus_hitU_x, m_EcalClus_hitU_y, m_EcalClus_hitU_z, m_EcalClus_hitU_E, m_EcalClus_hitU_tag,
          m_EcalClus_hitV_x, m_EcalClus_hitV_y, m_EcalClus_hitV_z, m_EcalClus_hitV_E, m_EcalClus_hitV_tag,
          m_EcalClus_trk_d0, m_EcalClus_trk_z0, m_EcalClus_trk_phi, m_EcalClus_trk_tanL, m_EcalClus_trk_omega, m_EcalClus_trk_kappa,
          m_EcalClus_truthMC_tag, m_EcalClus_truthMC_pid, m_EcalClus_truthMC_px, m_EcalClus_truthMC_py, m_EcalClus_truthMC_pz, m_EcalClus_truthMC_E,
          m_EcalClus_truthMC_EPx, m_EcalClus_truthMC_EPy, m_EcalClus_truthMC_EPz, m_EcalClus_truthMC_weight;
  FloatVec m_HcalClus_x, m_HcalClus_y, m_HcalClus_z, m_HcalClus_E, m_HcalClus_nTrk, m_HcalClus_pTrk, m_HcalClus_nHit,
          m_HcalClus_hit_x, m_HcalClus_hit_y, m_HcalClus_hit_z, m_HcalClus_hit_E, m_HcalClus_hit_tag,
          m_HcalClus_truthMC_tag, m_HcalClus_truthMC_pid, m_HcalClus_truthMC_px, m_HcalClus_truthMC_py, m_HcalClus_truthMC_pz, m_HcalClus_truthMC_E,
          m_HcalClus_truthMC_EPx, m_HcalClus_truthMC_EPy, m_HcalClus_truthMC_EPz, m_HcalClus_truthMC_weight;
  FloatVec m_SimpleHcalClus_x, m_SimpleHcalClus_y, m_SimpleHcalClus_z, m_SimpleHcalClus_E, m_SimpleHcalClus_nTrk, m_SimpleHcalClus_pTrk, m_SimpleHcalClus_nHit,
          m_SimpleHcalClus_hit_x, m_SimpleHcalClus_hit_y, m_SimpleHcalClus_hit_z, m_SimpleHcalClus_hit_E, m_SimpleHcalClus_hit_tag,
          m_SimpleHcalClus_truthMC_tag, m_SimpleHcalClus_truthMC_pid, m_SimpleHcalClus_truthMC_px, m_SimpleHcalClus_truthMC_py, m_SimpleHcalClus_truthMC_pz, m_SimpleHcalClus_truthMC_E,
          m_SimpleHcalClus_truthMC_EPx, m_SimpleHcalClus_truthMC_EPy, m_SimpleHcalClus_truthMC_EPz, m_SimpleHcalClus_truthMC_weight;

  TTree *t_Track;
  int m_Ntrk; 
  FloatVec m_trk_p, m_trk_px, m_trk_py, m_trk_pz, m_trk_truthweight; 
  FloatVec m_trkstate_d0, m_trkstate_z0, m_trkstate_phi, m_trkstate_tanL, m_trkstate_omega, m_trkstate_kappa;
  FloatVec m_trkstate_refx, m_trkstate_refy, m_trkstate_refz; 
  IntVec m_trkstate_tag, m_trkstate_location, m_type, m_Nhit, m_pid, m_pid_truth;
  FloatVec m_trkstate_x_ECAL, m_trkstate_y_ECAL, m_trkstate_z_ECAL,  m_trkstate_x_HCAL, m_trkstate_y_HCAL, m_trkstate_z_HCAL;
  IntVec m_trkstate_tag_ECAL, m_trkstate_tag_HCAL;

  // yyy: output PFO information
  TTree *t_PFO;
  IntVec pfo_tag, pfo_n_track, pfo_n_ecal_clus, pfo_n_hcal_clus, pfo_ecal_tag, pfo_hcal_tag;
  IntVec pfo_trk_location, pfo_trk_tag;
  FloatVec pfo_trk_d0, pfo_trk_z0, pfo_trk_phi, pfo_trk_tanL, pfo_trk_omega, pfo_trk_kappa;
  FloatVec  pfo_ecal_clus_x, pfo_ecal_clus_y, pfo_ecal_clus_z, pfo_ecal_clus_E, pfo_ecal_clus_Escale,
            pfo_hcal_clus_x, pfo_hcal_clus_y, pfo_hcal_clus_z, pfo_hcal_clus_E;


  void ClearMCParticle();
  void ClearBar();
  void ClearLocalMax();
  void ClearLayer();
  void ClearHough();
  void ClearCone();
  void ClearTrackAxis();
  void ClearAxis();
  void ClearHalfCluster();
  void ClearTower();
  void ClearCluster();
  void ClearTrack();
  void ClearPFO(); // yyy

  double GetParticleDepEnergy(edm4hep::MCParticle& mcp, std::vector<std::shared_ptr<Cyber::CaloUnit>>& barcol);
  double GetParticleDepEnergy(edm4hep::MCParticle& mcp, std::vector<std::shared_ptr<Cyber::CaloHit>>& hitcol);
};
#endif
