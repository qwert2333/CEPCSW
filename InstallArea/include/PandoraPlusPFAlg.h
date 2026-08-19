#ifndef PANDORAPLUS_ALG_H
#define PANDORAPLUS_ALG_H

#include <string>
#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/Segmentations.h>
#include "DetInterface/IGeomSvc.h"

#include "PandoraPlusDataCol.h"
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
//#include "Algorithm/ConeClusteringAlg.h"
#include "Algorithm/TrackExtrapolatingAlg.h"
#include "Algorithm/PFOCreatingAlg.h"

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

using namespace PandoraPlus;
using namespace std;
class PandoraPlusPFAlg : public GaudiAlgorithm
{
 
public:
 
  PandoraPlusPFAlg(const std::string& name, ISvcLocator* svcLoc);
 
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
  std::map<std::string, dd4hep::DDSegmentation::BitFieldCoder*> map_readout_decoder;

  //DataCollection: moved into execute() to ensure everything can be cleand after one event. 
  //PandoraPlusDataCol     m_DataCol; 


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
  PandoraPlus::AlgorithmManager m_algorithmManager; 


  //Readin collection names
  Gaudi::Property< std::string > name_MCParticleCol{ this, "MCParticleCollection", "MCParticle" };
  Gaudi::Property< std::string > name_MCPTrkAssoCol{this, "MCRecoTrackParticleAssociationCollection", "MarlinTrkAssociation"};
  Gaudi::Property< std::vector<std::string> > name_TrackCol{ this, "TrackCollections", {"MarlinTrkTracks"} };
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

  //Global parameters.
  Gaudi::Property<float>  m_BField{this,  "BField", 3., "Magnetic field"};
  Gaudi::Property<float>  m_seed{this,    "Seed", 2131, "Random Seed"};
  Gaudi::Property<int>    m_Debug{this,   "Debug", 0, "Debug level"};
  Gaudi::Property<int>    m_Nskip{this,   "SkipEvt", 0, "Skip event"};
  Gaudi::Property<std::string>   m_EcalType{this, "EcalType", "BarEcal", "ECAL type"};
  Gaudi::Property<bool>   m_useTruthTrk{this,   "UseTruthTrack", 1, "Use truth track or reconstructed track"};

  
  //Algorithms: 
  typedef std::vector<std::string> StringVector;
  Gaudi::Property< StringVector > name_Algs{ this, "AlgList", {} };
  Gaudi::Property< std::vector<StringVector> > name_AlgPars{ this, "AlgParNames", {} };
  Gaudi::Property< std::vector<StringVector> > type_AlgPars{ this, "AlgParTypes", {} };
  Gaudi::Property< std::vector<StringVector> > value_AlgPars{this, "AlgParValues", {} };


  // Output collections
  DataHandle<edm4hep::CalorimeterHitCollection>             w_RecCaloCol{"RecECALBarrel", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::ClusterCollection>                    w_ClusterCollection {"PandoraClusters",Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::ReconstructedParticleCollection>      w_ReconstructedParticleCollection {"PandoraPFOs"    ,Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::VertexCollection>                     w_VertexCollection {"PandoraPFANewStartVertices",Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoParticleAssociationCollection>  w_MCRecoParticleAssociationCollection {"pfoMCRecoParticleAssociation",Gaudi::DataHandle::Writer, this};


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
  FloatVec m_mcPx, m_mcPy, m_mcPz, m_mcEn, m_mcMass, m_mcCharge, m_mcEPx, m_mcEPy, m_mcEPz;

  //Raw bars and hits
  TTree* t_SimBar;
  float m_totE_EcalSim, m_totE_HcalSim;
  FloatVec m_simBar_x, m_simBar_y, m_simBar_z, m_simBar_T1, m_simBar_T2, m_simBar_Q1, m_simBar_Q2; 
  FloatVec m_simBar_truthMC_tag, m_simBar_truthMC_pid, m_simBar_truthMC_px, m_simBar_truthMC_py, m_simBar_truthMC_pz, m_simBar_truthMC_E, 
           m_simBar_truthMC_EPx, m_simBar_truthMC_EPy, m_simBar_truthMC_EPz, m_simBar_truthMC_weight;
  IntVec m_simBar_dlayer, m_simBar_part, m_simBar_stave, m_simBar_slayer, m_simBar_module, m_simBar_bar;
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

  //axis
  TTree *t_axis;
  int m_NaxisU, m_NaxisV;
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
  int towerID[3];
  int m_NclusU, m_NclusV; 
  float m_totEn; 


  //3D clusters
  TTree *t_Cluster;
  float m_totE_Ecal, m_totE_Hcal;
  int m_Nclus_Ecal, m_Nclus_Hcal; 
  FloatVec m_EcalClus_x, m_EcalClus_y, m_EcalClus_z, m_EcalClus_E, m_EcalClus_nTrk, m_EcalClus_pTrk, m_EcalClus_typeU, m_EcalClus_typeV,
           m_EcalClus_hitU_x, m_EcalClus_hitU_y, m_EcalClus_hitU_z, m_EcalClus_hitU_E, m_EcalClus_hitU_tag,
           m_EcalClus_hitV_x, m_EcalClus_hitV_y, m_EcalClus_hitV_z, m_EcalClus_hitV_E, m_EcalClus_hitV_tag,
           m_EcalClus_truthMC_tag, m_EcalClus_truthMC_pid, m_EcalClus_truthMC_px, m_EcalClus_truthMC_py, m_EcalClus_truthMC_pz, m_EcalClus_truthMC_E, 
           m_EcalClus_truthMC_EPx, m_EcalClus_truthMC_EPy, m_EcalClus_truthMC_EPz, m_EcalClus_truthMC_weight;
  FloatVec m_HcalClus_x, m_HcalClus_y, m_HcalClus_z, m_HcalClus_E, m_HcalClus_nTrk, m_HcalClus_pTrk,
           m_HcalClus_truthMC_tag, m_HcalClus_truthMC_pid, m_HcalClus_truthMC_px, m_HcalClus_truthMC_py, m_HcalClus_truthMC_pz, m_HcalClus_truthMC_E, 
           m_HcalClus_truthMC_EPx, m_HcalClus_truthMC_EPy, m_HcalClus_truthMC_EPz, m_HcalClus_truthMC_weight;
  FloatVec m_SimpleHcalClus_x, m_SimpleHcalClus_y, m_SimpleHcalClus_z, m_SimpleHcalClus_E, m_SimpleHcalClus_nTrk, m_SimpleHcalClus_pTrk,
           m_SimpleHcalClus_truthMC_tag, m_SimpleHcalClus_truthMC_pid, m_SimpleHcalClus_truthMC_px, m_SimpleHcalClus_truthMC_py, m_SimpleHcalClus_truthMC_pz, m_SimpleHcalClus_truthMC_E,
           m_SimpleHcalClus_truthMC_EPx, m_SimpleHcalClus_truthMC_EPy, m_SimpleHcalClus_truthMC_EPz, m_SimpleHcalClus_truthMC_weight;

  void ClearMCParticle();
  void ClearBar();
  void ClearLocalMax();
  void ClearLayer();
  void ClearAxis();
  void ClearHalfCluster();
  void ClearTower();
  void ClearCluster();
};
#endif
