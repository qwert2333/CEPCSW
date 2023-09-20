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
#include "Algorithm/LocalMaxFindingAlg.h"
#include "Algorithm/TrackMatchingAlg.h"
#include "Algorithm/HoughClusteringAlg.h"
#include "Algorithm/ConeClustering2DAlg.h"
#include "Algorithm/AxisMergingAlg.h"
#include "Algorithm/EnergySplittingAlg.h"
#include "Algorithm/EnergyTimeMatchingAlg.h"
//#include "Algorithm/ConeClusteringAlg.h"
#include "Algorithm/TrackExtrapolatingAlg.h"

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
  Gaudi::Property< std::string > name_MCPRecoCaloAssoCol{ this, "MCRecoCaloParticleAssociationCollection", "MCRecoCaloParticleAssociationCollection" };
  Gaudi::Property< std::vector<std::string> > name_TrackCol{ this, "TrackCollections", {"MarlinTrkTracks"} };
  Gaudi::Property< std::vector<std::string> > name_EcalHits{ this, "ECalCaloHitCollections", {"ECALBarrel"} };
  Gaudi::Property< std::vector<std::string> > name_EcalReadout{ this, "ECalReadOutNames", {"EcalBarrelCollection"} }; 
  Gaudi::Property< std::vector<std::string> > name_HcalHits{ this, "HCalCaloHitCollections", {"HCALBarrel"} };
  Gaudi::Property< std::vector<std::string> > name_HcalReadout{ this, "HCalReadOutNames", {"HcalBarrelCollection"} }; 
  

  //---Readin collections
  typedef DataHandle<edm4hep::TrackCollection>           TrackType; 
  typedef DataHandle<edm4hep::CalorimeterHitCollection>  CaloType; 
  DataHandle<edm4hep::MCParticleCollection>* r_MCParticleCol; 
  DataHandle<edm4hep::MCRecoCaloParticleAssociationCollection>* r_MCPRecoCaloAssoCol; 
  std::vector<TrackType*> r_TrackCols; 
  //std::vector<CaloType*>  r_ECalHitCols; 
  //std::vector<CaloType*>  r_HCalHitCols; 
  std::vector<CaloType*>  r_CaloHitCols; 


  //Global parameters.
  Gaudi::Property<float>  m_BField{this,  "BField", 3., "Magnetic field"};
  Gaudi::Property<float>  m_seed{this,    "Seed", 2131, "Random Seed"};
  Gaudi::Property<int>    m_Debug{this,   "Debug", 0, "Debug level"};
  Gaudi::Property<int>    m_Nskip{this,   "SkipEvt", 0, "Skip event"};
  Gaudi::Property<std::string>   m_EcalType{this, "EcalType", "BarEcal", "ECAL type"};

  
  //Algorithms: 
  typedef std::vector<std::string> StringVector;
  Gaudi::Property< StringVector > name_Algs{ this, "AlgList", {} };
  Gaudi::Property< std::vector<StringVector> > name_AlgPars{ this, "AlgParNames", {} };
  Gaudi::Property< std::vector<StringVector> > type_AlgPars{ this, "AlgParTypes", {} };
  Gaudi::Property< std::vector<StringVector> > value_AlgPars{this, "AlgParValues", {} };


  // Output collections
  DataHandle<edm4hep::CalorimeterHitCollection>         w_RecCaloCol{"RecECALBarrel", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::ClusterCollection>                w_ClusterCollection {"PandoraClusters",Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::ReconstructedParticleCollection>  w_ReconstructedParticleCollection {"PandoraPFOs"    ,Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::VertexCollection>                 w_VertexCollection {"PandoraPFANewStartVertices",Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoParticleAssociationCollection>  w_MCRecoParticleAssociationCollection {"pfoMCRecoParticleAssociation",Gaudi::DataHandle::Writer, this};


  //For Ana
  Gaudi::Property<bool>  m_WriteAna {this, "WriteAna", false, "Write Ntuples for analysis"};
  Gaudi::Property<std::string> m_filename{this, "AnaFileName", "testout.root", "Output file name"};

  typedef std::vector<float> FloatVec;
  typedef std::vector<int>   IntVec;

  TFile* m_wfile;
  TTree* t_SimBar;
  FloatVec m_simBar_x, m_simBar_y, m_simBar_z, m_simBar_T1, m_simBar_T2, m_simBar_Q1, m_simBar_Q2; 
  IntVec m_simBar_dlayer, m_simBar_part, m_simBar_stave, m_simBar_slayer, m_simBar_module, m_simBar_bar;
  FloatVec m_HcalHit_x, m_HcalHit_y, m_HcalHit_z, m_HcalHit_E;
  IntVec   m_HcalHit_layer;

  TTree *t_Layers;
  int m_NshowerU, m_NshowerV;
  FloatVec m_barShowerU_tag, m_barShowerU_x, m_barShowerU_y, m_barShowerU_z, m_barShowerU_E, m_barShowerU_T, m_barShowerU_module, m_barShowerU_part, m_barShowerU_stave, m_barShowerU_dlayer, m_barShowerU_slayer, m_barShowerU_bar;
  FloatVec m_barShowerV_tag, m_barShowerV_x, m_barShowerV_y, m_barShowerV_z, m_barShowerV_E, m_barShowerV_T, m_barShowerV_module, m_barShowerV_part, m_barShowerV_stave, m_barShowerV_dlayer, m_barShowerV_slayer, m_barShowerV_bar;

  // yyy: HalfClusters
  TTree *t_HalfCluster;
  FloatVec m_HalfClusterV_tag, m_HalfClusterV_x, m_HalfClusterV_y, m_HalfClusterV_z, m_HalfClusterV_E;
  FloatVec m_HalfClusterU_tag, m_HalfClusterU_x, m_HalfClusterU_y, m_HalfClusterU_z, m_HalfClusterU_E;

  // yyy: check Hough Algorithm
  TTree *t_Hough;
  FloatVec houghV_cluster_tag, houghV_axis_tag, houghV_axis_x, houghV_axis_y, houghV_axis_z, houghV_axis_E;
  FloatVec houghU_cluster_tag, houghU_axis_tag, houghU_axis_x, houghU_axis_y, houghU_axis_z, houghU_axis_E;

  // yyy: check TrackMatchingAlg
  TTree * t_Match;
  FloatVec matchV_cluster_tag, matchV_track_axis_tag, matchV_track_axis_x, matchV_track_axis_y, matchV_track_axis_z, matchV_track_axis_E; 
  FloatVec matchU_cluster_tag, matchU_track_axis_tag, matchU_track_axis_x, matchU_track_axis_y, matchU_track_axis_z, matchU_track_axis_E;

  // yyy: check ConeMatchingAlg
  TTree * t_Cone;
  FloatVec coneV_cluster_tag, coneV_axis_tag, coneV_axis_x, coneV_axis_y, coneV_axis_z, coneV_axis_E;
  FloatVec coneU_cluster_tag, coneU_axis_tag, coneU_axis_x, coneU_axis_y, coneU_axis_z, coneU_axis_E;

  // yyy: check MergedAxis (detailed information)
  TTree * t_Merge;
  IntVec mergeV_axis_index, mergeU_axis_index;
  FloatVec mergeV_axis_hit_x, mergeV_axis_hit_y, mergeV_axis_hit_z, mergeV_axis_hit_E;
  FloatVec mergeU_axis_hit_x, mergeU_axis_hit_y, mergeU_axis_hit_z, mergeU_axis_hit_E;
  FloatVec mergeV_cluster_E, mergeV_axis_E, mergeV_istrk, mergeV_type, mergeV_truthFracMax;
  FloatVec mergeU_cluster_E, mergeU_axis_E, mergeU_istrk, mergeU_type, mergeU_truthFracMax;
  IntVec mergeV_truthMaxPDGID, mergeU_truthMaxPDGID; // yyy: the PDG id of the particle with the max fragment in an axis

  TTree *t_axis;
  int m_NaxisU, m_NaxisV;
  IntVec m_axisU_Nhit, m_axisU_type, m_axisU_Nmcp, m_axisV_Nhit, m_axisV_type, m_axisV_Nmcp;
  FloatVec m_axisU_x, m_axisU_y, m_axisU_z, m_axisU_E, m_axisU_truthFracMax, m_axisV_x, m_axisV_y, m_axisV_z, m_axisV_E, m_axisV_truthFracMax;
  IntVec m_axisU_truthMaxPDGID, m_axisV_truthMaxPDGID;  

  TTree *t_Cluster;
  int m_Nclus, m_Nmc;
  float m_totE;
  IntVec m_Clus_Ntrk, m_Clus_Nhit;
  FloatVec m_Clus_x, m_Clus_y, m_Clus_z, m_Clus_E, m_Clus_Px, m_Clus_Py, m_Clus_Pz, m_Clus_Ptrk;
  FloatVec m_Clus_hitx, m_Clus_hity, m_Clus_hitz, m_Clus_hitE, m_Clus_hittag, m_Clus_hittag_trk;
  IntVec m_mcPdgid, m_mcStatus;
  FloatVec m_mcPx, m_mcPy, m_mcPz, m_mcEn, m_mcMass, m_mcCharge, m_mcEPx, m_mcEPy, m_mcEPz;


  TTree *t_Tower;
  int m_module, m_part, m_stave;
  IntVec m_HFClusU_Nhit, m_HFClusU_type, m_HFClusV_Nhit, m_HFClusV_type;
  FloatVec m_HFClusU_x, m_HFClusU_y, m_HFClusU_z, m_HFClusU_E, m_HFClusV_x, m_HFClusV_y, m_HFClusV_z, m_HFClusV_E;
  FloatVec m_HFClusUhit_tag, m_HFClusUhit_type, m_HFClusUhit_x, m_HFClusUhit_y, m_HFClusUhit_z, m_HFClusUhit_E;
  FloatVec m_HFClusVhit_tag, m_HFClusVhit_type, m_HFClusVhit_x, m_HFClusVhit_y, m_HFClusVhit_z, m_HFClusVhit_E;  

  TTree *t_Track;
  int m_Ntrk; 
  FloatVec m_trkstate_d0, m_trkstate_z0, m_trkstate_phi, m_trkstate_tanL, m_trkstate_omega, m_trkstate_kappa;
  FloatVec m_trkstate_refx, m_trkstate_refy, m_trkstate_refz; 
  IntVec m_trkstate_tag, m_trkstate_location, m_type;


  void ClearBar();
  void ClearLayer();
  void ClearHalfCluster(); // yyy
  void ClearHough(); // yyy
  void ClearMatch(); // yyy
  void ClearCone();  // yyy
  void ClearMerge(); // yyy
  void ClearTower();
  void ClearCluster();
  void ClearTrack();
  void ClearAxis();



};
#endif
