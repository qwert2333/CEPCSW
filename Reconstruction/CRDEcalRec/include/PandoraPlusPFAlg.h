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
  TTree* t_SimBar;
  float m_totE_EcalSim, m_totE_HcalSim;
  FloatVec m_simBar_x, m_simBar_y, m_simBar_z, m_simBar_T1, m_simBar_T2, m_simBar_Q1, m_simBar_Q2; 
  IntVec m_simBar_dlayer, m_simBar_part, m_simBar_stave, m_simBar_slayer, m_simBar_module, m_simBar_bar;
  FloatVec m_HcalHit_x, m_HcalHit_y, m_HcalHit_z, m_HcalHit_E;
  IntVec   m_HcalHit_layer;

  TTree *t_HalfCluster;
  float m_totE_HFClusU, m_totE_HFClusV;
  FloatVec m_HalfClusterV_x, m_HalfClusterV_y, m_HalfClusterV_z, m_HalfClusterV_E;
  FloatVec m_HalfClusterU_x, m_HalfClusterU_y, m_HalfClusterU_z, m_HalfClusterU_E;

  TTree *t_Cluster;
  int m_Nclus, m_Nmc;
  float m_totE_EcalRec, m_totE_HcalRec;
  IntVec m_Clus_Ntrk, m_Clus_Nhit, m_Clus_truthPDG;
  IntVec m_Clus_startLayer, m_Clus_endLayer, m_Clus_maxELayer, m_Clus_maxWidthLayer, m_Clus_typeU, m_Clus_typeV;
  FloatVec m_Clus_x, m_Clus_y, m_Clus_z, m_Clus_E, m_Clus_Px, m_Clus_Py, m_Clus_Pz, m_Clus_Ptrk, m_Clus_truthFrac;
  FloatVec m_Clus_width, m_Clus_ScndM, m_Clus_E1Etot, m_Clus_E2Etot, m_Clus_E5Etot, m_Clus_EhalfEtot, m_Clus_EaxisEtot;
  FloatVec m_Clus_hitx, m_Clus_hity, m_Clus_hitz, m_Clus_hitE, m_Clus_hittag, m_Clus_hittag_trk;
  IntVec m_mcPdgid, m_mcStatus;
  FloatVec m_mcPx, m_mcPy, m_mcPz, m_mcEn, m_mcMass, m_mcCharge, m_mcEPx, m_mcEPy, m_mcEPz;
  FloatVec m_Hcal_clus_x, m_Hcal_clus_y, m_Hcal_clus_z, m_Hcal_clus_E;
  FloatVec m_Hcal_hit_tag, m_Hcal_hit_x, m_Hcal_hit_y, m_Hcal_hit_z, m_Hcal_hit_E;


  TTree *t_Track;
  int m_Ntrk; 
  FloatVec m_trkstate_d0, m_trkstate_z0, m_trkstate_phi, m_trkstate_tanL, m_trkstate_omega, m_trkstate_kappa;
  FloatVec m_trkstate_refx, m_trkstate_refy, m_trkstate_refz; 
  IntVec m_trkstate_tag, m_trkstate_location, m_type;


  void ClearBar();
  void ClearCluster();
  void ClearTrack();
  void ClearHalfCluster();


};
#endif
