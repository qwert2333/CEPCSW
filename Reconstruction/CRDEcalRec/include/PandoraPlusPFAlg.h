#ifndef PANDORAPLUS_ALG_H
#define PANDORAPLUS_ALG_H

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
#include "Algorithm/HoughClusteringAlg.h"
#include "Algorithm/EnergySplittingAlg.h"
#include "Algorithm/EnergyTimeMatchingAlg.h"

#include "TVector3.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

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


  //DataCollection
  PandoraPlusDataCol     m_DataCol;   //TODO: decide if set it to const. 


  //Creators and their setting
  MCParticleCreator       *m_pMCParticleCreator;
  TrackCreator            *m_pTrackCreator;
  CaloHitsCreator         *m_pCaloHitsCreator;
  OutputCreator           *m_pOutputCreator;
  

  MCParticleCreator::Settings   m_pMCParticleCreatorSettings;
  TrackCreator::Settings        m_pTrackCreatorSettings;
  CaloHitsCreator::Settings     m_CaloHitsCreatorSettings;
  OutputCreator::Settings       m_OutputCreatorSettings;

  //Parameters for PFA algorithm
  Settings m_GlobalSettings; 
  //std::map<std::string, Settings> m_algorithmSettings; 

  //Algorithm for PFA
  PandoraPlus::AlgorithmManager m_algorithmManager; 


  //Readin collection names
  Gaudi::Property< std::string > name_MCParticleCol{ this, "MCParticleCollection", "MCParticle" };
  Gaudi::Property< std::vector<std::string> > name_TrackCol{ this, "TrackCollections", {"MarlinTrkTracks"} };
  Gaudi::Property< std::vector<std::string> > name_EcalHits{ this, "ECalCaloHitCollections", {"ECALBarrel"} };
  Gaudi::Property< std::vector<std::string> > name_EcalReadout{ this, "ECalReadOutNames", {"EcalBarrelCollection"} }; 
  Gaudi::Property< std::vector<std::string> > name_HcalHits{ this, "HCalCaloHitCollections", {"HCALBarrel"} };
  Gaudi::Property< std::vector<std::string> > name_HcalReadout{ this, "HCalReadOutNames", {"HcalBarrelCollection"} }; 

  //---Readin collections
  typedef DataHandle<edm4hep::TrackCollection>           TrackType; 
  typedef DataHandle<edm4hep::CalorimeterHitCollection>  CaloType; 
  DataHandle<edm4hep::MCParticleCollection>* r_MCParticleCol; 
  std::vector<TrackType*> r_TrackCols; 
  std::vector<CaloType*>  r_ECalHitCols; 
  std::vector<CaloType*>  r_HCalHitCols; 
  std::vector<CaloType*>  r_CaloHitCols; 


  //Global parameters.
  Gaudi::Property<float> m_BField{this,   "BField", 3., "Magnetic field"};
  Gaudi::Property<float> m_seed{this,   "Seed", 2131, "Random Seed"};
  Gaudi::Property<int>   m_Debug{this,   "Debug", 0, "Debug level"};
  Gaudi::Property<int>   m_Nskip{this,   "SkipEvt", 0, "Skip event"};


  
  //Algorithms: 
  typedef std::vector<double> FloatVector;
  typedef std::vector<std::string> StringVector;
  Gaudi::Property< StringVector > name_Algs{ this, "AlgList", {} };
  Gaudi::Property< std::vector<StringVector> > name_AlgPars{ this, "AlgParNames", {} };
  Gaudi::Property< std::vector<FloatVector> > value_AlgPars{this, "AlgParValues", {} };


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
  FloatVec m_simBar_x, m_simBar_y, m_simBar_z, m_simBar_T1, m_simBar_T2, m_simBar_Q1, m_simBar_Q2, m_simBar_dlayer, m_simBar_part, m_simBar_stave, m_simBar_slayer, m_simBar_module;

  TTree *t_Layers;
  int m_NshowerX, m_NshowerY;
  FloatVec m_barShowerX_x, m_barShowerX_y, m_barShowerX_z, m_barShowerX_E;
  FloatVec m_barShowerY_x, m_barShowerY_y, m_barShowerY_z, m_barShowerY_E;

  TTree *t_Cluster;
  int m_Nclus;
  FloatVec m_Clus_x, m_Clus_y, m_Clus_z, m_Clus_E;
  IntVec m_Nhit;

  //check neighbor clustering
  TTree* t_Clustering;
  int m_3dcluster, m_2dcluster, m_1dcluster, m_barcluster, m_bar; //efficiency
  FloatVec m_E_3dcluster, m_E_2dcluster, m_E_1dcluster, m_E_barcluster, m_E_bar; //resolution
  FloatVec m_bar_tag, m_bar_energy, m_bar_dlayer, m_bar_slayer, m_bar_x, m_bar_y, m_bar_z, m_bar_module, m_bar_part, m_bar_stave, m_bar_bar; //distribution check

  void ClearBar();
  void ClearLayer();
  void ClearCluster();
  void ClearClustering();
};
#endif
