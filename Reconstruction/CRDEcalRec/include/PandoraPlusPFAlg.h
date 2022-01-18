#ifndef PANDORAPLUS_ALG_H
#define PANDORAPLUS_ALG_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "PandoraPlusDataCol.h"
#include "CRDEcalSvc/ICRDEcalSvc.h"
#include "DetInterface/IGeomSvc.h"
#include "Tools/MCParticleCreator.h"
#include "Tools/TrackCreator.h"
#include "Tools/VertexCreator.h"
#include "Tools/EcalHitsCreator.h"
#include "Tools/HcalHitsCreator.h"
#include "Tools/PFOCreator.h"

#include "EcalClusterReconstruction.h"

#include "TVector3.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
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
  //float _seed = 1024;


  //Services
  SmartIF<ICRDEcalSvc> m_edmsvc;


  //DataCollection
  PandoraPlusDataCol     m_DataCol;


  //Creators and their setting
  MCParticleCreator       *m_pMCParticleCreator;
  TrackCreator            *m_pTrackCreator;
  VertexCreator           *m_pVertexCreator;
  EcalHitsCreator         *m_pEcalHitsCreator;
  HcalHitsCreator         *m_pHcalHitsCreator;
  PFOCreator              *m_pPfoCreator;   

  MCParticleCreator::Settings      *m_pMCParticleCreatorSettings;
  TrackCreator::Settings           *m_pTrackCreatorSettings;
  VertexCreator::Settings          *m_pVertexCreatorSettings;
  EcalHitsCreator::Settings        *m_EcalHitsCreatorSettings;
  HcalHitsCreator::Settings        *m_pHcalHitsCreatorSettings;
  PFOCreator::Settings             *m_pPfoCreatorSettings;

  //Algorithm for PFA
  EcalClusterReconstruction   *m_pEcalClusterRec;


  //Parameters for PFA algorithm
  EcalClusterReconstruction::Settings       *m_pEcalClusterRecSettings;



  //Input collection
  DataHandle<edm4hep::MCParticleCollection> r_MCParticleCol{"MCParticle", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackCollection>      r_MarlinTrkCol{"MarlinTrkTracks", Gaudi::DataHandle::Reader, this};

  //Parameters for algorithm settings.
  //TODO: use a xml file to readin all parameters. 
  mutable Gaudi::Property<std::string> _filename{this, "OutFileName", "testout.root", "Output file name"};
  mutable Gaudi::Property<float> _seed{this,   "Seed", 2131, "Random Seed"};
  mutable Gaudi::Property<int>  _Debug{this,   "Debug", 0, "Debug level"};
  mutable Gaudi::Property<int>  _Nskip{this,   "SkipEvt", 0, "Skip event"};
  //Gaudi::Property< std::vector<std::string> >  m_MCParticleCollections{ this, "MCParticleCollections", {"MCParticle"} };

  // Output collections
  // output: PFOs 




  //PFA input end here. 
  //-----------------------------------------------------------
  //Followings are for code testing. 
  typedef std::vector<float> FloatVec;
  typedef std::vector<int>   IntVec;
  TFile* m_wfile;
  //Stage0: check input simbars. 
  TTree* t_SimBar;
  FloatVec m_simBar_x, m_simBar_y, m_simBar_z, m_simBar_T1, m_simBar_T2, m_simBar_Q1, m_simBar_Q2, m_simBar_dlayer, m_simBar_part, m_simBar_stave, m_simBar_slayer, m_simBar_module;

  //Check local max
  TTree *t_Layers;
  int m_NshowerX, m_NshowerY;
  FloatVec m_barShowerX_x, m_barShowerX_y, m_barShowerX_z, m_barShowerX_E;
  FloatVec m_barShowerY_x, m_barShowerY_y, m_barShowerY_z, m_barShowerY_E; 
 
  //Check Hough clusters
  TTree *t_HoughClusters;
  int m_NclusX, m_NclusY;
  FloatVec m_clusX_x, m_clusX_y, m_clusX_z, m_clusX_E, m_clusX_alpha, m_clusX_rho, m_clusX_px, m_clusX_py, m_clusX_pz;
  FloatVec m_clusY_x, m_clusY_y, m_clusY_z, m_clusY_E, m_clusY_alpha, m_clusY_rho, m_clusY_px, m_clusY_py, m_clusY_pz;
  IntVec m_clusX_Nhit, m_clusY_Nhit;

  //Stage2: check reconstructed result
  TTree *t_recoPFO;
  int m_Npfo, m_Nmc, m_N3dclus; 
  FloatVec m_recPFO_px, m_recPFO_py, m_recPFO_pz, m_recPFO_En;
  IntVec m_recPFO_pid;
  FloatVec m_Clus_x, m_Clus_y, m_Clus_z, m_Clus_E; 
  IntVec m_N2dshInClus, m_Clus_type;
  IntVec m_mcPdgid, m_mcNdaughter, m_mcNparent, m_mcStatus; 
  FloatVec m_mcPx, m_mcPy, m_mcPz, m_mcEn;
  FloatVec m_scndM; 


  void ClearBar();
  void ClearLayer();
  void ClearCluster();
  void ClearRecPFO(); 

};
#endif
