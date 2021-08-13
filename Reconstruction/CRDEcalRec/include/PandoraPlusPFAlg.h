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

  TTree* t_ArborTree;
  int m_Nnode;
  float m_tree_centx, m_tree_centy, m_tree_centz, m_tree_minLayer, m_tree_maxLayer; 
  FloatVec m_nodex, m_nodey, m_nodez, m_nodeType, m_nodeE;

  TTree* t_dataColIter0;
  int      m_Iter0_Ngoodclus, m_Iter0_Nbadclus;
  IntVec   m_Iter0_clus_Nly, m_Iter0_clus_start, m_Iter0_clus_end; 
  FloatVec m_Iter0_clus_x, m_Iter0_clus_y, m_Iter0_clus_z, m_Iter0_clus_E;
  FloatVec m_Iter0_clus_px, m_Iter0_clus_py, m_Iter0_clus_pz, m_Iter0_clus_chi2, m_Iter0_clus_aveE, m_Iter0_clus_stdDevE;
  FloatVec m_Iter0_clus_alpha, m_Iter0_clus_beta, m_Iter0_clus_showermax;   
  FloatVec m_Iter0_gclus_2dshx, m_Iter0_gclus_2dshy, m_Iter0_gclus_2dshz, m_Iter0_gclus_2dshE;

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
  void ClearRecPFO(); 
  void ClearTree();
  void ClearIter();

};
#endif
