#ifndef PANDORAPLUS_ALG_H
#define PANDORAPLUS_ALG_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "PandoraPlusDataCol.h"
#include "CRDEcalEDMSvc/ICRDEcalEDMSvc.h"
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
  SmartIF<ICRDEcalEDMSvc> m_edmsvc;


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

  //Stage1: check intermediate info.
  TTree* t_dataColIter0; 
  FloatVec m_Iter0_ly_dlayer, m_Iter0_ly_part, m_Iter0_ly_stave, m_Iter0_ly_module, m_Iter0_ly_x, m_Iter0_ly_y, m_Iter0_ly_z;
  IntVec   m_Iter0_ly_NshowerX, m_Iter0_ly_NshowerY, m_Iter0_ly_NclusterX, m_Iter0_ly_NclusterY, m_Iter0_Nexpsh;
  FloatVec m_Iter0_barx, m_Iter0_bary, m_Iter0_barz, m_Iter0_barE; 
  IntVec   m_Iter0_barslayer;
  FloatVec m_Iter0_Expsh_x, m_Iter0_Expsh_y, m_Iter0_Expsh_z, m_Iter0_Expsh_E; 
  FloatVec m_Iter0_shower0E, m_Iter0_shower0X, m_Iter0_shower0Y, m_Iter0_shower0Z, m_Iter0_shower0T1, m_Iter0_shower0T2;
  FloatVec m_Iter0_shower1E, m_Iter0_shower1X, m_Iter0_shower1Y, m_Iter0_shower1Z, m_Iter0_shower1T1, m_Iter0_shower1T2;

  int      m_Iter0_N2dshower;
  IntVec   m_Iter0_2ds_dlayer, m_Iter0_2ds_part, m_Iter0_2ds_stave;
  FloatVec m_Iter0_2ds_x, m_Iter0_2ds_y, m_Iter0_2ds_z, m_Iter0_2ds_E;

  int      m_Iter0_Ngoodclus, m_Iter0_Nbadclus;
  IntVec   m_Iter0_clus_Nly; 
  FloatVec m_Iter0_clus_x, m_Iter0_clus_y, m_Iter0_clus_z, m_Iter0_clus_E;
  FloatVec m_Iter0_clus_px, m_Iter0_clus_py, m_Iter0_clus_pz; 


  TTree* t_dataColIter1;
  FloatVec m_Iter1_ly_dlayer, m_Iter1_ly_part, m_Iter1_ly_stave, m_Iter1_ly_module, m_Iter1_ly_x, m_Iter1_ly_y, m_Iter1_ly_z;
  IntVec   m_Iter1_ly_NshowerX, m_Iter1_ly_NshowerY, m_Iter1_ly_NclusterX, m_Iter1_ly_NclusterY, m_Iter1_Nexpsh;
  FloatVec m_Iter1_barx, m_Iter1_bary, m_Iter1_barz, m_Iter1_barE; 
  IntVec   m_Iter1_barslayer;
  FloatVec m_Iter1_Expsh_x, m_Iter1_Expsh_y, m_Iter1_Expsh_z, m_Iter1_Expsh_E; 
  FloatVec m_Iter1_shower0E, m_Iter1_shower0X, m_Iter1_shower0Y, m_Iter1_shower0Z, m_Iter1_shower0T1, m_Iter1_shower0T2;
  FloatVec m_Iter1_shower1E, m_Iter1_shower1X, m_Iter1_shower1Y, m_Iter1_shower1Z, m_Iter1_shower1T1, m_Iter1_shower1T2;

  int      m_Iter1_N2dshower;
  IntVec   m_Iter1_2ds_dlayer, m_Iter1_2ds_part, m_Iter1_2ds_stave;
  FloatVec m_Iter1_2ds_x, m_Iter1_2ds_y, m_Iter1_2ds_z, m_Iter1_2ds_E;

  int      m_Iter1_Ngoodclus, m_Iter1_Nbadclus;
  IntVec   m_Iter1_clus_Nly;
  FloatVec m_Iter1_clus_x, m_Iter1_clus_y, m_Iter1_clus_z, m_Iter1_clus_E;
  FloatVec m_Iter1_clus_px, m_Iter1_clus_py, m_Iter1_clus_pz;
  FloatVec m_Iter1_gclus_2dshx, m_Iter1_gclus_2dshy, m_Iter1_gclus_2dshz, m_Iter1_gclus_2dshE;
  FloatVec m_Iter1_bclus_2dshx, m_Iter1_bclus_2dshy, m_Iter1_bclus_2dshz, m_Iter1_bclus_2dshE;

  TTree* t_dataColIter2;
  FloatVec m_Iter2_ly_dlayer, m_Iter2_ly_part, m_Iter2_ly_stave, m_Iter2_ly_module, m_Iter2_ly_x, m_Iter2_ly_y, m_Iter2_ly_z;
  IntVec   m_Iter2_ly_NshowerX, m_Iter2_ly_NshowerY, m_Iter2_ly_NclusterX, m_Iter2_ly_NclusterY, m_Iter2_Nexpsh;
  FloatVec m_Iter2_barx, m_Iter2_bary, m_Iter2_barz, m_Iter2_barE;
  IntVec   m_Iter2_barslayer;
  FloatVec m_Iter2_Expsh_x, m_Iter2_Expsh_y, m_Iter2_Expsh_z, m_Iter2_Expsh_E;
  FloatVec m_Iter2_shower0E, m_Iter2_shower0X, m_Iter2_shower0Y, m_Iter2_shower0Z, m_Iter2_shower0T1, m_Iter2_shower0T2;
  FloatVec m_Iter2_shower1E, m_Iter2_shower1X, m_Iter2_shower1Y, m_Iter2_shower1Z, m_Iter2_shower1T1, m_Iter2_shower1T2;

  int      m_Iter2_N2dshower;
  IntVec   m_Iter2_2ds_dlayer, m_Iter2_2ds_part, m_Iter2_2ds_stave;
  FloatVec m_Iter2_2ds_x, m_Iter2_2ds_y, m_Iter2_2ds_z, m_Iter2_2ds_E;

  int      m_Iter2_Ngoodclus, m_Iter2_Nbadclus;
  IntVec   m_Iter2_clus_Nly;
  FloatVec m_Iter2_clus_x, m_Iter2_clus_y, m_Iter2_clus_z, m_Iter2_clus_E;
  FloatVec m_Iter2_clus_px, m_Iter2_clus_py, m_Iter2_clus_pz;
  FloatVec m_Iter2_gclus_2dshx, m_Iter2_gclus_2dshy, m_Iter2_gclus_2dshz, m_Iter2_gclus_2dshE;
  FloatVec m_Iter2_bclus_2dshx, m_Iter2_bclus_2dshy, m_Iter2_bclus_2dshz, m_Iter2_bclus_2dshE;


  //Stage2: check reconstructed result
  TTree *t_recoPFO;
  int m_Npfo, m_Nmc, m_N3dclus; 
  FloatVec m_recPFO_px, m_recPFO_py, m_recPFO_pz, m_recPFO_En;
  FloatVec m_2DShower_x, m_2DShower_y, m_2DShower_z, m_2DShower_E; 
  IntVec m_N2dshInClus;
  IntVec m_mcPdgid, m_mcNdaughter, m_mcNparent, m_mcStatus; 
  FloatVec m_mcPx, m_mcPy, m_mcPz, m_mcEn;

  void ClearBar();
  void ClearRecPFO(); 
  void ClearIter0();
  void ClearIter1();
  void ClearIter2();

};
#endif
