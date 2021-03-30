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



  //Parameters for algorithm settings.
  //TODO: use a xml file to readin all parameters. 
  mutable Gaudi::Property<std::string> _filename{this, "OutFileName", "testout.root", "Output file name"};
  mutable Gaudi::Property<float> _seed{this,   "Seed", 2131, "Random Seed"};
  mutable Gaudi::Property<int>  _Debug{this,   "Debug", 0, "Debug level"};
  mutable Gaudi::Property<int>  _Nskip{this,   "SkipEvt", 0, "Skip event"};


  // Output collections
  // output: PFOs 

  //PFA input end here. 
  //-----------------------------------------------------------
  //Followings are for code testing. 
  typedef std::vector<float> FloatVec;
  TFile* m_wfile;
  //Stage0: check input simbars. 
  TTree* t_SimBar;
  FloatVec m_simBar_x, m_simBar_y, m_simBar_z, m_simBar_T1, m_simBar_T2, m_simBar_Q1, m_simBar_Q2, m_simBar_dlayer, m_simBar_part, m_simBar_block, m_simBar_slayer, m_simBar_module;
  //Stage1: EnergySplittingAlg result: clustered bars. 
  TTree *t_PreRec;
  FloatVec m_PreRec_Bar0x, m_PreRec_Bar0y, m_PreRec_Bar0z, m_PreRec_Bar0E, m_PreRec_Bar1x, m_PreRec_Bar1y, m_PreRec_Bar1z, m_PreRec_Bar1E;
  FloatVec m_PreRec_shower0E, m_PreRec_shower0X, m_PreRec_shower0Y, m_PreRec_shower0Z, m_PreRec_shower0T1, m_PreRec_shower0T2;
  FloatVec m_PreRec_shower1E, m_PreRec_shower1X, m_PreRec_shower1Y, m_PreRec_shower1Z, m_PreRec_shower1T1, m_PreRec_shower1T2;
  int m_PreRec_NshowerX, m_PreRec_NshowerY, m_PreRec_NclusterX, m_PreRec_NclusterY;

  TTree *t_recoPFO;
  FloatVec m_recPFO_px, m_recPFO_py, m_recPFO_pz, m_recPFO_En;
  FloatVec m_2DShower_x, m_2DShower_y, m_2DShower_z, m_2DShower_E; 
  int m_Npfo, m_N2dshInClus; 

  void ClearBar();
  void ClearPreRec();
  void ClearRecPFO(); 

};
#endif
