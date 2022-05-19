#ifndef PANDORAPLUS_ALG_H
#define PANDORAPLUS_ALG_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "CRDEcalSvc/ICRDEcalSvc.h"
#include "DetInterface/IGeomSvc.h"
#include "PandoraPlusDataCol.h"

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
  PandoraPlusDataCol     m_DataCol;   //TODO: decide if set it to const. 


  //Creators and their setting
  MCParticleCreator       *m_pMCParticleCreator;
  TrackCreator            *m_pTrackCreator;
  EcalHitsCreator         *m_pEcalHitsCreator;
  //PFOCreator              *m_pPfoCreator;
  

  MCParticleCreator::Settings   m_pMCParticleCreatorSettings;
  TrackCreator::Settings        m_pTrackCreatorSettings;
  EcalHitsCreator::Settings     m_EcalHitsCreatorSettings;


  //Algorithm for PFA


  //Parameters for PFA algorithm



  //Readin collection names
  Gaudi::Property< std::string > name_MCParticleCol{ this, "MCParticleCollection", "MCParticle" };
  Gaudi::Property< std::vector<std::string> > name_TrackCol{ this, "TrackCollections", {"MarlinTrkTracks"} };
  Gaudi::Property< std::vector<std::string> > name_EcalHits{ this, "ECalCaloHitCollections", {"ECALBarrel"} };
  Gaudi::Property< std::vector<std::string> > name_EcalReadout{ this, "ECalReadOutNames", {"EcalBarrelCollection"} }; 


  //Parameters for algorithm settings.
  //TODO: use a xml file to readin all parameters. 
  Gaudi::Property<float> m_seed{this,   "Seed", 2131, "Random Seed"};
  Gaudi::Property<int>   m_Debug{this,   "Debug", 0, "Debug level"};
  Gaudi::Property<int>   m_Nskip{this,   "SkipEvt", 0, "Skip event"};
  Gaudi::Property<bool>  m_WriteAna {this, "WriteAna", false, "Write Ntuples for analysis"};

  // Output collections
  DataHandle<edm4hep::MCParticleCollection>     w_mcParCol  {"MCParticle", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::ClusterCollection>                w_ClusterCollection {"PandoraClusters",Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::ReconstructedParticleCollection>  w_ReconstructedParticleCollection {"PandoraPFOs"    ,Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::VertexCollection>                 w_VertexCollection {"PandoraPFANewStartVertices",Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoParticleAssociationCollection>  w_MCRecoParticleAssociationCollection {"pfoMCRecoParticleAssociation",Gaudi::DataHandle::Writer, this};

  //For Ana
  Gaudi::Property<std::string> _filename{this, "OutFileName", "testout.root", "Output file name"};




  //PFA input end here. 
  //-----------------------------------------------------------

};
#endif
