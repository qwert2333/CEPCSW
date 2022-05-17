#ifndef PANDORAPLUS_ALG_H
#define PANDORAPLUS_ALG_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "CRDEcalSvc/ICRDEcalSvc.h"
#include "DetInterface/IGeomSvc.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackCollection.h"


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

  //Creators and their setting

  //Algorithm for PFA


  //Parameters for PFA algorithm



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

};
#endif
