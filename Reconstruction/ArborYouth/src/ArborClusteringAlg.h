#ifndef ARBORYOUTHALG_H
#define ARBORYOUTHALG_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/Segmentations.h> 
#include "DetInterface/IGeomSvc.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TString.h"
#include "TH3.h"
#include "TH1.h"


class ArborClusteringAlg : public GaudiAlgorithm
{
 
public:
 
  ArborClusteringAlg(const std::string& name, ISvcLocator* svcLoc);
 
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


	void Clear();

protected:

  SmartIF<IGeomSvc> m_geosvc;
  //SmartIF<ICRDEcalSvc> m_edmsvc;

	int _nEvt ;
	TRandom3 rndm;
  int _seed = 1234;	

	dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
	dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

  //Input collections
  Gaudi::Property< std::string > name_MCParticleCol{ this, "MCParticleCollection", "MCParticle" };
  Gaudi::Property< std::string > name_TrackCol{ this, "TrackCollections", "MarlinTrkTracks" };
  Gaudi::Property< std::string > name_HcalHits{ this, "HCalCaloHitCollections", "HCALBarrel" };
  Gaudi::Property< std::string > name_HcalReadout{ this, "HCalReadOutNames", "HcalBarrelCollection" };

  DataHandle<edm4hep::MCParticleCollection>* r_MCParticleCol;
  DataHandle<edm4hep::TrackCollection>* r_TrackCols;
  DataHandle<edm4hep::CalorimeterHitCollection>* r_HCalHitCols;

  //Input parameters
  mutable Gaudi::Property<int>   _Nskip{this,  "SkipEvt", 0, "Skip event"};

  //Output collections
  DataHandle<edm4hep::ClusterCollection>                w_ClusterCollection {"ArborClusters",Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::ReconstructedParticleCollection>  w_ReconstructedParticleCollection {"ArborYouthPFOs"    ,Gaudi::DataHandle::Writer, this};

  //For Ana: Output tuples
  Gaudi::Property<std::string> m_filename{this, "AnaFileName", "testout.root", "Output file name"};
  typedef std::vector<float> FloatVec;
  typedef std::vector<int>   IntVec;

  TFile* m_wfile;
  TTree* t_digiHit;
  FloatVec m_digiHit_x, m_digiHit_y, m_digiHit_z, m_digiHit_E;
  FloatVec m_digiHit_layer;

  TTree * t_Track;
  int m_Ntrk;
  FloatVec m_trkstate_d0, m_trkstate_z0, m_trkstate_phi, m_trkstate_tanL, m_trkstate_omega, m_trkstate_kappa;
  FloatVec m_trkstate_refx, m_trkstate_refy, m_trkstate_refz;
  IntVec m_trkstate_location;

  TTree *t_Cluster;
  int m_Nclus;
  FloatVec m_Clus_x, m_Clus_y, m_Clus_z, m_Clus_E;
  IntVec m_Clus_Nhit;

  void ClearHit();
  void ClearTrack();
  void ClearCluster();  


};

#endif
