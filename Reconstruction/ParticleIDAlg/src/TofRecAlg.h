#ifndef TofRecAlg_h
#define TofRecAlg_h 1

#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/Algorithm.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

#include "edm4hep/TrackerHit.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/EDM4hepVersion.h"
#if edm4hep_VERSION >= EDM4HEP_VERSION(1, 0, 0)
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include "edm4hep/TrackMCParticleLinkCollection.h"
using CEPCSWTrackerHit3DCollection = edm4hep::TrackerHit3DCollection;
using CEPCSWTrackerHitSimTrackerHitLinkCollection = edm4hep::TrackerHitSimTrackerHitLinkCollection;
using CEPCSWTrackMCParticleLinkCollection = edm4hep::TrackMCParticleLinkCollection;
#else
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/MCRecoTrackParticleAssociationCollection.h"
using CEPCSWTrackerHit3DCollection = edm4hep::TrackerHitCollection;
using CEPCSWTrackerHitSimTrackerHitLinkCollection = edm4hep::MCRecoTrackerAssociationCollection;
using CEPCSWTrackMCParticleLinkCollection = edm4hep::MCRecoTrackParticleAssociationCollection;
#endif
#include "edm4hep/RecDqdx.h"
#include "edm4hep/RecDqdxCollection.h"
#include "edm4cepc/RecTofCollection.h"
#include "edm4hep/TrackState.h"
#include "SimplePIDSvc/ISimplePIDSvc.h"
#include "edm4hep/Vector3d.h"
#include <random>

class TofRecAlg : public Algorithm {
 public:
  // Constructor of this form must be provided
  TofRecAlg( const std::string& name, ISvcLocator* pSvcLocator );

  // Three mandatory member functions of any algorithm
  StatusCode initialize() override;
  StatusCode execute() override;
  StatusCode finalize() override;

 private:
  DataHandle<edm4hep::TrackCollection> m_FulTrkCol{"CompleteTracks", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::RecTofCollection> m_RecToFCol{"RecTofCollection", Gaudi::DataHandle::Writer, this};
  Gaudi::Property<std::string> m_method{this, "Method", "Simple"};
  SmartIF<ISimplePIDSvc> m_pid_svc;

  void FindToFHits(const edm4hep::Track& _track, bool& _hasFTDHit, bool& _hasSETHit, double& _Toft, double& _Tofx, double& _Tofy, double& _Tofz);
  
  double Tofx, Tofy, Tofz, Toft, Tofterr;
  double position0x, position0y, position0z, position1x, position1y, position1z;
  edm4hep::Vector3d first_hit_pos;
  edm4hep::Vector3d last_hit_pos;

  int _nEvt;

  const std::map<int, double> masses = {//masses in GeV for e mu pi K p
    {0, 0.000511},
    {1, 0.105658},
    {2, 0.139570},
    {3, 0.493677},
    {4, 0.938272},
  };

  float c_mm_ns = 2.99792458e2;//spead of light in mm/ns

  bool hasFTDHit = false;
  bool hasSETHit = false;

  float bunchcrossing = 0.02;//ns
  float tof_resolution = 0.05;//ns
  std::default_random_engine generator_tof;
  std::normal_distribution<float> normal_distribution_tof{0, tof_resolution};
  std::default_random_engine generator_bunch;
  std::normal_distribution<float> normal_distribution_bunch{0, bunchcrossing};
};
#endif