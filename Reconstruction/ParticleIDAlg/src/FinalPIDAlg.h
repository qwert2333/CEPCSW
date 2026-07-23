#ifndef FinalPIDAlg_h
#define FinalPIDAlg_h 1

#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/Algorithm.h"

#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4cepc/RecTof.h"
#include "edm4cepc/RecTofCollection.h"
#include "edm4hep/RecDqdx.h"
#include "edm4hep/RecDqdxCollection.h"
#include "edm4hep/ParticleIDCollection.h"
#include "edm4hep/EDM4hepVersion.h"
#if edm4hep_VERSION >= EDM4HEP_VERSION(1, 0, 0)
#include "edm4hep/TrackerHit3DCollection.h"
using CEPCSWTrackerHit3DCollection = edm4hep::TrackerHit3DCollection;
#else
#include "edm4hep/TrackerHitCollection.h"
using CEPCSWTrackerHit3DCollection = edm4hep::TrackerHitCollection;
#endif
#include "TVector3.h"

#include "FinalPIDSvc/IFinalPIDSvc.h"
#include "DetInterface/IGeomSvc.h"
#include <GaudiKernel/Service.h>

#include "FinalPIDSvc/WorkingPoint.h"

class FinalPIDAlg : public Algorithm {
 public:
  // Constructor of this form must be provided
  FinalPIDAlg( const std::string& name, ISvcLocator* pSvcLocator );

  // Three mandatory member functions of any algorithm
  StatusCode initialize() override;
  StatusCode execute() override;
  StatusCode finalize() override;

 private:
  DataHandle<edm4hep::ReconstructedParticleCollection> m_inPFOCol{"CyberPFO", Gaudi::DataHandle::Reader, this};

  DataHandle<edm4hep::RecTofCollection> m_inTofCol{"RecTofCollection", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::RecDqdxCollection> m_inDqdxCol{"DndxTracks", Gaudi::DataHandle::Reader, this};
  DataHandle<CEPCSWTrackerHit3DCollection> m_inputMuonBarrel{"MuonBarrelTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<CEPCSWTrackerHit3DCollection> m_inputMuonEndcap{"MuonEndcapTrackerHits", Gaudi::DataHandle::Reader, this};
  //DataHandle<edm4hep::ParticleIDCollection> m_PIDCol{"finalPID", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::ReconstructedParticleCollection> m_outPFOCol{"CyberPFOPID", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::ParticleIDCollection> m_ParticleID{"ParticleID", Gaudi::DataHandle::Writer, this};

  SmartIF<IFinalPIDSvc> m_pid_svc;

  // Gaudi::Property<std::string> m_method{this, "PIDMethod", "TPC+TOF+CALO"};


  // void FillTPCPID(const edm4hep::RecDqdxCollection* dqdxcol, edm4hep::MutableReconstructedParticle& pfo, std::array<double, 5>& chi2s);
  // void FillTOFPID(const edm4hep::RecTofCollection* tofcol, edm4hep::MutableReconstructedParticle& pfo, std::array<double, 5>& chi2s);
  // StatusCode FillCaloPID(edm4hep::MutableReconstructedParticle& pfo);

  int _nEvt;
  bool _hasTPC;
  bool _hasTOF;
  bool _hasMuonBarrel;
  bool _hasMuonEndcap;

  //Detector geometry size
  // const float EcalOuterR = 2130.;
  // const float EcalHalfZ = 3230.;
  // const float HcalOuterR = 3455.;
  // const float HcalHalfZ = 4575.;


  const std::map<int, int> PDGIDs = {
    {0, -11},
    {1, -13},
    {2, 211},
    {3, 321},
    {4, 2212},
    {5, 22},
    {6, 130}
  };
  //Particle mass from PDGLive 2024 edition [https://pdglive.lbl.gov/Viewer.action]
  const std::map<int, double> ParticleMass = {
    {11, 0.000511},
    {13, 0.105658},
    {211, 0.139570},
    {321, 0.493677},
    {2212, 0.938272},
    {22, 0.},
    //{130, 0.497611},
    {130, 0.} //Temp: set K_0 mass as 0, due to imperfect neutral hadron reconstruction in CyberPFA. 
  };

};
#endif
