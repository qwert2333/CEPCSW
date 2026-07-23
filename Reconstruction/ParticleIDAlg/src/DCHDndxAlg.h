#ifndef DCHDndxAlg_h
#define DCHDndxAlg_h 1

#include "k4FWCore/DataHandle.h"

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

#include "GaudiKernel/Algorithm.h"

#include "SimplePIDSvc/ISimplePIDSvc.h"
#include "DetInterface/IGeomSvc.h"

/**
 * @class DCHDndxAlg
 * @brief Algorithm to calculate dN/dx for DCH
 * @author Guang Zhao (zhaog@ihep.ac.cn)
*/

class DCHDndxAlg : public Algorithm {
public:
    DCHDndxAlg(const std::string& name, ISvcLocator* svcLoc);

    virtual StatusCode initialize();
    virtual StatusCode execute();
    virtual StatusCode finalize();

protected:
    Gaudi::Property<std::string> m_method{this, "Method", "Simple"};

    DataHandle<edm4hep::TrackCollection> _trackCol{"SDTRecTrackCollection", Gaudi::DataHandle::Reader, this};
    //DataHandle<edm4hep::MCRecoTrackerAssociationCollection> _hitassCol{"DCHitAssociationCollection", Gaudi::DataHandle::Reader, this};
    DataHandle<CEPCSWTrackMCParticleLinkCollection> _trkParAssCol{"SDTRecTrackCollectionParticleAssociation", Gaudi::DataHandle::Reader, this};
    DataHandle<edm4hep::RecDqdxCollection> _dndxCol{"DndxTracks", Gaudi::DataHandle::Writer, this};

private:
    SmartIF<ISimplePIDSvc> m_pid_svc;
    SmartIF<IGeomSvc> m_geom_svc;
    void getFirstAndLastHitsByZ(const podio::RelationRange<edm4hep::TrackerHit>& hitcol, int& first, int& last);
    void getFirstAndLastHitsByRadius(const podio::RelationRange<edm4hep::TrackerHit>& hitcol, int& first, int& last);
    int getDetTypeID(unsigned long long cellID) const;
    bool isCDCHit(edm4hep::TrackerHit* hit);
};

#endif
