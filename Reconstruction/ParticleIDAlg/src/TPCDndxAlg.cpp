#include "TPCDndxAlg.h"
#include "DataHelper/HelixClass.h"
#include "DataHelper/Navigation.h"
#include "DetIdentifier/CEPCConf.h"

#include "edm4hep/EDM4hepVersion.h"
#if edm4hep_VERSION >= EDM4HEP_VERSION(1, 0, 0)
// Hypothesis.h removed in EDM4hep v1.0.0
#else
#include "edm4hep/Hypothesis.h"
#endif
#include "edm4hep/MCParticle.h"
#include "edm4hep/Quantity.h"
#include "edm4hep/Vector3d.h"
#include "lcio.h"

#include <cmath>

using namespace edm4hep;

DECLARE_COMPONENT(TPCDndxAlg)

TPCDndxAlg::TPCDndxAlg(const std::string& name, ISvcLocator* svcLoc) : Algorithm(name, svcLoc) {
    // Input
    declareProperty("CompleteTracks", _trackCol, "handler of the input track collection");
    declareProperty("CompleteTracksParticleAssociation", _trkParAssCol, "handler of the input track particle association collection");

    // Output
    declareProperty("DndxTracks", _dndxCol, "handler of the collection of dN/dx tracks");
}

StatusCode TPCDndxAlg::initialize() {
    info() << "Initilize DndxAlg ..." << endmsg;

    if (m_method == "Simple") {
        m_pid_svc = service("SimplePIDSvc");
    }
    else {
        m_pid_svc = nullptr;
    }

    m_geosvc = service<IGeomSvc>("GeomSvc");
    if (!m_geosvc) {
        error() << "Could not find GeomSvc!" << endmsg;
        return StatusCode::FAILURE;
    }

    m_decoder = m_geosvc->getDecoder("TPCCollection");
    if (!m_decoder) {
        error() << "Could not find TPC decoder!" << endmsg;
        return StatusCode::FAILURE;
    }

    return Algorithm::initialize();
}

StatusCode TPCDndxAlg::execute() {
    info() << "Dndx reconstruction started" << endmsg;

    const edm4hep::TrackCollection* trkcol = nullptr;
    const CEPCSWTrackMCParticleLinkCollection* trkparasscol = nullptr;
    try {
        trkcol = _trackCol.get();
        trkparasscol = _trkParAssCol.get();
    }
    catch(...) {
        //
    }

    if (trkcol == nullptr || trkparasscol == nullptr) {
        return StatusCode::SUCCESS;
    }

    RecDqdxCollection* outCol = _dndxCol.createAndPut();

    // Navigation nav;
    // nav.Initialize();
    for (std::size_t i = 0; i < trkcol->size(); i++) {
        Track trk(trkcol->at(i));
        
        /// MC truth information
        int pdgid = -999;
        double p_truth = -999.;

        float max_weight = -999.;
        int max_weight_idx = -1;
        int ass_idx = 0;
        for (auto ass : *trkparasscol) {
#if edm4hep_VERSION >= EDM4HEP_VERSION(1, 0, 0)
            if (ass.getFrom() == trk) {
#else
            if (ass.getRec() == trk) {
#endif
                float weight = ass.getWeight();
                if (weight > max_weight) {
                    max_weight = weight;
                    max_weight_idx = ass_idx;
                }
            }
            ass_idx++;
        }
        if (max_weight_idx < 0) continue;
#if edm4hep_VERSION >= EDM4HEP_VERSION(1, 0, 0)
        pdgid = trkparasscol->at(max_weight_idx).getTo().getPDG();
        double px = trkparasscol->at(max_weight_idx).getTo().getMomentum()[0];
        double py = trkparasscol->at(max_weight_idx).getTo().getMomentum()[1];
        double pz = trkparasscol->at(max_weight_idx).getTo().getMomentum()[2];
#else
        pdgid = trkparasscol->at(max_weight_idx).getSim().getPDG();
        double px = trkparasscol->at(max_weight_idx).getSim().getMomentum()[0];
        double py = trkparasscol->at(max_weight_idx).getSim().getMomentum()[1];
        double pz = trkparasscol->at(max_weight_idx).getSim().getMomentum()[2];
#endif
        p_truth = sqrt(px*px + py*py + pz*pz);

        /// Track parameters
        podio::RelationRange<edm4hep::TrackerHit> hitcol = trk.getTrackerHits();
        int nhits = hitcol.size();
        if (nhits == 0) return StatusCode::SUCCESS;

        int ihit_first = -1;
        int ihit_last = -1;
        TrackState trk_par;
        for (auto it = trk.trackStates_begin(); it != trk.trackStates_end(); it++) {
            if (it->location == TrackState::AtFirstHit) {
                trk_par = *it;
                break;
            }
        }
        if (trk_par.tanLambda > 0.1) {
            getFirstAndLastHitsByZ(hitcol, ihit_first, ihit_last);
        }
        else {
            getFirstAndLastHitsByRadius(hitcol, ihit_first, ihit_last);
        }
        if (ihit_first < 0 || ihit_last < 0 || ihit_first == ihit_last) continue;
    
        const Vector3d& first_hit_pos = hitcol[ihit_first].getPosition();
        const Vector3d& last_hit_pos = hitcol[ihit_last].getPosition();
        double len, p, cos;
        m_pid_svc->getTrackPars(first_hit_pos, last_hit_pos, trk, TrackState::AtFirstHit, len, p, cos);

        /// dN/dx reconstruction
        if (m_method.value() == "Simple") { // Track level implementation
            int particle_type = m_pid_svc->getParticleType(pdgid);
            double bg = m_pid_svc->getBetaGamma(p_truth, particle_type);
            double dndx_mean = m_pid_svc->getDndxMean(bg, cos);
            double dndx_sigma = m_pid_svc->getDndxSigma(bg, cos, len);
            double dndx_meas = m_pid_svc->getDndx(dndx_mean, dndx_sigma);

            Quantity q;
            q.value = dndx_meas;
            q.error = dndx_sigma;

#if edm4hep_VERSION >= EDM4HEP_VERSION(1, 0, 0)
            throw std::runtime_error("Hypothesis not available in EDM4hep v1.0.0");
#else
            std::array<Hypothesis, 5> hypotheses;
            for (int pid = 0; pid < 5; pid++) {
                bg = m_pid_svc->getBetaGamma(p_truth, pid);
                dndx_mean = m_pid_svc->getDndxMean(bg, cos);
                dndx_sigma = m_pid_svc->getDndxSigma(bg, cos, len);

                Hypothesis h;
                h.chi2 = m_pid_svc->getChi2(dndx_meas, dndx_mean, dndx_sigma);
                h.expected = dndx_mean;
                h.sigma = dndx_sigma;

                hypotheses[pid] = h;
            }

            MutableRecDqdx dndx_track(q, particle_type, 0, hypotheses);
            dndx_track.setTrack(trk);
            outCol->push_back(dndx_track);
#endif
        }
        else if (m_method.value() == "Full") {
            // Hit level implementation, loop over hits ...
        }
        else {
            return StatusCode::FAILURE;
        }
    }


    return StatusCode::SUCCESS;

}

StatusCode TPCDndxAlg::finalize() {

    return StatusCode::SUCCESS;
}

void TPCDndxAlg::getFirstAndLastHitsByRadius(const podio::RelationRange<edm4hep::TrackerHit>& hitcol, int& first, int& last) {
    bool first_hit = true;
    double rmin(-999.), rmax(-999.);
    for (size_t ihit = 0; ihit < hitcol.size(); ihit++) {
        auto hit = hitcol[ihit];
        auto cellid = hit.getCellID();
        auto sysid = m_decoder->get(cellid, "system");
        if (sysid != CEPCConf::DetID::TPC) {
            continue;
        }

        double x = hit.getPosition()[0];
        double y = hit.getPosition()[1];
        double r = sqrt(x*x + y*y);
        if (first_hit) {
            rmin = r;
            rmax = r;
            first = ihit;
            last = ihit;
            first_hit = false;
        }
        else {
            if (r < rmin) {
                rmin = r;
                first = ihit;
            }
            if (r > rmax) {
                rmax = r;
                last = ihit;
            }
        }
    }
}

void TPCDndxAlg::getFirstAndLastHitsByZ(const podio::RelationRange<edm4hep::TrackerHit>& hitcol, int& first, int& last) {
    bool first_hit = true;
    double zmin(-999.), zmax(-999.);
    for (size_t ihit = 0; ihit < hitcol.size(); ihit++) {
        auto hit = hitcol[ihit];
        auto cellid = hit.getCellID();
        auto sysid = m_decoder->get(cellid, "system");
        if (sysid != CEPCConf::DetID::TPC) {
            continue;
        }

        double z = hit.getPosition()[2];
        if (first_hit) {
            zmin = z;
            zmax = z;
            first = ihit;
            last = ihit;
            first_hit = false;
        }
        else {
            if (z < zmin) {
                zmin = z;
                first = ihit;
            }
            if (z > zmax) {
                zmax = z;
                last = ihit;
            }
        }
    }
}
