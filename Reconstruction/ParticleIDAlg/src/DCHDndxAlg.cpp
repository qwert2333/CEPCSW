#include "DCHDndxAlg.h"
#include "DataHelper/HelixClass.h"
#include "DataHelper/Navigation.h"
#include "DetInterface/IGeomSvc.h"
#include "UTIL/ILDConf.h"

#include "DD4hep/Detector.h"
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
#include <fstream>

using namespace edm4hep;

DECLARE_COMPONENT(DCHDndxAlg)

DCHDndxAlg::DCHDndxAlg(const std::string& name, ISvcLocator* svcLoc) : Algorithm(name, svcLoc) {
    // Input
    declareProperty("SDTRecTrackCollection", _trackCol, "handler of the input track collection");
    declareProperty("SDTRecTrackCollectionParticleAssociation", _trkParAssCol, "handler of the input track particle association collection");

    // Output
    declareProperty("DndxTracks", _dndxCol, "handler of the collection of dN/dx tracks");
}

StatusCode DCHDndxAlg::initialize() {
    info() << "Initilize DndxAlg ..." << endmsg;

    if (m_method == "Simple") {
        m_pid_svc = service("SimplePIDSvc");
    }
    else {
        m_pid_svc = nullptr;
    }

    m_geom_svc = service("GeomSvc");

    return Algorithm::initialize();
}

StatusCode DCHDndxAlg::execute() {
    info() << "Dndx reconstruction started" << endmsg;
    // static std::ofstream check_file("log.txt");
    // static int ievent = 0;
    // check_file << "--------                 Event " << ievent++ << " --------" << std::endl;

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

    // check_file << "before track loop" << std::endl;

    for (std::size_t i = 0; i < trkcol->size(); i++) {
        // check_file << "  Track " << i << std::endl;
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

        // check_file << "before track par" << std::endl;

        /// Track parameters
        podio::RelationRange<edm4hep::TrackerHit> hitcol = trk.getTrackerHits();
        if (hitcol.size() == 0) return StatusCode::SUCCESS;

        // check_file << "before get track state" << std::endl;
        int ihit_first = -1;
        int ihit_last = -1;
        TrackState trk_par;
        for (auto it = trk.trackStates_begin(); it != trk.trackStates_end(); it++) {
            if (it->location == 0) {
                trk_par = *it;
                break;
            }
        }
        if (trk_par.tanLambda > 0.1) {
            getFirstAndLastHitsByZ(hitcol, ihit_first, ihit_last);
        }
        else {
            // check_file << "before get first and last by r" << std::endl;
            getFirstAndLastHitsByRadius(hitcol, ihit_first, ihit_last);
        }
        if (ihit_first < 0 || ihit_last < 0 || ihit_first == ihit_last) continue;

        // check_file << "before conversion unit. ihit_first, ihit_last = " << ihit_first << ", " << ihit_last << std::endl;
        
        // convert the position from cm to mm
        Vector3d first_hit_pos = Vector3d(hitcol[ihit_first].getPosition().x*10.0, hitcol[ihit_first].getPosition().y*10.0, hitcol[ihit_first].getPosition().z*10.0);
        Vector3d last_hit_pos = Vector3d(hitcol[ihit_last].getPosition().x*10.0, hitcol[ihit_last].getPosition().y*10.0, hitcol[ihit_last].getPosition().z*10.0);
        double len, p, cos;
        m_pid_svc->getTrackPars(first_hit_pos, last_hit_pos, trk, 0, len, p, cos);

        // check_file << "before dndx calc" << std::endl;

        /// dN/dx reconstruction
        if (m_method.value() == "Simple") { // Track level implementation
            int particle_type = m_pid_svc->getParticleType(pdgid);
            double bg = m_pid_svc->getBetaGamma(p_truth, particle_type);
            double dndx_mean = m_pid_svc->getDndxMean(bg, cos);
            double dndx_sigma = m_pid_svc->getDndxSigma(bg, cos, len);
            double dndx_meas = m_pid_svc->getDndx(dndx_mean, dndx_sigma);
            // std::cout << "pdgid: " << pdgid << " bg: " << bg << " dndx_mean: " << dndx_mean << " dndx_sigma: " << dndx_sigma << " dndx_meas: " << dndx_meas << std::endl;

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

StatusCode DCHDndxAlg::finalize() {

    return StatusCode::SUCCESS;
}

void DCHDndxAlg::getFirstAndLastHitsByZ(const podio::RelationRange<edm4hep::TrackerHit>& hitcol, int& first, int& last) {
    bool first_dc_hit = true;
    double zmin, zmax;
    for (size_t ihit = 0; ihit < hitcol.size(); ihit++) {
        auto hit = hitcol[ihit];
        double z = hit.getPosition()[2];
        if (!isCDCHit(&hit)) {
            continue;
        }
        if (first_dc_hit) {
            zmin = z;
            zmax = z;
            first = ihit;
            last = ihit;
            first_dc_hit = false;
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

void DCHDndxAlg::getFirstAndLastHitsByRadius(const podio::RelationRange<edm4hep::TrackerHit>& hitcol, int& first, int& last) {
    bool first_dc_hit = true;
    double rmin, rmax;
    for (size_t ihit = 0; ihit < hitcol.size(); ihit++) {
        auto hit = hitcol[ihit];
        double x = hit.getPosition()[0];
        double y = hit.getPosition()[1];
        double r = sqrt(x*x + y*y);
        if (!isCDCHit(&hit)) {
            continue;
        }
        if (first_dc_hit) {
            rmin = r;
            rmax = r;
            first = ihit;
            last = ihit;
            first_dc_hit = false;
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

int DCHDndxAlg::getDetTypeID(unsigned long long cellID) const {
    UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string);
    encoder.setValue(cellID);
    return encoder[lcio::ILDCellID0::subdet];
}

bool DCHDndxAlg::isCDCHit(edm4hep::TrackerHit* hit){
    return m_geom_svc->lcdd()->constant<int>("DetID_DC")==
        getDetTypeID(hit->getCellID());
}