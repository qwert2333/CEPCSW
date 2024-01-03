#ifndef TEST_EDM4HEP_WRITE_ALG_H
#define TEST_EDM4HEP_WRITE_ALG_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/TrackCollection.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

namespace edm4hep {
    class EventHeaderCollection;
    class MCParticleCollection;
    class SimCalorimeterHitCollection;
    class CaloHitContributionCollection;
}

class Edm4hepReadAlg : public GaudiAlgorithm
{

    public :

        Edm4hepReadAlg(const std::string& name, ISvcLocator* svcLoc);

        virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();

    private :

        DataHandle<edm4hep::EventHeaderCollection> m_headerCol{"EventHeader", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCParticleCollection> m_mcParCol{"MCParticle", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimCalorimeterHitCollection> m_calorimeterCol{"SimCalorimeterCol", 
                Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::CaloHitContributionCollection> m_caloContribCol{"SimCaloContributionCol", 
                Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackCollection> m_trkCol{"MarlinTrkTracks", Gaudi::DataHandle::Reader, this};

    TFile *m_wfile;
    TTree *m_wtree;
    int m_Nmc;
    float m_deltaTheta_yy;
    std::vector<int> m_mcPdgid, m_mcStatus;
    std::vector<float> m_mcPx, m_mcPy, m_mcPz, m_mcEn;

    int m_Ntrk;
    std::vector<int> m_trk_type, m_trk_Nhit;
    std::vector<float> m_trk_aveOmega;

    void Clear();

};

#endif  // TEST_EDM4HEP_WRITE_ALG_H
