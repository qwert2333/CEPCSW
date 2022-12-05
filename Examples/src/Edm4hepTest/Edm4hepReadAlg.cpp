#include "Edm4hepReadAlg.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"

DECLARE_COMPONENT(Edm4hepReadAlg)

Edm4hepReadAlg::Edm4hepReadAlg(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc)
{
    declareProperty("HeaderCol", m_headerCol);
    declareProperty("MCParticleCol", m_mcParCol, "MCParticle collection (input)");
    declareProperty("SimCalorimeterHitCol", m_calorimeterCol, "MCParticle collection (input)");
}

StatusCode Edm4hepReadAlg::initialize()
{
    debug() << "begin initialize Edm4hepReadAlg" << endmsg;

    m_wfile = new TFile("Edm4hepRead.root", "recreate");
    m_wtree = new TTree("MCInfo", "MCInfo");
    m_wtree->Branch("Nmc", &m_Nmc);
    m_wtree->Branch("deltaTheta_yy", &m_deltaTheta_yy);
    m_wtree->Branch("mcPdgid",     &m_mcPdgid);
    m_wtree->Branch("mcStatus",    &m_mcStatus);
    m_wtree->Branch("mcPx", &m_mcPx);
    m_wtree->Branch("mcPy", &m_mcPy);
    m_wtree->Branch("mcPz", &m_mcPz);
    m_wtree->Branch("mcEn", &m_mcEn);    

    return GaudiAlgorithm::initialize();
}

StatusCode Edm4hepReadAlg::execute()
{
    debug() << "begin execute Edm4hepReadAlg" << endmsg;

    Clear();
    auto mcCol = m_mcParCol.get();
    m_Nmc = mcCol->size();
    std::vector<TVector3> p_gam; p_gam.clear();
    for ( auto p : *mcCol ) {
        //info() << p.getObjectID().index << " : [";
        //for ( auto it = p.daughters_begin(), end = p.daughters_end(); it != end; ++it ) {
        //    info() << " " << it->getObjectID().index;
        //}
        //info() << " ]; ";
        m_mcPdgid.push_back( p.getPDG() );
        m_mcStatus.push_back( p.getGeneratorStatus() );
        m_mcPx.push_back( p.getMomentum()[0] );
        m_mcPy.push_back( p.getMomentum()[1] );
        m_mcPz.push_back( p.getMomentum()[2] );
        m_mcEn.push_back( p.getEnergy() );        

        if(p.getPDG()==22){
          TVector3 gam(p.getMomentum()[0], p.getMomentum()[1], p.getMomentum()[2]);
          p_gam.push_back(gam);
        }
    }
    if(p_gam.size()==2) m_deltaTheta_yy=p_gam[0].Angle(p_gam[1]);
    m_wtree->Fill();

    //info() << "}" << endmsg;

/*    auto caloCol = m_calorimeterCol.get();
    for (auto calohit : *caloCol) {
        unsigned int contrib_size = calohit.contributions_size();
        info() << " contributions_size: " 
               << contrib_size
               << endmsg;
        for (unsigned int i = 0; i < contrib_size; ++i) {
            auto contrib = calohit.getContributions(i);
            auto primary_particle = contrib.getParticle();

            info() << " - #" << i << ": "
                   << " track with "
                   << " PDG: " << contrib.getPDG() // current track
                   << ". "
                   << " primary track with "
                   << " PDG: " << primary_particle.getPDG()
                   << endmsg;


        }

    }
*/

    return StatusCode::SUCCESS;
}

StatusCode Edm4hepReadAlg::finalize()
{
    m_wfile->cd();
    m_wtree->Write();
    m_wfile->Close();
    delete m_wfile, m_wtree;
    debug() << "begin finalize Edm4hepReadAlg" << endmsg;
    return GaudiAlgorithm::finalize();
}

void Edm4hepReadAlg::Clear(){
  m_Nmc=-99;
  m_deltaTheta_yy = -99;
  m_mcPdgid.clear();
  m_mcStatus.clear();
  m_mcPx.clear();
  m_mcPy.clear();
  m_mcPz.clear();
  m_mcEn.clear();
}
