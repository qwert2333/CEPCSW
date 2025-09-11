#include "FinalPIDAlg.h"

#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "DetInterface/IGeomSvc.h"

using namespace edm4hep;

DECLARE_COMPONENT( FinalPIDAlg )

//------------------------------------------------------------------------------
FinalPIDAlg::FinalPIDAlg( const std::string& name, ISvcLocator* pSvcLocator )
    : Algorithm( name, pSvcLocator ) {
  // input
  declareProperty("PFOCollection", m_inPFOCol, "Reconstructed particles to be PIDed");
  declareProperty("TOFPIDCollection", m_inTofCol, "TOF PID information");
  declareProperty("TPCPIDCollection", m_inDqdxCol, "TPC PID information");
  declareProperty("BarrelTrackerHits", m_inputMuonBarrel, "handler of the input barrel tracker hit collection");
  declareProperty("EndcapTrackerHits", m_inputMuonEndcap, "handler of the input endcap tracker hit collection");

  // output
  declareProperty("OutputPFOName", m_outPFOCol, "Reconstructed particles with PID information");
  declareProperty("ParticleIDName", m_ParticleID, "Particle ID collection");
}

//------------------------------------------------------------------------------
StatusCode FinalPIDAlg::initialize(){
  debug() << "Initializing..." << endmsg;
  _nEvt = 0;
  m_pid_svc = service("FinalPIDSvc");
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode FinalPIDAlg::execute(){
  
  const edm4hep::ReconstructedParticleCollection* inpfocol = nullptr;
  const edm4hep::RecTofCollection* tofcol = nullptr;
  const edm4hep::RecDqdxCollection* dqdxcol = nullptr;
  const edm4hep::TrackerHitCollection* barrelhits = nullptr;
  const edm4hep::TrackerHitCollection* endcaphits = nullptr;

  _hasTOF = true;
  _hasTPC = true;
  _hasMuonBarrel = true;
  _hasMuonEndcap = true;

  //auto outpfocol = m_outPFOCol.createAndPut();
  edm4hep::ReconstructedParticleCollection* outpfocol = m_outPFOCol.createAndPut();
  edm4hep::ParticleIDCollection* ParticleID = m_ParticleID.createAndPut();

  try {
    inpfocol = m_inPFOCol.get();
  }
  catch ( GaudiException &e ) {
    debug() << " Reconstructed particle " << m_inPFOCol.fullKey() << " is unavailable in event " << _nEvt << endmsg;
    _nEvt++;
    return StatusCode::SUCCESS;
  }
  try {
    tofcol = m_inTofCol.get();
  }
  catch ( GaudiException &e ) {
    _hasTOF = false;
  }
  try {
    dqdxcol = m_inDqdxCol.get();
  }
  catch ( GaudiException &e ) {
    _hasTPC = false;
  }
  try{
    barrelhits = m_inputMuonBarrel.get();
  }
  catch ( GaudiException &e ) {
    _hasMuonBarrel = false;
  }

  try{
    endcaphits = m_inputMuonEndcap.get();
  }
  catch ( GaudiException &e ) {
    _hasMuonEndcap = false;
  }

  debug()<<" has TOF : "<<_hasTOF<<" has TPC : "<<_hasTPC<<" TOF size : "<<tofcol->size()<<" dqdx size : "<<dqdxcol->size()<<endmsg;
  debug()<<" has Muon Barrel : "<<_hasMuonBarrel<<" has Muon Endcap : "<<_hasMuonEndcap<<" Muon Barrel size : "<<barrelhits->size()<<" Muon Endcap size : "<<endcaphits->size()<<endmsg;

  m_pid_svc->SetCollections(barrelhits, endcaphits, tofcol, dqdxcol, inpfocol);
  m_pid_svc->MatchMuonHitsToTracks();
  m_pid_svc->SetWP_mu(WP::Best);
  m_pid_svc->SetWP_ele(WP::Best);

  debug()<<"Begin loop over PFOs"<<endmsg;

  for ( auto pfo : *inpfocol ){
    auto outpfo = pfo.clone();
    bool load=m_pid_svc->LoadPFO(pfo);
    if (!load) continue;

    if (outpfo.getCharge()==0) {
      for (unsigned int i_fl=0;i_fl<2;i_fl++) {
        int PDG=PDGIDs.at(i_fl+5);
        edm4hep::MutableParticleID pid(PDG, PDG, 0, m_pid_svc->GetProb(i_fl));
        ParticleID->push_back(pid);
        outpfo.addToParticleIDs(pid);
      }
      int PDG=m_pid_svc->GetType();
      outpfo.setType( PDG );
      outpfo.setMass( ParticleMass.at( abs(PDG) ) );
      outpfo.setEnergy( sqrt(outpfo.getMomentum()[0]*outpfo.getMomentum()[0] + outpfo.getMomentum()[1]*outpfo.getMomentum()[1] + outpfo.getMomentum()[2]*outpfo.getMomentum()[2] + outpfo.getMass()*outpfo.getMass()) );
    }
    else {
      for (unsigned int i_fl=0;i_fl<5;i_fl++) {
        int PDG=PDGIDs.at(i_fl)*outpfo.getCharge();
        edm4hep::MutableParticleID pid(PDG, PDG, 0, m_pid_svc->GetProb(i_fl));
        ParticleID->push_back(pid);
        outpfo.addToParticleIDs(pid);
      }
      int PDG=m_pid_svc->GetType();
      if (PDG==11 || PDG==13) PDG*=-1;
      PDG*=outpfo.getCharge();
      outpfo.setType( PDG );
      outpfo.setMass( ParticleMass.at( abs(PDG) ) );
      outpfo.setEnergy( sqrt(outpfo.getMomentum()[0]*outpfo.getMomentum()[0] + outpfo.getMomentum()[1]*outpfo.getMomentum()[1] + outpfo.getMomentum()[2]*outpfo.getMomentum()[2] + outpfo.getMass()*outpfo.getMass()) );
    }

    outpfocol->push_back(outpfo);

    debug()<<"Pdgid: "<<outpfo.getType()<<endmsg;
  }

  _nEvt++;

  return StatusCode::SUCCESS;
}// end execute

//------------------------------------------------------------------------------
StatusCode FinalPIDAlg::finalize(){
  debug() << "Finalizing..." << endmsg;
  return StatusCode::SUCCESS;
}