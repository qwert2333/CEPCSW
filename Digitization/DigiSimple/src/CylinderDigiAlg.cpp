#include "CylinderDigiAlg.h"

#include "DetIdentifier/CEPCConf.h"

#include "edm4hep/Vector3f.h"

#include "DD4hep/Detector.h"
#include <DD4hep/Objects.h>
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/Vector3D.h"

#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IRndmGen.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

DECLARE_COMPONENT( CylinderDigiAlg )

CylinderDigiAlg::CylinderDigiAlg(const std::string& name, ISvcLocator* svcLoc)
: Algorithm(name, svcLoc){
  // Input collections
  declareProperty("SimTrackHitCollection", m_inputColHdls, "Handle of the Input SimTrackerHit collection");

  // Output collections
  declareProperty("TrackerHitCollection", m_outputColHdls, "Handle of the output TrackerHit collection");
  declareProperty("TrackerHitAssociationCollection", m_assColHdls, "Handle of the Association collection between SimTrackerHit and TrackerHit");
}

StatusCode CylinderDigiAlg::initialize(){
  m_geosvc = service<IGeomSvc>("GeomSvc");
  if(!m_geosvc){
    error() << "Failed to get the GeomSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  auto detector = m_geosvc->lcdd();
  if(!detector){
    error() << "Failed to get the Detector from GeomSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  std::string name = m_inputColHdls.objKey(); 
  debug() << "Readout name: " << name << endmsg;
  m_decoder = m_geosvc->getDecoder(name);
  if(!m_decoder){
    error() << "Failed to get the decoder. " << endmsg;
    return StatusCode::FAILURE;
  }

  info() << "CylinderDigiAlg::initialized" << endmsg;
  return Algorithm::initialize();
}


StatusCode CylinderDigiAlg::execute(){
  auto trkhitVec = m_outputColHdls.createAndPut();
  auto assVec = m_assColHdls.createAndPut();

  const edm4hep::SimTrackerHitCollection* STHCol = nullptr;
  try {
    STHCol = m_inputColHdls.get();
  }
  catch(GaudiException &e){
    debug() << "Collection " << m_inputColHdls.fullKey() << " is unavailable in event " << m_nEvt << endmsg;
    return StatusCode::SUCCESS;
  }
  if(STHCol->size()==0) return StatusCode::SUCCESS;
  debug() << m_inputColHdls.fullKey() << " has SimTrackerHit "<< STHCol->size() << endmsg;
  
  for(auto simhit : *STHCol){
#if edm4hep_VERSION >= EDM4HEP_VERSION(1, 0, 0)
    auto particle = simhit.getParticle();
#else
    auto particle = simhit.getMCParticle();
#endif
    if(!particle.isAvailable()) continue;
    
    auto& mom0 = particle.getMomentum();
    double pt = sqrt(mom0.x*mom0.x+mom0.y*mom0.y);
    if(particle.isCreatedInSimulation()&&pt<0.01&&particle.isStopped()) continue;

    auto cellId = simhit.getCellID();
    int system  = m_decoder->get(cellId, "system");
    int layer   = m_decoder->get(cellId, "layer"  );
    int module  = m_decoder->get(cellId, "module");
    int sensor  = m_decoder->get(cellId, "sensor"  );
    auto& pos   = simhit.getPosition();
    auto& mom   = simhit.getMomentum();
    
    double phi = atan2(pos[1], pos[0]);
    double r   = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
    double dphi = m_resRPhi/r;
    phi += randSvc()->generator(Rndm::Gauss(0, dphi))->shoot();
    double smearedX = r*cos(phi);
    double smearedY = r*sin(phi);
    double smearedZ = pos[2] + randSvc()->generator(Rndm::Gauss(0, m_resZ))->shoot();

    auto trkHit = trkhitVec->create();
    trkHit.setCellID(cellId);
    trkHit.setTime(simhit.getTime());
    trkHit.setEDep(simhit.getEDep());
    trkHit.setPosition (edm4hep::Vector3d(smearedX, smearedY, smearedZ));
    trkHit.setCovMatrix(std::array<float, 6>{m_resRPhi*m_resRPhi/2, 0, m_resRPhi*m_resRPhi/2, 0, 0, m_resZ*m_resZ});
    trkHit.setType(1<<CEPCConf::TrkHitTypeBit::CYLINDER);
#if edm4hep_VERSION >= EDM4HEP_VERSION(1, 0, 0)
    //throw std::runtime_error("The RAWHIT related interfaces are removed from TrackerHit");
#else
    trkHit.addToRawHits(simhit.getObjectID());
#endif
    debug() << "Hit " << simhit.id() << ": " << pos << " -> " << trkHit.getPosition() << "s:" << system << " l:" << layer << " m:" << module << " s:" << sensor
	    << " pt = " << pt << " " << mom.x << " " << mom.y << " " << mom.z << endmsg;

    auto ass = assVec->create();

    float weight = 1.0;

    debug() <<" Set relation between " << " sim hit " << simhit.id() << " to tracker hit " << trkHit.id() << " with a weight of " << weight << endmsg;
#if edm4hep_VERSION >= EDM4HEP_VERSION(1, 0, 0)
    ass.setTo(simhit);
    ass.setFrom(trkHit);
#else
    ass.setSim(simhit);
    ass.setRec(trkHit);
#endif
    ass.setWeight(weight);
  }

  m_nEvt++;

  return StatusCode::SUCCESS;
}

StatusCode CylinderDigiAlg::finalize(){
  info() << "Processed " << m_nEvt << " events " << endmsg;
  return Algorithm::finalize();
}
