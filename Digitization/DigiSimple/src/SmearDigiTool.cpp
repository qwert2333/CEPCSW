#include "SmearDigiTool.h"

#include "DataHelper/TrackerHitHelper.h"
#include "DetIdentifier/CEPCConf.h"
#include "DetInterface/IGeomSvc.h"

#include "edm4hep/Vector3f.h"

#include "DD4hep/Detector.h"
#include <DD4hep/Objects.h>
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/ISurface.h"

#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IRndmGen.h"
#include "GaudiKernel/RndmGenerators.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

DECLARE_COMPONENT(SmearDigiTool)

StatusCode SmearDigiTool::initialize() {
  if (m_parameterize) {
    if (m_parU.size()!=10 || m_parV.size()!=10) {
      fatal() << "parameters number must be 10! now " << m_parU.size() << " for U and " << m_parV.size() << " for V" << endmsg;
      return StatusCode::FAILURE;
    }
  }
  else {
    if (m_resU.size() != m_resV.size()) {
      fatal() << "Inconsistent number of resolutions given for U and V coordinate: "
              << "ResolutionU  :" << m_resU.size() << " != ResolutionV : " << m_resV.size()
              << endmsg;
      return StatusCode::FAILURE;
    }
  }

  m_geosvc = service<IGeomSvc>("GeomSvc");
  if(!m_geosvc){
    error() << "Failed to get the GeomSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_detName=="") {
    std::string toolName = name();
    m_detName = toolName.substr(toolName.find(".")+1);
    debug() << toolName << " --> " << m_detName.value() << endmsg;
  }
  m_surfaces = m_geosvc->getSurfaceMap(m_detName.value());
  debug() << "Surface map size of " << m_detName.value() << ": " << m_surfaces->size() << endmsg;
  if (msgLevel(MSG::VERBOSE)) {
    unsigned long old = 0; 
    for (const auto& pair : *m_surfaces) {
      verbose() << pair.first << " : " << pair.second << endmsg;
      if (old==pair.first) error() << old << " repeat!" << endmsg;
      old = pair.first;
    }
  }

  if (m_readoutName=="") m_readoutName = m_detName.value() + "Collection";
  debug() << "Readout name: " << m_readoutName.value() << endmsg;
  m_decoder = m_geosvc->getDecoder(m_readoutName.value());
  if(!m_decoder){
    error() << "Failed to get the decoder. " << endmsg;
    return StatusCode::FAILURE;
  }

  m_randSvc = service<IRndmGenSvc>("RndmGenSvc");

  info() << "initialized" << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode SmearDigiTool::Call(const edm4hep::SimTrackerHitCollection* simCol, IDigiTool::CEPCSWTrackerHit3DCollection* hitCol,
			       IDigiTool::CEPCSWTrackerHitSimTrackerHitLinkCollection* assCol) {
  for (auto simhit : *simCol) {
    StatusCode sc = Call(simhit, hitCol, assCol);
    if (sc.isFailure()) return sc;
  }

  if (m_appendNoise) {
    StatusCode sc = AppendNoise(hitCol);
    if (sc.isFailure()) return sc;
  }

  return StatusCode::SUCCESS;
}

StatusCode SmearDigiTool::Call(edm4hep::SimTrackerHit simhit, IDigiTool::CEPCSWTrackerHit3DCollection* hitCol, IDigiTool::CEPCSWTrackerHitSimTrackerHitLinkCollection* assCol) {
  if (!simhit.isAvailable()) {
    error() << "input SimTrackerHit not available!" << endmsg;
    return StatusCode::SUCCESS;
  }

  auto e = simhit.getEDep();
  if (e < m_eThreshold) return StatusCode::SUCCESS;
  if (m_randSvc->generator(Rndm::Flat(0, 1))->shoot() > m_efficiency) return StatusCode::SUCCESS;
  auto t = simhit.getTime();

  auto cellId = simhit.getCellID();
  int system  = m_decoder->get(cellId, CEPCConf::DetCellID::system);
  int side    = m_decoder->get(cellId, CEPCConf::DetCellID::side);
  int layer   = m_decoder->get(cellId, CEPCConf::DetCellID::layer);
  int module  = m_decoder->get(cellId, CEPCConf::DetCellID::module);
  int sensor  = m_decoder->get(cellId, CEPCConf::DetCellID::sensor);

  auto& pos   = simhit.getPosition();
  auto& mom   = simhit.getMomentum();

  debug() << "Hit " << simhit.id() << " cell: " << cellId << " d:" << system << " l:" << layer << " m:" << module << " s:" << sensor
	  << " p = " << mom.x << ", " << mom.y << ", " << mom.z << " e:" << e << " t:" << t << " isBG: " << simhit.isOverlay() << endmsg;

  dd4hep::rec::ISurface* surface = nullptr;
  auto it = m_surfaces->find(cellId);
  if (it != m_surfaces->end()) {
    surface = it->second;
    if (!surface) {
      fatal() << "found surface for cell id " << cellId << ", but NULL" << endmsg;
      return StatusCode::FAILURE;
    }
  }
  else {
    fatal() << "not found surface for cell id " << cellId << endmsg;
    return StatusCode::FAILURE;
  }

  // CLHEP::mm is divided while Edm4hepWriterAnaElemTool write, so pos without unit
  dd4hep::rec::Vector3D oldPos(pos.x*dd4hep::mm, pos.y*dd4hep::mm, pos.z*dd4hep::mm);
  dd4hep::rec::Vector3D uVec = surface->u(oldPos);
  dd4hep::rec::Vector3D vVec = surface->v(oldPos);
  float u_direction[2];
  u_direction[0] = uVec.theta();
  u_direction[1] = uVec.phi();

  float v_direction[2];
  v_direction[0] = vVec.theta();
  v_direction[1] = vVec.phi();

  debug() << " U: " << uVec << endmsg;
  debug() << " V: " << vVec << endmsg;
  debug() << " N: " << surface->normal() << endmsg;
  debug() << " O: " << surface->origin() << endmsg;

  float resU(0), resV(0), resT(0);
  if (!m_parameterize) {
    if ((m_resU.size() > 1 && layer >= (int)m_resU.size()) || (m_resV.size() > 1 && layer >= (int)m_resV.size())) {
      fatal() << "layer exceeds resolution vector, please check input parameters ResolutionU and ResolutionV" << endmsg;
      return StatusCode::FAILURE;
    }

    resU = (m_resU.size() > 1 ? m_resU.value().at(layer) : m_resU.value().at(0));
    resV = (m_resV.size() > 1 ? m_resV.value().at(layer) : m_resV.value().at(0));
  }
  else { // Riccardo's parameterized model
    CLHEP::Hep3Vector momVec(mom[0], mom[1], mom[2]);
    CLHEP::Hep3Vector uVecCLHEP(uVec[0], uVec[1], uVec[2]);
    CLHEP::Hep3Vector vVecCLHEP(vVec[0], vVec[1], vVec[2]);
    const double alpha = uVecCLHEP.azimAngle(momVec, vVecCLHEP);
    const double cotanAlpha = 1./tan(alpha);
    // TODO: title angle (PI/2), magnetic field (3)
    const double tanLorentzAngle = (side==0) ? 0. : 0.053 * 3 * cos(M_PI/2.);
    const double x = fabs(-cotanAlpha - tanLorentzAngle);
    resU = m_parU[0] + m_parU[1] * x + m_parU[2] * exp(-m_parU[9] * x) * cos(m_parU[3] * x + m_parU[4])
      + m_parU[5] * exp(-0.5 * pow(((x - m_parU[6]) / m_parU[7]), 2)) + m_parU[8] * pow(x, 0.5);

    const double beta = vVecCLHEP.azimAngle(momVec, uVecCLHEP);
    const double cotanBeta = 1./tan(beta);
    const double y = fabs(-cotanBeta);
    resV = m_parV[0] + m_parV[1] * y + m_parV[2] * exp(-m_parV[9] * y) * cos(m_parV[3] * y + m_parV[4])
      + m_parV[5] * exp(-0.5 * pow(((y - m_parV[6]) / m_parV[7]), 2)) + m_parV[8] * pow(y, 0.5);
  }
  resU *= dd4hep::mm;
  resV *= dd4hep::mm;
  // parameterize only for position now, todo
  resT = (m_resT.size() > 1 ? m_resT.value().at(layer) : m_resT.value().at(0));
  // from ps (input unit) to ns (record unit, Geant4)
  resT *= CLHEP::ps/CLHEP::ns;
  debug() << " --- will smear hit with resU = " << resU/dd4hep::mm << " resV = " << resV/dd4hep::mm << " resT = " << resT << endmsg;

  auto& typeSurface = surface->type();
  debug() << " Surface id: " << surface->id() << " type:" << endmsg;
  debug() << " " << typeSurface << endmsg;
  verbose() << "      count: " << m_surfaces->count(cellId) << " inner: " << surface->innerMaterial().name()
	    << " outer: " << surface->outerMaterial().name() << endmsg;

  if (typeSurface.isPlane() || typeSurface.isCylinder()) {
    // scale to cov at edge
    double scale = 1.0;

    dd4hep::rec::Vector2D localPoint = surface->globalToLocal(oldPos);

    // for planar, same calculation by local is ok, but reduce repeat
    if (typeSurface.isCylinder()) {
      dd4hep::rec::Vector3D mom_ddrec(mom.x*dd4hep::GeV, mom.y*dd4hep::GeV, mom.z*dd4hep::GeV);
      double                length    = simhit.getPathLength()*dd4hep::mm;
      dd4hep::rec::Vector3D pre       = oldPos - (0.5*length)*mom_ddrec.unit();
      dd4hep::rec::Vector3D post      = oldPos + (0.5*length)*mom_ddrec.unit();
      dd4hep::rec::Vector2D localPre  = surface->globalToLocal(pre);
      dd4hep::rec::Vector2D localPost = surface->globalToLocal(post);
      localPoint = dd4hep::rec::Vector2D(0.5*(localPre.u()+localPost.u()), 0.5*(localPre.v()+localPost.v()));
      debug() << "pre: (" << pre.x() << " " << pre.y() << " " << pre.z() << " ) local (" << localPre.u() << ", " << localPre.v() << " ) "
	      << "post: (" << post.x() << " " << post.y() << " " << post.z() << " ) local (" << localPost.u() << ", " << localPost.v() << " ) " << endmsg;
    }
    //dd4hep::rec::Vector3D local3D(localPoint.u(), localPoint.v(), 0);
    // A small check, if the hit is in the boundaries:
    if (!surface->insideBounds(oldPos)) {
      double dSToHit = surface->distance(oldPos);
      debug() << " global: (" << oldPos.x()/dd4hep::mm << " " << oldPos.y()/dd4hep::mm << " " << oldPos.z()/dd4hep::mm
	      << ") local: (" << localPoint.u()/dd4hep::mm << ", " << localPoint.v()/dd4hep::mm << " )"
	      << " distance: " << dSToHit/dd4hep::mm
	      << " is not within boundaries." << endmsg;

      // FIXME: change position to path at center plane? or enlarge cov? other tracker
      if (system == CEPCConf::DetID::OTKBarrel || system == CEPCConf::DetID::OTKEndcap) {
#if 0
	dd4hep::rec::Vector3D global3D = oldPos;
	dd4hep::rec::Vector3D mom_ddrec(mom.x*dd4hep::GeV, mom.y*dd4hep::GeV, mom.z*dd4hep::GeV);
	double                length = simhit.getPathLength()*dd4hep::mm;
	dd4hep::rec::Vector3D pre    = oldPos - (0.5*length)*mom_ddrec.unit();
	dd4hep::rec::Vector3D post   = oldPos + (0.5*length)*mom_ddrec.unit();
	double dSToPre   = surface->distance(pre);
	double dPreToHit = dSToPre - dSToHit;
	global3D         = pre + (length/dPreToHit*dSToPre)*mom_ddrec.unit();
	localPoint       = surface->globalToLocal(global3D);
	scale = 1;

	debug() << " global: (" << global3D.x()/dd4hep::mm << " " << global3D.y()/dd4hep::mm << " " << global3D.z()/dd4hep::mm
		<< ") local: (" << localPoint.u()/dd4hep::mm << ", " << localPoint.v()/dd4hep::mm << ", 0 )"
		<< " distance: " << surface->distance(global3D)/dd4hep::mm
		<< " is not within boundaries." << endmsg;
#endif
      }
      //else return StatusCode::SUCCESS;
    }
    dd4hep::rec::Vector3D globalPointSmeared;//CLHEP::Hep3Vector globalPoint(pos[0],pos[1],pos[2]);
    dd4hep::rec::Vector2D localPointSmeared;

    debug() << std::setprecision(8) << " Before smearing global: (" << pos[0] << " " << pos[1] << " " << pos[2] << ") "
	    << "local: (" << localPoint.u()/dd4hep::mm << " " << localPoint.v()/dd4hep::mm << ")" << endmsg;

    unsigned tries = 0;

    bool accept_hit = false;

    // Now try to smear the hit and make sure it is still on the surface
    while (tries < 100) {
      if (tries > 0) {
	debug() << "retry smearing for side" << side << " layer"<< layer<< " module" << module
		<< " sensor" << sensor << " : retries " << tries << endmsg;
      }

      double du = m_randSvc->generator(Rndm::Gauss(0, resU))->shoot();
      double dv = m_randSvc->generator(Rndm::Gauss(0, resV))->shoot();
      localPointSmeared.u() = localPoint.u() + du;
      localPointSmeared.v() = localPoint.v() + dv;

      dd4hep::rec::Vector3D local3DSmeared(localPointSmeared.u(), localPointSmeared.v(), 0);
      globalPointSmeared = surface->localToGlobal(localPointSmeared);

      //check if hit is in boundaries
      if (surface->insideBounds(globalPointSmeared) && fabs(du) <= m_maxPull*resU && fabs(dv) <= m_maxPull*resV) {
	accept_hit = true;
	break;
      }
      else if (system == CEPCConf::DetID::ITKEndcap && fabs(du) <= m_maxPull*resU && fabs(dv) <= m_maxPull*resV) {
	accept_hit = true;
	break;
      }
      tries++;
    }

    if (accept_hit == false) {
      debug() << "hit could not be smeared within ladder after 100 tries: hit dropped"  << endmsg;
      return StatusCode::SUCCESS;
    }

    // for 1D strip measurements: set v to 0! Only the measurement in u counts!
    if(m_isStrip || (resU != 0 && resV == 0)) localPointSmeared.v() = 0. ;
    // convert back to global position for TrackerHit
    globalPointSmeared = surface->localToGlobal(localPointSmeared);

    debug() << "  After smearing global: ("
	    << globalPointSmeared.x()/dd4hep::mm <<" "<< globalPointSmeared.y()/dd4hep::mm <<" "<< globalPointSmeared.z()/dd4hep::mm << ") "
	    << "local: ("
	    << localPointSmeared.u()/dd4hep::mm << " " << localPointSmeared.v()/dd4hep::mm << ")" << endmsg;

    auto outhit = hitCol->create();

    outhit.setCellID(cellId);

    edm4hep::Vector3d smearedPos(globalPointSmeared.x()/dd4hep::mm,
				 globalPointSmeared.y()/dd4hep::mm,
				 globalPointSmeared.z()/dd4hep::mm);
    outhit.setPosition(smearedPos);
    // recover CLHEP/Geant4 unit
    resU /= dd4hep::mm;
    resV /= dd4hep::mm;

    std::bitset<32> type;
    if (typeSurface.isPlane() && m_usePlanarTag) {
      std::array<float, 6> cov;
      cov[0] = u_direction[0];
      cov[1] = u_direction[1];
      cov[2] = resU*scale;
      cov[3] = v_direction[0];
      cov[4] = v_direction[1];
      cov[5] = resV*scale;
      outhit.setCovMatrix(cov);

      type.set(CEPCConf::TrkHitTypeBit::PLANAR);
    }
    else if (typeSurface.isPlane()) {
      outhit.setCovMatrix(CEPC::ConvertToCovXYZ(resU*scale, u_direction[0], u_direction[1], resV*scale, v_direction[0], v_direction[1]));
    }
    else if (typeSurface.isCylinder()) {
      outhit.setCovMatrix(std::array<float, 6>{resU*resU*scale*scale/2., 0, resU*resU*scale*scale/2, 0, 0, resV*resV*scale*scale});
      type.set(CEPCConf::TrkHitTypeBit::CYLINDER);
    }

    if(m_isStrip || (resU != 0 && resV == 0)){
      type.set(CEPCConf::TrkHitTypeBit::ONE_DIMENSIONAL);
    }
    outhit.setType((int)type.to_ulong());
    outhit.setEDep(e);
    float dt = m_randSvc->generator(Rndm::Gauss(0, resT))->shoot();
    outhit.setTime(simhit.getTime() + dt);
    // make the relation
    auto ass = assCol->create();

    float weight = 1.0;

    debug() <<" Set relation between "
	    << " sim hit " << simhit.id()
	    << " to tracker hit " << outhit.id()
	    << " with a weight of " << weight
	    << endmsg;
#if edm4hep_VERSION >= EDM4HEP_VERSION(1, 0, 0)
    //throw std::runtime_error("The RAWHIT related interfaces are removed from TrackerHit");
    ass.setTo(simhit);
    ass.setFrom(outhit);
#else
    outhit.addToRawHits(simhit.getObjectID());
    ass.setSim(simhit);
    ass.setRec(outhit);
#endif
    ass.setWeight(weight);

    debug() << "-------------------------------------------------------" << endmsg;
  }
  else {
    fatal() << "Not plane and cylinder: " << typeSurface << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

StatusCode SmearDigiTool::AppendNoise(IDigiTool::CEPCSWTrackerHit3DCollection* hitCol) {
  for (const auto& pair/*[cellId, surface]*/ : *m_surfaces) {
    auto cellId = pair.first;
    auto surface = pair.second;
    int layer   = m_decoder->get(cellId, CEPCConf::DetCellID::layer);
    if ((m_resU.size() > 1 && layer > m_resU.size()-1) || (m_resV.size() > 1 && layer > m_resV.size()-1)) {
      fatal() << "layer exceeds resolution vector, please check input parameters ResolutionU and ResolutionV" << endmsg;
      return StatusCode::FAILURE;
    }
    double resU = (m_resU.size() > 1 ? m_resU.value().at(layer) : m_resU.value().at(0));
    double resV = (m_resV.size() > 1 ? m_resV.value().at(layer) : m_resV.value().at(0));

    auto& typeSurface = surface->type();
    double width  = surface->length_along_u();
    double length = surface->length_along_v();
    double area   = (width/dd4hep::cm)*(length/dd4hep::cm);
    debug() << "plane " << typeSurface.isPlane() << " cylinder " << typeSurface.isCylinder() << " cellId = " << cellId << " width = " << width << " length = " << length << endmsg;

    int nnoise = std::ceil(m_noiseLevel*area);
    double factor = sqrt(nnoise/m_noiseLevel/area);
    for (int inoise=0; inoise<nnoise; inoise++) {
      double u = m_randSvc->generator(Rndm::Flat(-0.5*width*factor, 0.5*width*factor))->shoot();
      double v = m_randSvc->generator(Rndm::Flat(-0.5*length*factor, 0.5*length*factor))->shoot();
      dd4hep::rec::Vector2D local(u, v);
      dd4hep::rec::Vector3D global = surface->localToGlobal(local);
      if (!surface->insideBounds(global)) continue;

      auto trkHit = hitCol->create();
      trkHit.setCellID(cellId);
      edm4hep::Vector3d pos(global.x()/dd4hep::mm, global.y()/dd4hep::mm, global.z()/dd4hep::mm);
      trkHit.setPosition(pos);

      dd4hep::rec::Vector3D uVec = surface->u(global);
      dd4hep::rec::Vector3D vVec = surface->v(global);
      float u_direction[2];
      u_direction[0] = uVec.theta();
      u_direction[1] = uVec.phi();

      float v_direction[2];
      v_direction[0] = vVec.theta();
      v_direction[1] = vVec.phi();

      debug() << " U: " << uVec << endmsg;
      debug() << " V: " << vVec << endmsg;
      debug() << " N: " << surface->normal() << endmsg;
      debug() << " O: " << surface->origin() << endmsg;

      std::bitset<32> type;
      if (typeSurface.isPlane() && m_usePlanarTag) {
	std::array<float, 6> cov;
	cov[0] = u_direction[0];
	cov[1] = u_direction[1];
	cov[2] = resU;
	cov[3] = v_direction[0];
	cov[4] = v_direction[1];
	cov[5] = resV;
	trkHit.setCovMatrix(cov);

	type.set(CEPCConf::TrkHitTypeBit::PLANAR);
      }
      else if (typeSurface.isPlane()) {
	trkHit.setCovMatrix(CEPC::ConvertToCovXYZ(resU, u_direction[0], u_direction[1], resV, v_direction[0], v_direction[1]));
      }
      else if (typeSurface.isCylinder()) {
	trkHit.setCovMatrix(std::array<float, 6>{resU*resU/2., 0, resU*resU/2, 0, 0, resV*resV});
	type.set(CEPCConf::TrkHitTypeBit::CYLINDER);
      }

      if(/*m_isStrip ||*/ (resU != 0 && resV == 0)){
	type.set(CEPCConf::TrkHitTypeBit::ONE_DIMENSIONAL);
      }
      trkHit.setType((int)type.to_ulong());

      debug() <<"Position of noise hit global: ( "
	      << global.x() <<" "<< global.y() <<" "<< global.z() << " ) "
	      << "local: ( "
	      << local.u() << " " << local.v() << " )"
	      << endmsg;
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode SmearDigiTool::finalize(){
  StatusCode sc;
  return sc;
}
