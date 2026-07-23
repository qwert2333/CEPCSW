#include "TofRecAlg.h"

#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "DetInterface/IGeomSvc.h"
#include "DataHelper/HelixClass.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include <cmath>
#include "UTIL/ILDConf.h"
#include "DetIdentifier/CEPCConf.h"
#include "UTIL/CellIDEncoder.h"

using namespace edm4hep;

DECLARE_COMPONENT( TofRecAlg )

//------------------------------------------------------------------------------
TofRecAlg::TofRecAlg( const std::string& name, ISvcLocator* pSvcLocator )
    : Algorithm( name, pSvcLocator ) {
  // input
  declareProperty("CompleteTracks", m_FulTrkCol, "handler of the input track particle association collection");
  // output
  declareProperty("RecTofCollection", m_RecToFCol, "handler of the collection of reconstructed tof");
}

//------------------------------------------------------------------------------
StatusCode TofRecAlg::initialize(){
  debug() << "Initializing..." << endmsg;
  if (m_method == "Simple") {
    m_pid_svc = service("SimplePIDSvc");
  }
  else {
    m_pid_svc = nullptr;
  }
  _nEvt = 0;
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode TofRecAlg::execute(){
  
  const edm4hep::TrackCollection* fultrkcol = nullptr;
  auto rectofcol = m_RecToFCol.createAndPut();

  try {
    fultrkcol = m_FulTrkCol.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Complete track particle association collection " << m_FulTrkCol.fullKey() << " is unavailable in event " << _nEvt << endmsg;
    _nEvt++;
    return StatusCode::SUCCESS;
  }

  for(auto track : *fultrkcol){
      hasFTDHit = false;
      hasSETHit = false;

      FindToFHits( track, hasFTDHit, hasSETHit, Toft, Tofx, Tofy, Tofz );

      if ( !hasSETHit && !hasFTDHit ) {
        debug() << "No SET or FTD hit found for this track" << endmsg;
        continue;
      }

      for( std::vector<edm4hep::TrackState>::const_iterator it = track.trackStates_begin(); it != track.trackStates_end(); it++ ){
          edm4hep::TrackState trackState = *it;
          if ( trackState.location == TrackState::AtIP ){
            position0x = trackState.referencePoint[0];
            position0y = trackState.referencePoint[1];
            position0z = trackState.referencePoint[2];
            debug()<<"position0x = "<<position0x<<" position0y = "<<position0y<<" position0z = "<<position0z<<endmsg;
            first_hit_pos = {position0x, position0y, position0z};//for trajectory length calculation
          }
          if ( trackState.location == TrackState::AtLastHit ){
            position1x = trackState.referencePoint[0];
            position1y = trackState.referencePoint[1];
            position1z = trackState.referencePoint[2];
            debug()<<"position1x = "<<position1x<<" position1y = "<<position1y<<" position1z = "<<position1z<<endmsg;
            last_hit_pos = {position1x, position1y, position1z};
          }
      }// find track states at IP and OTK

      // smearing the measured time by gaussian 0.05 for intrinsic and 0.02 for bunchcrossing, giving the measured terr a constant
      Tofterr = std::sqrt( tof_resolution * tof_resolution + bunchcrossing * bunchcrossing );
      Toft += normal_distribution_tof( generator_tof );
      Toft += normal_distribution_bunch( generator_bunch );

      // length of the trajectory using simple pid svc
      double _length, _p_atIP, _costheta;
      m_pid_svc->getTrackPars( first_hit_pos, last_hit_pos, track, TrackState::AtIP, _length, _p_atIP, _costheta );
      debug()<<"_length = "<<_length<<" _p_atIP = "<<_p_atIP<<" _costheta = "<<_costheta<<endmsg;

      std::array<float, 5> exp_times = { -99.9, -99.9, -99.9, -99.9, -99.9 };
      // evaluate expected times and errors for e, mu, pi, k, proton hypotheses
      for (auto it = masses.begin(); it != masses.end(); ++it){
        int idx = it->first;
        double assum = it->second;//mass in GeV
        double energy = std::sqrt( _p_atIP * _p_atIP + assum * assum );
        double velocity = ( _p_atIP / energy ) * c_mm_ns;//mm/ns
        double t_exp = _length / velocity;//ns
        exp_times[idx] = t_exp;
        debug()<<"    idx : "<< idx <<" assum: "<< assum <<" length: "<<_length<<" exp_t: "<<t_exp<<" obs_t: "<<Toft<<endmsg;
      }//end loop over masses

      auto rectof = rectofcol->create();
      rectof.setTime(Toft);
      rectof.setTimeExp( exp_times );
      rectof.setPosition({Tofx, Tofy, Tofz});
      rectof.setSigma(Tofterr);
      float _flength = float(_length);//avoid the warning of narrowing conversion
      rectof.setPathLength({_flength, _flength, _flength, _flength, _flength});//all particles have the same path length for now, becasue the path length is calculated by the simple pid svc with helix.
      rectof.setTrack(track);

  }//end loop over tracks
  _nEvt++;
  return StatusCode::SUCCESS;
}// end execute

//------------------------------------------------------------------------------
StatusCode TofRecAlg::finalize(){
  debug() << "Finalizing..." << endmsg;
  return StatusCode::SUCCESS;
}

void TofRecAlg::FindToFHits(const edm4hep::Track& _track, bool& _hasFTDHit, bool& _hasSETHit, double& _Toft, double& _Tofx, double& _Tofy, double& _Tofz){
  int NHits = _track.trackerHits_size();
  int det_id = 0;
  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string );      
  for ( int hit=0; hit<NHits; hit++ ) {
    edm4hep::TrackerHit SimTHit = _track.getTrackerHits(hit);
    encoder.setValue(SimTHit.getCellID());
    det_id  = encoder[lcio::ILDCellID0::subdet];
    if ( det_id == lcio::ILDDetID::SET ){
      _hasSETHit = true;
      _Toft = SimTHit.getTime();
      _Tofx = SimTHit.getPosition()[0];
      _Tofy = SimTHit.getPosition()[1];
      _Tofz = SimTHit.getPosition()[2];
    }
    else if( det_id == lcio::ILDDetID::ETD ){
      _Toft = SimTHit.getTime();
      _Tofx = SimTHit.getPosition()[0];
      _Tofy = SimTHit.getPosition()[1];
      _Tofz = SimTHit.getPosition()[2];
      _hasFTDHit = true;
    }//find etd hit
  }//end track hits loop
}
