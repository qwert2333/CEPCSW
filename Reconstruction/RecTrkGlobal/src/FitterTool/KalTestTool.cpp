#include "KalTestTool.h"

#include "TrackSystemSvc/ITrackSystemSvc.h"
#include "TrackSystemSvc/MarlinTrkUtils.h"
#include "TrackSystemSvc/IMarlinTrack.h"
#include "DetInterface/IGeomSvc.h"
#include "DetIdentifier/CEPCConf.h"

#include "DataHelper/Navigation.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"

#include "edm4hep/TrackerHit.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/MutableTrack.h"
#include "edm4hep/EDM4hepVersion.h"
#include "podio/podioVersion.h"

DECLARE_COMPONENT(KalTestTool)

StatusCode KalTestTool::initialize() {
  StatusCode sc;
  always() << m_fitterName << endmsg;
  if (m_fitterName=="KalTest" || m_fitterName=="DDKalTest") {
    auto _trackSystemSvc = service<ITrackSystemSvc>("TrackSystemSvc");
    if (!_trackSystemSvc) {
      error() << "Failed to find TrackSystemSvc ..." << endmsg;
      return StatusCode::FAILURE;
    }
    m_factoryMarlinTrk = _trackSystemSvc->getTrackSystem(this, m_fitterName.value());
    m_factoryMarlinTrk->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useQMS, m_useQMS);
    m_factoryMarlinTrk->setOption(MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, m_usedEdx);
    m_factoryMarlinTrk->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, m_useSmoothing);
    m_factoryMarlinTrk->init();
  }
  else {
    error() << "fitter " << m_fitterName << " has not been imported" << endmsg;
    return StatusCode::FAILURE;
  }

  auto _geomSvc = service<IGeomSvc>("GeomSvc");
  if ( !_geomSvc ) {
    error() << "Failed to find GeomSvc ..." << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_magneticField.value()==0) {
    const dd4hep::Direction& field = _geomSvc->lcdd()->field().magneticField(dd4hep::Position(0,0,0));
    double Bz = field.z()/dd4hep::tesla;
    if (Bz==0) {
      error() << "magnetic field = 0, KalmanFilter cannot run" << endmsg;
      return StatusCode::FAILURE;
    }
    m_magneticField = Bz;
  }

  return sc;
}

StatusCode KalTestTool::finalize() {
  StatusCode sc;
  return sc;
}

int KalTestTool::Fit(edm4hep::MutableTrack track, std::vector<edm4hep::TrackerHit>& trackHits,
			  const decltype(edm4hep::TrackState::covMatrix)& covMatrix, double maxChi2perHit, bool backward) {
  if (m_hitsInFit.size()!=0 || m_outliers.size()!=0) {
    error() << "Important! vector not clear, still store the data of last event!" << endmsg;
    return 0;
  }

  if (m_fitterName=="KalTest" || m_fitterName=="DDKalTest") {
    debug() << "start..." << endmsg;
    std::shared_ptr<MarlinTrk::IMarlinTrack> marlinTrack(m_factoryMarlinTrk->createTrack());
    debug() << "created MarlinKalTestTrack" << endmsg;
    int status = this->createFinalisedTrack(marlinTrack.get(), trackHits, &track, backward, covMatrix, m_magneticField, maxChi2perHit);

    marlinTrack->getHitsInFit(m_hitsInFit);
    marlinTrack->getOutliers(m_outliers);
    return status;
  }
  
  error() << "Don't support the Fitter " << m_fitterName << endmsg;
  return 0;
}

int KalTestTool::Fit(edm4hep::MutableTrack track, std::vector<edm4hep::TrackerHit>& trackHits,
                          edm4hep::TrackState trackState, double maxChi2perHit, bool backward) {
  if (m_hitsInFit.size()!=0 || m_outliers.size()!=0) {
    error() << "Important! vector not clear, still store the data of last event!" << endmsg;
    return 0;
  }

  if (m_fitterName=="KalTest" || m_fitterName=="DDKalTest") {
    debug() << "start..." << endmsg;
    std::shared_ptr<MarlinTrk::IMarlinTrack> marlinTrack(m_factoryMarlinTrk->createTrack());
    debug() << "created MarlinKalTestTrack" << endmsg;
    int status = this->createFinalisedTrack(marlinTrack.get(), trackHits, &track, backward, &trackState, m_magneticField, maxChi2perHit);

    marlinTrack->getHitsInFit(m_hitsInFit);
    marlinTrack->getOutliers(m_outliers);
    return status;
  }

  error() << "Don't support the Fitter " << m_fitterName << endmsg;
  return 0;
}

// from MarlinTrk
int KalTestTool::createFinalisedTrack(MarlinTrk::IMarlinTrack* marlinTrk, std::vector<edm4hep::TrackerHit>& hit_list, edm4hep::MutableTrack* track, bool fit_backwards,
				      const decltype(edm4hep::TrackState::covMatrix)& initial_cov_for_prefit, float bfield_z, double maxChi2Increment) {
  ///////////////////////////////////////////////////////
  // check inputs
  ///////////////////////////////////////////////////////
  if ( hit_list.empty() ) return MarlinTrk::IMarlinTrack::bad_intputs;

  if( track == 0 ){
    throw std::runtime_error( "MarlinTrk::finaliseLCIOTrack: TrackImpl == NULL " ) ;
  }

  int return_error = 0;
  ///////////////////////////////////////////////////////
  // produce prefit parameters
  ///////////////////////////////////////////////////////

  edm4hep::TrackState pre_fit ;

  //std::cout << "debug:=====================before createPrefit" << std::endl;
  return_error = MarlinTrk::createPrefit(hit_list, &pre_fit, bfield_z, fit_backwards);
  //std::cout << "debug:=====================after createPrefit return code=" << return_error << std::endl;
  pre_fit.covMatrix = initial_cov_for_prefit;

  ///////////////////////////////////////////////////////
  // use prefit parameters to produce Finalised track
  ///////////////////////////////////////////////////////

  if( return_error == 0 ) {

    return_error = createFinalisedTrack( marlinTrk, hit_list, track, fit_backwards, &pre_fit, bfield_z, maxChi2Increment);

  }
  else {
    warning() << "createFinalisedLCIOTrack : Prefit failed error = " << MarlinTrk::errorCode(return_error) << endmsg;
  }
  return return_error;
}

// from MarlinTrk
int KalTestTool::createFinalisedTrack(MarlinTrk::IMarlinTrack* marlinTrk, std::vector<edm4hep::TrackerHit>& hit_list, edm4hep::MutableTrack* track, bool fit_backwards,
				      edm4hep::TrackState* pre_fit, float bfield_z, double maxChi2Increment) {

  ///////////////////////////////////////////////////////
  // check inputs
  ///////////////////////////////////////////////////////
  if ( hit_list.empty() ) return MarlinTrk::IMarlinTrack::bad_intputs;

  if( track == 0 ){
    throw std::runtime_error("MarlinTrk::finaliseLCIOTrack: TrackImpl == NULL ");
  }

  if( pre_fit == 0 ){
    throw std::runtime_error("MarlinTrk::finaliseLCIOTrack: TrackStateImpl == NULL ");
  }


  int fit_status = MarlinTrk::createFit(hit_list, marlinTrk, pre_fit, bfield_z, fit_backwards, maxChi2Increment);

  if (fit_status != MarlinTrk::IMarlinTrack::success) {
    debug() << "createFinalisedTrack fit failed: fit_status = " << MarlinTrk::errorCode(fit_status) << endmsg;
    return fit_status;
  }

  int error = finaliseTrack(marlinTrk, track, hit_list, fit_backwards);
  debug() << "finaliseTrack. status = " << MarlinTrk::errorCode(error) << endmsg;

  return error;
}

// from MarlinTrk
int KalTestTool::finaliseTrack(MarlinTrk::IMarlinTrack* marlintrk, edm4hep::MutableTrack* track, std::vector<edm4hep::TrackerHit>& hit_list, bool fit_backwards,
			       edm4hep::TrackState* atLastHit, edm4hep::TrackState* atCaloFace) {

  ///////////////////////////////////////////////////////
  // check inputs
  ///////////////////////////////////////////////////////
  if (marlintrk == 0) {
    throw std::runtime_error("MarlinTrk::finaliseLCIOTrack: IMarlinTrack == NULL ");
  }

  if (track == 0) {
    throw std::runtime_error("MarlinTrk::finaliseLCIOTrack: TrackImpl == NULL ");
  }

  if (atCaloFace && atLastHit == 0) {
    throw std::runtime_error("MarlinTrk::finaliseLCIOTrack: atLastHit == NULL ");
  }

  if (atLastHit && atCaloFace == 0) {
    throw std::runtime_error("MarlinTrk::finaliseLCIOTrack: atCaloFace == NULL ");
  }

  ///////////////////////////////////////////////////////
  // error to return if any
  ///////////////////////////////////////////////////////
  int return_error = 0;

  int    ndf  = 0;
  double chi2 = -DBL_MAX;

  /////////////////////////////////////////////////////////////
  // First check NDF to see if it make any sense to continue.
  // The track will be dropped if the NDF is less than 0
  /////////////////////////////////////////////////////////////

  return_error = marlintrk->getNDF(ndf);

  if (return_error != MarlinTrk::IMarlinTrack::success) {
    debug() << "getNDF returns " << MarlinTrk::errorCode(return_error) << endmsg;
    return return_error;
  }
  else if (ndf < 0) {
    debug() << "number of degrees of freedom less than 0 track dropped : NDF = " << ndf << endmsg;
    return MarlinTrk::IMarlinTrack::error;
  }
  else {
    debug() << "NDF = " << ndf << endmsg;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // get the list of hits used in the fit
  // add these to the track, add spacepoints as long as at least on strip hit is used.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<std::pair<edm4hep::TrackerHit, double> > hits_in_fit;
  std::vector<std::pair<edm4hep::TrackerHit, double> > outliers;
  std::vector<edm4hep::TrackerHit> used_hits;

  hits_in_fit.reserve(300);
  outliers.reserve(300);

  marlintrk->getHitsInFit(hits_in_fit);
  marlintrk->getOutliers(outliers);

  ///////////////////////////////////////////////
  // now loop over the hits provided for fitting
  // we do this so that the hits are added in the
  // order in which they have been fitted
  ///////////////////////////////////////////////

  for ( unsigned ihit = 0; ihit < hit_list.size(); ++ihit) {

    edm4hep::TrackerHit trkHit = hit_list[ihit];

    std::bitset<32> type = trkHit.getType();
    if (type.test(CEPCConf::TrkHitTypeBit::COMPOSITE_SPACEPOINT)) { //it is a composite spacepoint
      //std::cout << "Error: space point is not still valid! pelease wait updating..." <<std::endl;
      //exit(1);
      // get strip hits
#if edm4hep_VERSION >= EDM4HEP_VERSION(1, 0, 0)
      throw std::runtime_error("The RAWHIT related interfaces are removed from TrackerHit");
#else
      int nRawHit = trkHit.rawHits_size();
      for( unsigned k=0; k< nRawHit; k++ ){
	auto rawHitOpt = Navigation::Instance()->GetTrackerHit(trkHit.getRawHits(k));

        if (!rawHitOpt) {
            throw std::runtime_error("Failed to find raw TrackerHit from ObjectID");
        }

        edm4hep::TrackerHit rawHit = *rawHitOpt;
	bool is_outlier = false;
	// here we loop over outliers as this will be faster than looping over the used hits
	for ( unsigned ohit = 0; ohit < outliers.size(); ++ohit) {
	  if ( rawHit.id() == outliers[ohit].first.id() ) {
	    is_outlier = true;
	    break; // break out of loop over outliers
	  }
	}
	if (is_outlier == false) {
	  used_hits.push_back(hit_list[ihit]);
	  track->addToTrackerHits(used_hits.back());
	  break; // break out of loop over rawObjects
	}
      }
#endif
    }
    else {
      bool is_outlier = false;
      // here we loop over outliers as this will be faster than looping over the used hits
      for ( unsigned ohit = 0; ohit < outliers.size(); ++ohit) {
	if ( trkHit == outliers[ohit].first ) {
	  is_outlier = true;
	  break; // break out of loop over outliers
	}
      }
      if (is_outlier == false) {
	used_hits.push_back(hit_list[ihit]);
	track->addToTrackerHits(used_hits.back());
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // We now need to find out at which point the fit is constrained
  // and therefore be able to provide well formed (pos. def.) cov. matrices
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////
  // first hit
  ///////////////////////////////////////////////////////

  edm4hep::TrackState* trkStateAtFirstHit = new edm4hep::TrackState() ;
  edm4hep::TrackerHit firstHit = (fit_backwards == MarlinTrk::IMarlinTrack::backward) ? hits_in_fit.back().first : hits_in_fit.front().first;

  ///////////////////////////////////////////////////////
  // last hit
  ///////////////////////////////////////////////////////

  edm4hep::TrackState* trkStateAtLastHit = new edm4hep::TrackState() ;

  edm4hep::TrackerHit lastHit = (fit_backwards == MarlinTrk::IMarlinTrack::backward) ? hits_in_fit.front().first : hits_in_fit.back().first;

#if PODIO_BUILD_VERSION < PODIO_VERSION(0, 17, 4)
  edm4hep::TrackerHit last_constrained_hit(0);
#else
  auto last_constrained_hit = edm4hep::TrackerHit::makeEmpty();
#endif

  marlintrk->getTrackerHitAtPositiveNDF(last_constrained_hit);

  debug() << "finaliseLCIOTrack:  firstHit: " << firstHit.getCellID() << " " << firstHit.getPosition()
	  << " lastHit: " << lastHit.getCellID() << " " << lastHit.getPosition()
	  << " last constrained hit: " << last_constrained_hit.getCellID() << " " << last_constrained_hit.getPosition()
	  << " fit direction is backward : " << fit_backwards << endmsg;

  return_error = marlintrk->smooth(lastHit);
  //return_error = marlintrk->smooth(last_constrained_hit);

  if (return_error != MarlinTrk::IMarlinTrack::success) {
    debug() << "return_code for smoothing to " << lastHit << " = " << MarlinTrk::errorCode(return_error) << " NDF = " << ndf << endmsg;
    delete trkStateAtFirstHit;
    delete trkStateAtLastHit;
    return return_error ;
  }

  ///////////////////////////////////////////////////////
  // first create trackstate at IP
  ///////////////////////////////////////////////////////
  const edm4hep::Vector3d point; // nominal IP

  edm4hep::TrackState* trkStateIP = new edm4hep::TrackState();

  ///////////////////////////////////////////////////////
  // make sure that the track state can be propagated to the IP
  ///////////////////////////////////////////////////////

  //FIXME: AidaTT not used now, if to use AidaTT, these code need to update and apply
  //MarlinTrk::IMarlinTrkSystem* trksystem =  MarlinTrk::Factory::getCurrentMarlinTrkSystem() ;
  //bool usingAidaTT = ( trksystem->name() == "AidaTT" ) ;
  bool usingAidaTT = false;
  if (fit_backwards == MarlinTrk::IMarlinTrack::backward ||  usingAidaTT ) {
    return_error = marlintrk->propagate(point, firstHit, *trkStateIP, chi2, ndf ) ;
  }
  else {
    // if we fitted forward, we start from the last_constrained hit
    // and then add the last inner hits with a Kalman step ...

    // create a temporary IMarlinTrack
    auto mTrk = std::shared_ptr<MarlinTrk::IMarlinTrack>( m_factoryMarlinTrk->createTrack() );

    edm4hep::TrackState ts;

    auto hI = hits_in_fit.rbegin();
    
    double chi2Tmp = 0;
    int    ndfTmp  = 0;
    //return_error = marlintrk->getTrackState( last_constrained_hit, ts, chi2, ndf);
    return_error = marlintrk->getTrackState( lastHit, ts, chi2, ndf);

    debug() << "-- TrackState at last constrained hit : " << ts << endmsg;

    //need to add a dummy hit to the track
    //mTrk->addHit(last_constrained_hit);
    mTrk->addHit(lastHit);

    double _bfield = m_magneticField;
    // fixme: the implementation for DDKalTest does no longer need this value but the IMarlinTrk interface is not yet changed
    mTrk->initialise(ts, _bfield, fit_backwards);

    //while (hI->first.id() != last_constrained_hit.id()) {
    while (hI->first.id() != lastHit.id()) {
      debug() << "-- hit in reverse_iterator : " << hI->first.getCellID() << " " << hI->first.getPosition() << endmsg;
      ++hI;
    }

    ++hI;

    while (hI != hits_in_fit.rend()) {
      auto hit = (*hI).first;

      double deltaChi;
      double maxChi2Increment = DBL_MAX; // not apply chisquare cut, since hits_in_fit have past in filter

      int addHit = mTrk->addAndFit(hit, deltaChi, maxChi2Increment);

      debug() << "-- hit id: " << hit.id() << " cellId: " << hit.getCellID() << " pos: " << hit.getPosition()
	      << "  added : " << MarlinTrk::errorCode(addHit) << " deltaChi2: " << deltaChi << endmsg;

      if (addHit !=  MarlinTrk::IMarlinTrack::success) {
	debug() << "-- could not add inner hit to track !!! " << maxChi2Increment << endmsg;
      }

      ++hI;
    }//------------------------------------

    debug() << "-- temporary kaltest track for track state at the IP: " << mTrk->toString() << endmsg;

    // now propagate the temporary track to the IP
    return_error = mTrk->propagate(point, firstHit, *trkStateIP, chi2Tmp, ndfTmp);

    debug() << "-- propagated temporary track from first hit to IP : " <<  (*trkStateIP) << endmsg;
    //FIXME: if forward, better by add the last inner hits with a Kalman step from the last_constrained hit
    //return_error = marlintrk->propagate(point, firstHit, *trkStateIP, chi2, ndf ) ;
  }
  if (return_error != MarlinTrk::IMarlinTrack::success) {
    debug() << "-- return_code for propagation = " << MarlinTrk::errorCode(return_error) << " NDF = " << ndf << std::endl;
    delete trkStateIP;
    delete trkStateAtFirstHit;
    delete trkStateAtLastHit;

    return return_error ;
  }

  trkStateIP->location = edm4hep::TrackState::AtIP;
  track->addToTrackStates(*trkStateIP);
  track->setChi2(chi2);
  track->setNdf(ndf);

  ///////////////////////////////////////////////////////
  // set the track states at the first and last hits
  ///////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////
  // @ first hit
  ///////////////////////////////////////////////////////

  return_error = marlintrk->getTrackState(firstHit, *trkStateAtFirstHit, chi2, ndf ) ;

  if ( return_error == 0 ) {
    trkStateAtFirstHit->location = edm4hep::TrackState::AtFirstHit;
    track->addToTrackStates(*trkStateAtFirstHit);
    debug() << ">>>> create TrackState AtFirstHit" << (*trkStateAtFirstHit) << endmsg;
  }
  else {
    warning() << ">>>> could not get TrackState at First Hit " << firstHit.getCellID() << " " << firstHit.getPosition() << endmsg;
  }

  double r_first = firstHit.getPosition()[0]*firstHit.getPosition()[0] + firstHit.getPosition()[1]*firstHit.getPosition()[1];

  //#if podio_VERSION >= PODIO_VERSION(1, 0, 0)
#if edm4hep_VERSION >= EDM4HEP_VERSION(1, 0, 0)
  //throw std::runtime_error("The setRadiusOfInnermostHit interface is removed from Track");
#else
  track->setRadiusOfInnermostHit(sqrt(r_first));
#endif

  if ( atLastHit == 0 && atCaloFace == 0 ) {

    ///////////////////////////////////////////////////////
    // @ last hit
    ///////////////////////////////////////////////////////

    debug() << ">>>> create TrackState AtLastHit : using trkhit " << lastHit.getCellID() << " " << lastHit.getPosition()
	    << " from last_constrained_hit " << last_constrained_hit.getCellID() << " " << last_constrained_hit.getPosition() << endmsg;

    edm4hep::Vector3d last_hit_pos(lastHit.getPosition());

    return_error = marlintrk->propagate(last_hit_pos, last_constrained_hit, *trkStateAtLastHit, chi2, ndf);

    if ( return_error == 0 ) {
      trkStateAtLastHit->location = edm4hep::TrackState::AtLastHit;
      track->addToTrackStates(*trkStateAtLastHit);

      debug() << ">>>> createTrackStateAtLastHit OK: " << (*trkStateAtLastHit) << endmsg;
    }
    else {
      error() << ">>>> could not get TrackState at Last Hit " << lastHit.getCellID() << " " << lastHit.getPosition()
	      << " from last_constrained_hit " << last_constrained_hit.getCellID() << " " << last_constrained_hit.getPosition() << endmsg;
      //delete trkStateAtLastHit;
    }

    ///////////////////////////////////////////////////////
    // set the track state at Calo Face
    ///////////////////////////////////////////////////////

    edm4hep::TrackState trkStateCalo;
    bool tanL_is_positive = trkStateIP->tanLambda >0 ;

    return_error = MarlinTrk::createTrackStateAtCaloFace(marlintrk, &trkStateCalo, last_constrained_hit, tanL_is_positive);

    if (return_error == 0) {
      trkStateCalo.location = edm4hep::TrackState::AtCalorimeter;
      track->addToTrackStates(trkStateCalo);

      debug() << ">>>> createTrackStateAtCaloFace OK: " << trkStateCalo << endmsg;
    }
    else {
      if (msgLevel(MSG::DEBUG)) {
	debug() << ">>>> could not get TrackState at Calo Face " << endmsg;
	if (last_constrained_hit.isAvailable()) {
	  auto pos = last_constrained_hit.getPosition();
	  debug() << ">>>>   last_constrained_hit = " << pos.x << "," << pos.y << "," << pos.z << endmsg;
	}
	else {
	  debug() << ">>>>   last_constrained_hit not Available" << endmsg;
	}
      }
      //delete trkStateCalo;
    }
  }
  else {
    track->addToTrackStates(*atLastHit);
    track->addToTrackStates(*atCaloFace);
    //delete trkStateAtLastHit;
  }

  if(trkStateAtFirstHit) delete trkStateAtFirstHit;
  if(trkStateAtLastHit)  delete trkStateAtLastHit;
  if(trkStateIP)         delete trkStateIP;
  ///////////////////////////////////////////////////////
  // done
  ///////////////////////////////////////////////////////
  return return_error;
}
