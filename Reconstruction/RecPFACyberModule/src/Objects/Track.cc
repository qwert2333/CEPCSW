#ifndef TRACK_C
#define TRACK_C

#include "Objects/Track.h"
#include <cmath>
#include "TGraph.h"

namespace Cyber{

  const double Track::B = 3.;

  int Track::getTrackerHits() const {
    if(m_track.isAvailable()) return m_track.trackerHits_size(); 
    else return 0;
  }

  int Track::trackStates_size(std::string name) const{
    std::vector<TrackState> emptyCol; emptyCol.clear(); 
    if(m_trackStates.find(name)!=m_trackStates.end()) emptyCol = m_trackStates.at(name);
    return emptyCol.size();
  }

  int Track::trackStates_size() const{
    std::vector<TrackState> emptyCol; emptyCol.clear();
    for(auto iter: m_trackStates) emptyCol.insert(emptyCol.end(), iter.second.begin(), iter.second.end());
    return emptyCol.size();
  }

  std::vector<TrackState> Track::getTrackStates(std::string name) const{
    std::vector<TrackState> emptyCol; emptyCol.clear(); 
    if(m_trackStates.find(name)!=m_trackStates.end()) emptyCol = m_trackStates.at(name);
    return emptyCol;
  }

  std::vector<TrackState> Track::getAllTrackStates() const{
    std::vector<TrackState> emptyCol; emptyCol.clear();
    for(auto iter: m_trackStates) emptyCol.insert(emptyCol.end(), iter.second.begin(), iter.second.end());
    return emptyCol;
  }

  float Track::getD0() const{
    std::vector<TrackState> trkStates = getTrackStates("Input");
    float d0 = -99.;
    for(auto it: trkStates){
      if(it.location==Cyber::TrackState::AtIP){
        d0 = it.D0;
      }
    }

    return d0;
  }  

  float Track::getZ0() const{
    std::vector<TrackState> trkStates = getTrackStates("Input");
    float z0 = -99.;
    for(auto it: trkStates){
      if(it.location==Cyber::TrackState::AtIP){
        z0 = it.Z0;
      }
    }

    return z0;
  }


  float Track::getPt() const{
    std::vector<TrackState> trkStates = getTrackStates("Input");
    float pt = -99.;
    for(auto it: trkStates){
      if(it.location==Cyber::TrackState::AtIP){ 
        pt = 1./fabs(it.Kappa);
      }
    }

    return pt;
  }

  float Track::getPz() const{
    std::vector<TrackState> trkStates = getTrackStates("Input");
    float pz = -99.;
    for(auto it: trkStates){
      if( it.location==Cyber::TrackState::AtIP ){
        pz = it.tanLambda/fabs(it.Kappa);
      }
    }

    return pz;
  }

  TVector3 Track::getP3() const{
    std::vector<TrackState> trkStates = getTrackStates("Input");
    float phi = -99;
    float pt = -99.;
    float pz = -99.;
    for(auto it: trkStates){
      if(it.location==Cyber::TrackState::AtIP){ 
        pt = 1./fabs(it.Kappa);
        phi = it.phi0;
        pz = it.tanLambda/fabs(it.Kappa);
      }
    }
  
    TVector3 p3(pt*cos(phi), pt*sin(phi), pz);
    return p3; 
  }

  float Track::getCharge() const{
    std::vector<TrackState> trkStates = getTrackStates("Input");
    float omega = -99.;
    for(auto it: trkStates){
      if(it.location==Cyber::TrackState::AtIP){ 
        omega = it.Omega;
      }
    }

    return omega/fabs(omega);
  }

  TVector3 Track::getStartPoint() const{
    std::vector<TrackState> trkStates = getTrackStates("Input");
    TVector3 startpoint (0.,0.,0.);
    for(auto it: trkStates)
      if(it.location==2) startpoint = it.referencePoint;

    return startpoint;
  }

  TVector3 Track::getEndPoint() const{
    std::vector<TrackState> trkStates = getTrackStates("Input");
    TVector3 endpoint (0.,0.,0.);
    for(auto it: trkStates)
      if(it.location==3) endpoint = it.referencePoint;

    return endpoint;
  }

  edm4hep::MCParticle Track::getLeadingMCP() const{
    float maxWeight = -1.;
    edm4hep::MCParticle mcp;
    for(auto& iter: MCParticleWeight){
      if(iter.second>maxWeight){
        mcp = iter.first;
        maxWeight = iter.second;
      }
    }

    return mcp;
  }

  float Track::getLeadingMCPweight() const{
    float maxWeight = -1.;
    for(auto& iter: MCParticleWeight){
      if(iter.second>maxWeight){
        maxWeight = iter.second;
      }
    }  
    return maxWeight;
  }

};
#endif
