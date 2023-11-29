#ifndef TRACK_C
#define TRACK_C

#include "Objects/Track.h"
#include <cmath>
#include "TGraph.h"

namespace PandoraPlus{

  const double Track::B = 3.;

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

  float Track::getPt() const{
    std::vector<TrackState> trkStates = getTrackStates("Input");
    float pt = -99.;
    for(auto it: trkStates){
      if(it.location==4){
        pt = 1./it.Kappa;
      }
    }

    return pt;
  }

  float Track::getPz() const{
    std::vector<TrackState> trkStates = getTrackStates("Input");
    float pz = -99.;
    for(auto it: trkStates){
      if(it.location==4){
        pz = it.tanLambda/it.Kappa;
      }
    }

    return pz;
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
