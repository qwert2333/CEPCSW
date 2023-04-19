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


};
#endif