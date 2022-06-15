#ifndef _TRACK_H
#define _TRACK_H

#include "Objects/TrackState.h"
#include "TVector3.h"

namespace PandoraPlus {

  class Track{
  public:
    Track() {}
    ~Track() {}; 
  void Clear() { m_trackStates.clear(); }

  int trackStates_size() const { return m_trackStates.size(); }
  TrackState getTrackStates(int _i) const { return m_trackStates[_i]; }
  std::vector<TrackState> getTrackStates() const { return m_trackStates; }
  
  void setTrackStates( std::vector<TrackState>& _states ) { m_trackStates=_states; }
  void AddTrackState( TrackState& _state ) { m_trackStates.push_back(_state); }

  private:
    static const double B ; //direction of magnetic field and charge need to check

    std::vector<TrackState> m_trackStates; 


  };

};
#endif
