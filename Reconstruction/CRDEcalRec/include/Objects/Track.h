#ifndef _TRACK_H
#define _TRACK_H

#include "Objects/TrackState.h"
#include "TVector3.h"

namespace CRDEcalEDM {

  class Track{
  public:
    Track() {}
    ~Track() {}; 
  void Clear() { m_trackStates.clear(); }

  int trackStates_size() const { return m_trackStates.size(); }
  CRDEcalEDM::TrackState getTrackStates(int _i) const { return m_trackStates[_i]; }
  std::vector<CRDEcalEDM::TrackState> getTrackStates() const { return m_trackStates; }
  
  void setTrackStates( std::vector<CRDEcalEDM::TrackState>& _states ) { m_trackStates=_states; }
  void AddTrackState( CRDEcalEDM::TrackState& _state ) { m_trackStates.push_back(_state); }

  private:
    static const double B ; //direction of magnetic field and charge need to check

    std::vector<CRDEcalEDM::TrackState> m_trackStates; 


  };

};
#endif
