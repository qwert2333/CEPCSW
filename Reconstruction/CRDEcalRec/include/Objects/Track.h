#ifndef _TRACK_H
#define _TRACK_H

#include "Objects/TrackState.h"
#include "Objects/CaloHalfCluster.h"
#include "TVector3.h"

namespace PandoraPlus {

  class CaloHalfCluster;
  class Track{
  public:
    Track() {}
    ~Track() { Clear(); }; 
  void Clear() { m_trackStates.clear(); m_halfClusterUCol.clear(); m_halfClusterVCol.clear(); }

  int trackStates_size() const { return m_trackStates.size(); }
  TrackState getTrackStates(int _i) const { return m_trackStates[_i]; }
  std::vector<TrackState> getTrackStates() const { return m_trackStates; }
  std::vector<PandoraPlus::CaloHalfCluster*> getAssociatedHalfClustersU() const { return m_halfClusterUCol; }  
  std::vector<PandoraPlus::CaloHalfCluster*> getAssociatedHalfClustersV() const { return m_halfClusterVCol; }  

  void setTrackStates( std::vector<TrackState>& _states ) { m_trackStates=_states; }
  void AddTrackState( TrackState& _state ) { m_trackStates.push_back(_state); }
  void addAssociatedHalfClusterU( PandoraPlus::CaloHalfCluster* _cl ) { m_halfClusterUCol.push_back(_cl); }
  void addAssociatedHalfClusterV( PandoraPlus::CaloHalfCluster* _cl ) { m_halfClusterVCol.push_back(_cl); }

  void setType(int _type) { m_type=_type; }
  int getType() const { return m_type; }

  private:
    static const double B ; //direction of magnetic field and charge need to check

    std::vector<TrackState> m_trackStates; 
    std::vector<PandoraPlus::CaloHalfCluster*> m_halfClusterUCol; 
    std::vector<PandoraPlus::CaloHalfCluster*> m_halfClusterVCol; 

    int m_type;
  };

};
#endif
