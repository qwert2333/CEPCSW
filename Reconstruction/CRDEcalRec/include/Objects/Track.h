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

  int trackStates_size(std::string name) const;
  int trackStates_size() const ;
  std::vector<TrackState> getTrackStates(std::string name) const;
  std::vector<TrackState> getAllTrackStates() const;
  std::map<std::string, std::vector<TrackState> > getTrackStatesMap() const { return m_trackStates; }
  std::vector<PandoraPlus::CaloHalfCluster*> getAssociatedHalfClustersU() const { return m_halfClusterUCol; }  
  std::vector<PandoraPlus::CaloHalfCluster*> getAssociatedHalfClustersV() const { return m_halfClusterVCol; }  

  void setTrackStates( std::string name, std::vector<TrackState>& _states ) { m_trackStates[name]=_states; }
  void addAssociatedHalfClusterU( PandoraPlus::CaloHalfCluster* _cl ) { m_halfClusterUCol.push_back(_cl); }
  void addAssociatedHalfClusterV( PandoraPlus::CaloHalfCluster* _cl ) { m_halfClusterVCol.push_back(_cl); }

  void setType(int _type) { m_type=_type; }
  int getType() const { return m_type; }

  private:
    static const double B ; //direction of magnetic field and charge need to check

    std::map<std::string, std::vector<TrackState> > m_trackStates; // name = Input, Ecal, Hcal
    std::vector<PandoraPlus::CaloHalfCluster*> m_halfClusterUCol; 
    std::vector<PandoraPlus::CaloHalfCluster*> m_halfClusterVCol; 

    int m_type;
  };

};
#endif
