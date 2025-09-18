#ifndef _TRACK_H
#define _TRACK_H

#include "edm4hep/Track.h"
#include "Objects/TrackState.h"
#include "Objects/CaloHalfCluster.h"
#include "TVector3.h"

namespace Cyber {

  class CaloHalfCluster;
  class Track{
  public:
    Track() {}
    ~Track() { Clear(); }; 
  void Clear() { m_trackStates.clear(); m_halfClusterUCol.clear(); m_halfClusterVCol.clear(); }

  void setOriginTrack(edm4hep::Track& _trk) { m_track = _trk; }
  int getTrackerHits() const;
  int trackStates_size(std::string name) const;
  int trackStates_size() const ;
  edm4hep::Track getOriginTrack() const { return m_track; }
  std::vector<TrackState> getTrackStates(std::string name) const;
  std::vector<TrackState> getAllTrackStates() const;
  std::map<std::string, std::vector<TrackState> > getTrackStatesMap() const { return m_trackStates; }
  std::vector<Cyber::CaloHalfCluster*> getAssociatedHalfClustersU() const { return m_halfClusterUCol; }  
  std::vector<Cyber::CaloHalfCluster*> getAssociatedHalfClustersV() const { return m_halfClusterVCol; }  
  std::vector< std::pair<edm4hep::MCParticle, float> > getLinkedMCP() const { return MCParticleWeight; }
  edm4hep::MCParticle getLeadingMCP() const;
  float getD0() const;
  float getZ0() const;
  float getLeadingMCPweight() const;
  float getPt() const;
  float getPz() const;
  float getMomentum() const { return sqrt( getPt()*getPt() + getPz()*getPz() ); } 
  TVector3 getP3() const; 
  float getCharge() const;
  TVector3 getStartPoint() const;
  TVector3 getEndPoint() const;

  void setTrackStates( std::string name, std::vector<TrackState>& _states ) { m_trackStates[name]=_states; }
  void addAssociatedHalfClusterU( Cyber::CaloHalfCluster* _cl ) { m_halfClusterUCol.push_back(_cl); }
  void addAssociatedHalfClusterV( Cyber::CaloHalfCluster* _cl ) { m_halfClusterVCol.push_back(_cl); }
  void addLinkedMCP( std::pair<edm4hep::MCParticle, float> _pair ) { MCParticleWeight.push_back(_pair); }
  void setLinkedMCP( std::vector<std::pair<edm4hep::MCParticle, float>> _pairVec ) { MCParticleWeight.clear(); MCParticleWeight = _pairVec; }

  void setPID(int _pid) {m_pid = _pid; }
  void setType(int _type) { m_type=_type; }
  int getType() const { return m_type; }
  int getPID() const { return m_pid; }

  private:
    edm4hep::Track m_track;
    static const double B ; //direction of magnetic field and charge need to check

    std::map<std::string, std::vector<TrackState> > m_trackStates; // name = Input, Ecal, Hcal
    std::vector<Cyber::CaloHalfCluster*> m_halfClusterUCol; 
    std::vector<Cyber::CaloHalfCluster*> m_halfClusterVCol; 

    int m_type;
    int m_pid;
    std::vector< std::pair<edm4hep::MCParticle, float> > MCParticleWeight;
  };

};
#endif
