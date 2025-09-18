#ifndef PFOBJECT_H
#define PFOBJECT_H

#include "Objects/Calo3DCluster.h"
#include "Objects/Track.h"

namespace Cyber{
  class PFObject{
  public:
    PFObject () {};
    ~PFObject() { Clear(); }

    void Clear();
    std::shared_ptr<Cyber::PFObject> Clone() const; 

    void addTrack(const Track* _track);
    void addECALCluster(const Calo3DCluster* _ecal_cluster);
    void addHCALCluster(const Calo3DCluster* _hcal_cluster);

    void setPID(int _pid) { m_pid = _pid; }
    void setTrack(std::vector<const Track*> _trkCol) { m_tracks = _trkCol; } 
    void setECALCluster( std::vector<const Calo3DCluster*> _cluCol ) { m_ecal_clusters = _cluCol; }
    void setHCALCluster( std::vector<const Calo3DCluster*> _cluCol ) { m_hcal_clusters = _cluCol; }

    std::vector<const Track*> getTracks() const { return m_tracks; }
    std::vector<const Calo3DCluster*> getECALClusters() const { return m_ecal_clusters; }
    std::vector<const Calo3DCluster*> getHCALClusters() const { return m_hcal_clusters; }

    int getPID() const { return m_pid; }
    double getECALClusterEnergy() const;
    double getHCALClusterEnergy() const;
    double getTrackMomentum() const;


  private:
    int m_pid; 
    std::vector<const Track*> m_tracks;
    std::vector<const Calo3DCluster*> m_ecal_clusters;
    std::vector<const Calo3DCluster*> m_hcal_clusters;

  };

};
#endif
