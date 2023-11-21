#ifndef PFOBJECT_H
#define PFOBJECT_H

#include "Objects/Calo3DCluster.h"
#include "Objects/Track.h"

namespace PandoraPlus{
  class PFObject{
  public:
    PFObject () {};
    ~PFObject() { Clear(); }

    void Clear();

    void addTrack(const Track* _track);
    void addECALCluster(const Calo3DCluster* _ecal_cluster);
    void addHCALCluster(const Calo3DCluster* _hcal_cluster);

    std::vector<const Track*> getTracks() const { return m_tracks; }
    std::vector<const Calo3DCluster*> getECALClusters() const { return m_ecal_clusters; }
    std::vector<const Calo3DCluster*> getHCALClusters() const { return m_hcal_clusters; }

  private:
    std::vector<const Track*> m_tracks;
    std::vector<const Calo3DCluster*> m_ecal_clusters;
    std::vector<const Calo3DCluster*> m_hcal_clusters;

  };

};
#endif