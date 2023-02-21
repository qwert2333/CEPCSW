#ifndef CALO_HALFCLUSTER_H
#define CALO_HALFCLUSTER_H

#include "Objects/CaloUnit.h"
#include "Objects/Calo1DCluster.h"
#include "Objects/LongiCluster.h"
#include "Tools/TrackFitInEcal.h"

namespace PandoraPlus {

  class CaloHalfCluster {
  public: 
    CaloHalfCluster() {};
    ~CaloHalfCluster() { Clear(); };

    void Clear();
    void Clean();
    void Check(); 

    inline bool operator == (const CaloHalfCluster &x) const{
      return m_1dclusters==x.getCluster();
    }

    bool isNeighbor(const PandoraPlus::Calo1DCluster* m_1dcluster) const; 

    int getSlayer() const { return m_slayer; }

    std::vector<const CaloUnit*> getBars() const;
    std::vector<const Calo1DCluster*> getCluster() const { return m_1dclusters;};

    double getEnergy() const; 
    // TVector3 getPos() const; 

    void addUnit(const Calo1DCluster* _1dcluster);
    void setLocalMax( std::string name, std::vector<const Calo1DCluster*>& _col)
    { map_localMax[name]=_col;  }
    std::vector<const Calo1DCluster*> getLocalMaxCol(std::string name) const;

    void setLongiClusters( std::string name, std::vector<const PandoraPlus::LongiCluster*>& _cl) { map_longiClusCol[name]=_cl; }
    std::vector<const LongiCluster*> getLongiClusterCol(std::string name) const;

  private:
    int m_slayer;
    TVector3 axis;
    double trk_dr;
    double trk_dz;
    std::vector<const Calo1DCluster*> m_1dclusters; 
    std::map<std::string, std::vector<const PandoraPlus::Calo1DCluster*> > map_localMax;
    std::map<std::string, std::vector<const PandoraPlus::CaloHalfCluster*> > map_longiClusCol;  

    TrackFitInEcal* track = new TrackFitInEcal();
  };

};
#endif
