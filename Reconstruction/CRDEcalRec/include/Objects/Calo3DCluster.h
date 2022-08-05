#ifndef CALO_3DCLUSTER_H
#define CALO_3DCLUSTER_H

#include "Objects/Calo2DCluster.h"
#include "Objects/CaloTower.h"

namespace PandoraPlus {

  class Calo3DCluster {
  public: 
    Calo3DCluster() {};
    ~Calo3DCluster() { Clear(); };

    void Clear();
    void Clean();
    void Check(); 
    void Clean2DClusters();
    void CleanLongiClusters();

    inline bool operator == (const Calo3DCluster &x) const{
      return ( m_2dclusters == x.getCluster()  );
    }

    std::vector<int> getModules() const { return m_modules; }
    std::vector<int> getParts() const { return m_parts; }
    std::vector<int> getStaves() const { return m_staves; }

    bool isNeighbor(const PandoraPlus::Calo2DCluster* m_2dcluster) const; 
    std::vector<const Calo2DCluster*> getCluster() const { return m_2dclusters; }

    std::vector<const CaloBarShower*> getLocalMaxUCol(std::string name) const;
    std::vector<const CaloBarShower*> getLocalMaxVCol(std::string name) const;
    std::vector<const LongiCluster*> getLongiClusterUCol(std::string name) const;
    std::vector<const LongiCluster*> getLongiClusterVCol(std::string name) const; 
    std::map<std::string, std::vector<const PandoraPlus::LongiCluster*> > getLongiClusterUMap() const { return map_longiClusUCol; }
    std::map<std::string, std::vector<const PandoraPlus::LongiCluster*> > getLongiClusterVMap() const { return map_longiClusVCol; }

    std::vector<const CaloTower*> getTowers() const {return m_towers}
    std::vector<const CaloUnit*> getBars() const;
    double getEnergy() const; 
  
    void createTower(); 
    void addCluster(const Calo2DCluster* _2dcluster);
    void setLongiClusters( std::string name1, std::vector<const LongiCluster*>& _clU, 
                           std::string name2, std::vector<const LongiCluster*>& _clV )
    { map_longiClusUCol[name1]=_clU; map_longiClusVCol[name2]=_clV;  }

    void setLocalMax( std::string name1, std::vector<const CaloBarShower*>& _colU,
                      std::string name2, std::vector<const CaloBarShower*>& _colV )
    { map_localMaxU[name1]=_colU; map_localMaxV[name2]=_colV; }

    void addLongiClusterU( std::string name, const LongiCluster* _clU ) { map_longiClusUCol[name].push_back(_clU); }
    void addLongiClusterV( std::string name, const LongiCluster* _clV ) { map_longiClusVCol[name].push_back(_clV); }
    void addLocalMaxU( std::string name, const CaloBarShower* _shU ) { map_localMaxU[name].push_back(_shU); }
    void addLocalMaxV( std::string name, const CaloBarShower* _shV ) { map_localMaxV[name].push_back(_shV); }

  private:
    std::vector<int> m_modules;
    std::vector<int> m_parts;
    std::vector<int> m_staves;  

    std::vector<const PandoraPlus::Calo2DCluster*> m_2dclusters;

    std::map<std::string, std::vector<const PandoraPlus::CaloBarShower*> > map_localMaxU;
    std::map<std::string, std::vector<const PandoraPlus::CaloBarShower*> > map_localMaxV;

    std::map<std::string, std::vector<const PandoraPlus::LongiCluster*> > map_longiClusUCol;
    std::map<std::string, std::vector<const PandoraPlus::LongiCluster*> > map_longiClusVCol;

    std::vector<const CaloTower*> m_towers;

    //static const int m_module = 7;
    //static const int m_modulestart = 0;
    //static const int m_part = 4;
    //static const int m_stave = 11;
    //static const int m_superlayer = 14;
    //static const int m_startnumber = 1;
    //static const int m_phibarnumber = 60;
    //static const int m_zbarnumber = 47;
  };

};
#endif
