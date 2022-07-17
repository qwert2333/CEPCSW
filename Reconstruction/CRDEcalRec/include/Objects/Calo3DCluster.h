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

    inline bool operator == (const Calo3DCluster &x) const{
      return ( m_2dclusters == x.getCluster()  );
    }

    std::vector<int> getModules() const { return m_modules; }
    std::vector<int> getParts() const { return m_parts; }
    std::vector<int> getStaves() const { return m_staves; }

    bool isNeighbor(const PandoraPlus::Calo2DCluster* m_2dcluster) const; 
    void addCluster(const Calo2DCluster* _2dcluster);
    std::vector<const Calo2DCluster*> getCluster() const { return m_2dclusters; }
    
    std::vector<const CaloUnit*> getBars() const;
    double getEnergy() const; 

  private:
    std::vector<int> m_modules;
    std::vector<int> m_parts;
    std::vector<int> m_staves;  

    std::vector<const Calo2DCluster*> m_2dclusters;
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
