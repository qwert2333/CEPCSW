#ifndef CALO_2DCLUSTER_H
#define CALO_2DCLUSTER_H

#include "Objects/CaloBar.h"
#include "Objects/Calo1DCluster.h"
#include "Objects/CaloBarShower.h"

namespace PandoraPlus {

  class Calo2DCluster {
  public: 
    Calo2DCluster() {};
    ~Calo2DCluster() { Clear(); };

    void Clear();
    void ClearShower();
    void Clean();
    void Check(); 

    inline bool operator == (const Calo2DCluster &x) const{
      return ( m_v1dclusters == x.getClusterV() &&   m_u1dclusters == x.getClusterU());
    }

    bool isNeighbor(const PandoraPlus::Calo1DCluster* m_1dcluster) const; 

    std::vector<int> getModules() const { return m_modules; }
    std::vector<int> getParts() const { return m_parts; }
    std::vector<int> getStaves() const { return m_staves; }

    std::vector<const CaloBar *> getBarXCol() const { return barXCol; }
    std::vector<const CaloBar *> getBarYCol() const { return barYCol; }
    std::vector<const CaloBarShower*> getShowerXCol() const {return barShowerXCol;}
    std::vector<const CaloBarShower*> getShowerYCol() const {return barShowerYCol;}

    void addBar(const CaloBar* _bar) { if(_bar->getSlayer()==0) barXCol.push_back(_bar); if(_bar->getSlayer()==1) barYCol.push_back(_bar); }
    void setShowerXCol(std::vector<const CaloBarShower*> _sh) { barShowerXCol=_sh; }
    void setShowerYCol(std::vector<const CaloBarShower*> _sh) { barShowerYCol=_sh; }
    void addCluster(const Calo1DCluster* _1dcluster);
    std::vector<const Calo1DCluster*> getClusterV() const { return m_v1dclusters; }
    std::vector<const Calo1DCluster*> getClusterU() const { return m_u1dclusters; }
    std::vector<const Calo1DCluster*> getCluster() const; 

    std::vector<const CaloBar*> getBars() const;
    double getEnergy() const; 

	
  private:
    std::vector<int> m_modules;
    std::vector<int> m_parts;
    std::vector<int> m_staves;    

    std::vector<const CaloBar *> barXCol;  //slayer == 0.
    std::vector<const CaloBar *> barYCol;  //slayer == 1.
    std::vector<const CaloBarShower *> barShowerXCol; 
    std::vector<const CaloBarShower *> barShowerYCol; 
    std::vector<const Calo1DCluster*> m_v1dclusters;
    std::vector<const Calo1DCluster*> m_u1dclusters;

    static const int m_module = 7;
    static const int m_modulestart = 0;
    static const int m_part = 4;
    static const int m_stave = 11;
    static const int m_superlayer = 14;
    static const int m_startnumber = 1;
    static const int m_phibarnumber = 60;
    static const int m_zbarnumber = 47;
  };

};
#endif