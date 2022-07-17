#ifndef CALO_2DCLUSTER_H
#define CALO_2DCLUSTER_H

#include "Objects/CaloUnit.h"
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
      return ( barUCol == x.getBarUCol() &&   barVCol == x.getBarVCol());
    }

    bool isNeighbor(const PandoraPlus::Calo1DCluster* m_1dcluster) const; 

    int getDlayer() const { if(barUCol.size()>0) return barUCol[0]->getDlayer(); else if(barVCol.size()>0) return barVCol[0]->getDlayer(); else return -99; }
    std::vector<int> getModules() const { return m_modules; }
    std::vector<int> getParts() const { return m_parts; }
    std::vector<int> getStaves() const { return m_staves; }

    std::vector<const CaloUnit *> getBarUCol() const { return barUCol; }
    std::vector<const CaloUnit *> getBarVCol() const { return barVCol; }
    std::vector<const CaloBarShower*> getShowerUCol() const {return barShowerUCol;}
    std::vector<const CaloBarShower*> getShowerVCol() const {return barShowerVCol;}

    void addBar(const CaloUnit* _bar) { if(_bar->getSlayer()==0) barUCol.push_back(_bar); if(_bar->getSlayer()==1) barVCol.push_back(_bar); }
    void setShowerUCol(std::vector<const CaloBarShower*> _sh) { barShowerUCol=_sh; }
    void setShowerVCol(std::vector<const CaloBarShower*> _sh) { barShowerVCol=_sh; }
    void addCluster(const Calo1DCluster* _1dcluster);
    std::vector<const Calo1DCluster*> getClusterU() const { return barClusterUCol; }
    std::vector<const Calo1DCluster*> getClusterV() const { return barClusterVCol; }
    std::vector<const Calo1DCluster*> getCluster() const; 

    std::vector<const CaloUnit*> getBars() const;
    double getEnergy() const; 

	
  private:
    std::vector<int> m_modules;
    std::vector<int> m_parts;
    std::vector<int> m_staves;    

    std::vector<const CaloUnit *> barUCol;  //slayer == 0.
    std::vector<const CaloUnit *> barVCol;  //slayer == 1.
    std::vector<const CaloBarShower *> barShowerUCol; 
    std::vector<const CaloBarShower *> barShowerVCol; 
    std::vector<const Calo1DCluster*>  barClusterUCol;
    std::vector<const Calo1DCluster*>  barClusterVCol;

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
