#ifndef CALO_2DCLUSTER_H
#define CALO_2DCLUSTER_H

#include "Objects/CaloHit.h"
#include "Objects/CaloUnit.h"
#include "Objects/Calo1DCluster.h"

namespace PandoraPlus {

  class CaloHit;
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
    std::vector< std::vector<int> > getTowerID() const { return towerID; }

    std::vector<const CaloUnit*> getBars() const;
    std::vector<const CaloHit*> getCaloHits() const { return hits; }
    std::vector<const CaloUnit *> getBarUCol() const { return barUCol; }
    std::vector<const CaloUnit *> getBarVCol() const { return barVCol; }
    std::vector<const Calo1DCluster*> getShowerUCol() const {return barShowerUCol;}
    std::vector<const Calo1DCluster*> getShowerVCol() const {return barShowerVCol;}
    std::vector<const Calo1DCluster*> getCluster() const;
    double getEnergy() const; 
    TVector3 getPos() const; 

    void setCaloHits( std::vector<const CaloHit*> _hits) { hits = _hits; }
    void addBar(const CaloUnit* _bar) { if(_bar->getSlayer()==0) barUCol.push_back(_bar); if(_bar->getSlayer()==1) barVCol.push_back(_bar); }
    void setBarUCol( std::vector<const CaloUnit*> _bars ) { barUCol=_bars; }
    void setBarVCol( std::vector<const CaloUnit*> _bars ) { barVCol=_bars; }
    void addShowerU( const Calo1DCluster* _sh) { barShowerUCol.push_back(_sh); }
    void addShowerV( const Calo1DCluster* _sh) { barShowerVCol.push_back(_sh); }
    void setShowerUCol(std::vector<const Calo1DCluster*> _sh) { barShowerUCol=_sh; }
    void setShowerVCol(std::vector<const Calo1DCluster*> _sh) { barShowerVCol=_sh; }
    void addUnit(const Calo1DCluster* _1dcluster);
    void addTowerID(int _m, int _p, int _s) { std::vector<int> id(3); id[0] = _m; id[1] = _p; id[2] = _s; towerID.push_back(id); }
    void addTowerID(std::vector<int> id) { towerID.push_back(id); }
	
  private:
    std::vector< std::vector<int> > towerID; //[module, part, stave]

    std::vector<const CaloHit*> hits; 
    std::vector<const CaloUnit*> barUCol;  //slayer == 0.
    std::vector<const CaloUnit*> barVCol;  //slayer == 1.
    std::vector<const Calo1DCluster*> barShowerUCol; 
    std::vector<const Calo1DCluster*> barShowerVCol; 

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
