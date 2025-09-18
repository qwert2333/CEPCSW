#ifndef CALO_2DCLUSTER_H
#define CALO_2DCLUSTER_H

#include "Objects/CaloHit.h"
#include "Objects/CaloUnit.h"
#include "Objects/Calo1DCluster.h"

namespace Cyber {

  class CaloHit;
  class Calo2DCluster {
  public: 
    Calo2DCluster() { pos.SetXYZ(0., 0., 0.); };
    ~Calo2DCluster() { Clear(); };

    void Clear();
    void ClearShower();
    void Clean();
    void Check(); 

    inline bool operator == (const Calo2DCluster &x) const{
      return ( barUCol == x.getBarUCol() &&   barVCol == x.getBarVCol());
    }

    bool isNeighbor(const Cyber::Calo1DCluster* m_1dcluster) const; 

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
    void addTowerID(int _sys, int _m,  int _s, int _p) { std::vector<int> id(4); id[0] = _sys; id[1] = _m; id[2] = _s; id[3] = _p; towerID.push_back(id); }
    void addTowerID(std::vector<int> id) { towerID.push_back(id); }
    void setTowerID(std::vector<int> id) { towerID.clear(); towerID.push_back(id); }
    void setPos( TVector3 _vec ) { pos = _vec; }
    void setPos( double _x, double _y, double _z ) { pos.SetXYZ(_x, _y, _z); }

  private:
    std::vector< std::vector<int> > towerID; //[system, module, stave, part]

    TVector3 pos = TVector3(0.,0.,0.);
    std::vector<const CaloHit*> hits; 
    std::vector<const CaloUnit*> barUCol;  //slayer == 0.
    std::vector<const CaloUnit*> barVCol;  //slayer == 1.
    std::vector<const Calo1DCluster*> barShowerUCol; 
    std::vector<const Calo1DCluster*> barShowerVCol; 

  };

};
#endif
