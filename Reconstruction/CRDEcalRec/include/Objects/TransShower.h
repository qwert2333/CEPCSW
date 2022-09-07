#ifndef TRANSSHOWER_H
#define TRANSSHOWER_H

#include "Objects/Calo1DCluster.h"
#include "Objects/CaloHit.h"
#include <vector>

namespace PandoraPlus{

  class CaloHit; 
  class TransShower{
  public: 	
    TransShower(const Calo1DCluster* _bars1, const Calo1DCluster* _bars2, std::vector<const CaloHit*> _hits)
      : barShowerU(_bars1), 
        barShowerV(_bars2), 
        hits(_hits)
    {};
 
    TransShower() {};
    ~TransShower() { Clear(); };

    void Clear(){
      barShowerU = NULL; barShowerV = NULL; hits.clear();
    }
    void Clean(); 


    TVector3 getPos() const;
    //double getHitsE() const; 
    double getShowerE() const { return barShowerU->getEnergy() + barShowerV->getEnergy();}

    int getModule() const { return module; }
    int getStave()  const { return stave;  }
    int getDlayer() const { return dlayer; }
    int getPart()   const { return part;   }
    double getWidthU() const { return barShowerU->getWidth(); }
    double getWidthV() const { return barShowerV->getWidth(); }
    //double getHitsWidth() const; 
    const Calo1DCluster* getShowerU() const { return barShowerU; }
    const Calo1DCluster* getShowerV() const { return barShowerV; }
    std::vector<const CaloHit*> getCaloHits() const { return hits; }

    void setBarShowers( const Calo1DCluster* _bar1, const Calo1DCluster* _bar2){ barShowerU=_bar1; barShowerV=_bar2; }
    void setCaloHits( std::vector<const CaloHit*> _hits) { hits = _hits; }
    void setIDInfo( int _m, int _s, int _d, int _p ){ module=_m; stave=_s; dlayer=_d; part=_p; } 

  private:
    int module;
    int stave;
    int dlayer;
    int part;
    const Calo1DCluster* barShowerU;
    const Calo1DCluster* barShowerV;
    std::vector<const CaloHit*> hits;
  };

};
#endif
