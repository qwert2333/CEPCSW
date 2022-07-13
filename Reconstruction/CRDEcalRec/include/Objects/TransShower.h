#ifndef TRANSSHOWER_H
#define TRANSSHOWER_H

#include "Objects/CaloBarShower.h"
#include "Objects/CaloHit.h"
#include <vector>

namespace PandoraPlus{

  class CaloHit; 
  class TransShower{
  public: 	
    TransShower(const CaloBarShower* _bars1, const CaloBarShower* _bars2, std::vector<const CaloHit*> _hits)
      : barShowerX(_bars1), 
        barShowerY(_bars2), 
        hits(_hits)
    {};
 
    TransShower() {};
    ~TransShower() { Clear(); };

    void Clear(){
      barShowerX = NULL; barShowerY = NULL; hits.clear();
    }
    void Clean(); 


    TVector3 getPos() const;
    //double getHitsE() const; 
    double getShowerE() const { return barShowerX->getE() + barShowerY->getE();}

    int getModule() const { return module; }
    int getStave()  const { return stave;  }
    int getDlayer() const { return dlayer; }
    int getPart()   const { return part;   }
    double getWidthX() const { return barShowerX->getWidth(); }
    double getWidthY() const { return barShowerY->getWidth(); }
    //double getHitsWidth() const; 
    const CaloBarShower* getShowerX() const { return barShowerX; }
    const CaloBarShower* getShowerY() const { return barShowerY; }
    std::vector<const CaloHit*> getCaloHits() const { return hits; }

    void setBarShowers( const CaloBarShower* _bar1, const CaloBarShower* _bar2){ barShowerX=_bar1; barShowerY=_bar2; }
    void setCaloHits( std::vector<const CaloHit*> _hits) { hits = _hits; }
    void setIDInfo( int _m, int _s, int _d, int _p ){ module=_m; stave=_s; dlayer=_d; part=_p; } 

  private:
    int module;
    int stave;
    int dlayer;
    int part;
    const CaloBarShower* barShowerX;
    const CaloBarShower* barShowerY;
    std::vector<const CaloHit*> hits;
  };

};
#endif
