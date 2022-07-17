#ifndef CALOBARSHOWER_H
#define CALOBARSHOWER_H

#include <vector>
#include "Objects/CaloUnit.h"

namespace PandoraPlus {

  class CaloBarShower{

  public: 
    CaloBarShower() {}; 
    ~CaloBarShower() { Clear(); }

    inline bool operator == (const CaloBarShower &x) const{
      return Bars == x.Bars;
    }

    void Clear();
    void Clean();
    void Check();

    //bool isNeighbor(PandoraPlus::CaloUnit* iBar); 
    //bool inShower(PandoraPlus::CaloUnit* iBar); 

    TVector3 getPos() const;
    double getE()  const;
    double getT1() const;
    double getT2() const;
    double getWidth() const; 
    int getModule() const { return module; }
    int getStave()  const { return stave;  }
    int getDlayer() const { return dlayer; }
    int getPart()   const { return part;   }
    int getSlayer() const { return slayer; }
    std::vector<const PandoraPlus::CaloUnit*> getBars() const { return Bars; }
    const PandoraPlus::CaloUnit* getSeed() const { return Seed; }

    void setIDInfo( int _m, int _s, int _d, int _p, int _sl ){ module=_m; stave=_s; dlayer=_d; part=_p; slayer=_sl; }
    void setIDInfo(); 
    void setBars(std::vector<const PandoraPlus::CaloUnit*> _bars){ Bars = _bars; }
    void addBar(const PandoraPlus::CaloUnit* _bar) { Bars.push_back(_bar); }
    void setSeed(const PandoraPlus::CaloUnit* _seed ) { Seed = _seed; }

  private:
    int module;
    int stave;
    int part;
    int dlayer;
    int slayer;

    std::vector<const PandoraPlus::CaloUnit*> Bars;
    const PandoraPlus::CaloUnit* Seed;
  };


};
#endif
