#ifndef CALOBARSHOWER_H
#define CALOBARSHOWER_H

#include <vector>
#include "Objects/CaloBar.h"

namespace PandoraPlus {

  class CaloBarShower{

  public: 
    CaloBarShower() {}; 
    ~CaloBarShower() { Clean(); }

    inline bool operator == (const CaloBarShower &x) const{
      return Bars == x.Bars;
    }

    void Clear();
    void Clean();
    void Check();

    //bool isNeighbor(PandoraPlus::CaloBar* iBar); 
    //bool inShower(PandoraPlus::CaloBar* iBar); 

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
    std::vector<const PandoraPlus::CaloBar*> getBars() const { return Bars; }
    const PandoraPlus::CaloBar* getSeed() const { return Seed; }

    void setIDInfo( int _m, int _s, int _d, int _p, int _sl ){ module=_m; stave=_s; dlayer=_d; part=_p; slayer=_sl; }
    void setIDInfo(); 
    void setBars(std::vector<const PandoraPlus::CaloBar*> _bars){ Bars = _bars; }
    void addBar(const PandoraPlus::CaloBar* _bar) { Bars.push_back(_bar); }
    void setSeed(const PandoraPlus::CaloBar* _seed ) { Seed = _seed; }

  private:
    int module;
    int stave;
    int part;
    int dlayer;
    int slayer;

    std::vector<const PandoraPlus::CaloBar*> Bars;
    const PandoraPlus::CaloBar* Seed;
  };


};
#endif
