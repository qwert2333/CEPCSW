#ifndef CALOBLOCK_H
#define CALOBLOCK_H

#include "Objects/CaloUnit.h"
#include "Objects/CaloBarShower.h"
#include "Objects/Track.h"
namespace PandoraPlus {

  class CaloBlock {
  public: 
    CaloBlock(int _module, int _stave, int _dlayer, int _part): module(_module), stave(_stave), dlayer(_dlayer), part(_part) {};
    CaloBlock() {};
    ~CaloBlock() { Clear(); };

    void Clear();
    void ClearShower();
    void Clean();
    void Check(); 


    int getModule() const { return module; }
    int getStave()  const { return stave;  }
    int getDlayer() const { return dlayer; }
    int getPart()   const { return part;   }
    std::vector<const CaloUnit *> getBarXCol() const { return barXCol; }
    std::vector<const CaloUnit *> getBarYCol() const { return barYCol; }
    std::vector<const CaloBarShower*> getShowerXCol() const {return barShowerXCol;}
    std::vector<const CaloBarShower*> getShowerYCol() const {return barShowerYCol;}

    void setIDInfo( int _m, int _s, int _d, int _p ) { module=_m; stave=_s; dlayer=_d; part=_p; }
    void addBar(const CaloUnit* _bar) { if(_bar->getSlayer()==0) barXCol.push_back(_bar); if(_bar->getSlayer()==1) barYCol.push_back(_bar); }
    void setShowerXCol(std::vector<const CaloBarShower*> _sh) { barShowerXCol=_sh; }
    void setShowerYCol(std::vector<const CaloBarShower*> _sh) { barShowerYCol=_sh; }

  private:
    int module;
    int stave;
    int dlayer;
    int part;

    std::vector<const CaloUnit *> barXCol;  //slayer == 0.
    std::vector<const CaloUnit *> barYCol;  //slayer == 1.
    std::vector<const CaloBarShower *> barShowerXCol; 
    std::vector<const CaloBarShower *> barShowerYCol; 

  };

};
#endif
