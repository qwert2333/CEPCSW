#ifndef CALOBLOCK_C
#define CALOBLOCK_C

#include "Objects/CaloBlock.h"
#include <cmath>

namespace PandoraPlus{

  void CaloBlock::Clear() {
    barXCol.clear();
    barYCol.clear();
    barShowerXCol.clear();
    barShowerYCol.clear();
  }

  void CaloBlock::ClearShower() {
    barShowerXCol.clear(); 
    barShowerYCol.clear();
  }

  void CaloBlock::Check(){
    for(int i=0; i<barXCol.size(); i++)
      if(!barXCol[i]) { barXCol.erase(barXCol.begin()+i); i--; }
    for(int i=0; i<barYCol.size(); i++)
      if(!barYCol[i]) { barYCol.erase(barYCol.begin()+i); i--; }
    for(int i=0; i<barShowerXCol.size(); i++)
      if(!barShowerXCol[i]) { barShowerXCol.erase(barShowerXCol.begin()+i); i--; }
    for(int i=0; i<barShowerYCol.size(); i++)
      if(!barShowerYCol[i]) { barShowerYCol.erase(barShowerYCol.begin()+i); i--; }
  }

  void CaloBlock::Clean(){
    for(int i=0; i<barXCol.size(); i++) { /*delete barXCol[i];*/ barXCol[i]=NULL; }
    for(int i=0; i<barYCol.size(); i++) { /*delete barYCol[i];*/ barYCol[i]=NULL; }
    for(int i=0; i<barShowerXCol.size(); i++) { /*delete barShowerXCol[i];*/ barShowerXCol[i]=NULL; }
    for(int i=0; i<barShowerYCol.size(); i++) { /*delete barShowerYCol[i];*/ barShowerYCol[i]=NULL; }
    Clear();
  }

};
#endif
