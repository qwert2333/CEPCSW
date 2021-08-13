#ifndef _CRD_CALOBAR_C
#define _CRD_CALOBAR_C

#include "Objects/CRDCaloBar.h"
#include <cmath>

namespace CRDEcalEDM{

  bool CRDCaloBar::isNeighbor(CRDCaloBar &x){
    if(x.getcellID() != cellID &&
       x.getSystem() == system &&
       x.getModule() == module &&
       x.getDlayer() == dlayer &&
       x.getPart()   == part   &&
       x.getStave()  == stave  &&
       x.getSlayer() == slayer &&
       (x.getBar()==bar+1 || x.getBar()==bar-1)
      ) return true;
      return false;
  }

};
#endif
