#ifndef _CRD_CALOBAR_C
#define _CRD_CALOBAR_C

#include "CRDEcalEDMSvc/CRDCaloBar.h"
#include <cmath>

namespace CRDEcalEDM{

  bool CRDCaloBar::isNeighbor(CRDCaloBar &x){
    if(x.getcellID() != cellID &&
       x.getSystem() == system &&
       x.getModele() == module &&
       x.getDlayer() == dlayer &&
       x.getPart()   == part   &&
       x.getBlock()  == block  &&
       x.getSlayer() == slayer &&
       (x.getBar()==bar+1 || x.getBar()==bar-1)
      ) return true;
      return false;
  }

};
#endif
