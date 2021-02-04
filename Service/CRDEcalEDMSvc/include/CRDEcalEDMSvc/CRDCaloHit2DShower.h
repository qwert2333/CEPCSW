#ifndef _CRD_CALOHIT_2DSHOWER_
#define _CRD_CALOHIT_2DSHOWER_
#include <vector>
#include "CRDEcalEDMSvc/CRDCaloBarShower.h"
#include "edm4hep/CalorimeterHit.h"

namespace CRDEcalDigiEDM {

  class CRDCaloHit2DShower{
  public: 	
    void Clear(){
      barShowerX.Clear(); barShowerY.Clear(); CaloHits.clear();
    }
    CRDEcalDigiEDM::CRDCaloBarShower barShowerX;
    CRDEcalDigiEDM::CRDCaloBarShower barShowerY;
    std::vector<edm4hep::ConstCalorimeterHit> CaloHits;
  };

};
#endif
