#ifndef _CRD_CALOLAYER_
#define _CRD_CALOLAYER_

#include <vector>
#include "CRDEcalEDMSvc/CRDCaloBar.h"
#include "CRDEcalEDMSvc/CRDCaloBarShower.h"
#include "CRDEcalEDMSvc/CRDCaloBarCluster.h"

namespace CRDEcalEDM {

  class CRDCaloLayer{

  public: 
    CRDCaloLayer(){}; 

    void Clear() { 
      ClearBars();
      ClearShowers();
      ClearClusters();
    }
    void ClearBars() { barXCol.clear(); barYCol.clear(); }
    void ClearShowers() { barShowerXCol.clear(); barShowerYCol.clear(); }
    void ClearClusters() { barClusXCol.clear(); barClusYCol.clear(); }

    int getDlayer() const {
      if(barXCol.size()>0) return barXCol[0].getDlayer();
      if(barYCol.size()>0) return barYCol[0].getDlayer();
      return -1;
    }

    std::vector<CRDEcalEDM::CRDCaloBar> barXCol;
    std::vector<CRDEcalEDM::CRDCaloBar> barYCol;
    std::vector<CRDEcalEDM::CRDCaloBarShower> barShowerXCol;
    std::vector<CRDEcalEDM::CRDCaloBarShower> barShowerYCol;
    std::vector<CRDEcalEDM::CRDCaloBarCluster> barClusXCol;
    std::vector<CRDEcalEDM::CRDCaloBarCluster> barClusYCol;

  };
};
#endif
