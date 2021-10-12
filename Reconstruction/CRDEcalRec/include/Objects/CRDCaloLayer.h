#ifndef _CRD_CALOLAYER_
#define _CRD_CALOLAYER_

#include <vector>
#include "Objects/CRDCaloBar.h"
#include "Objects/CRDCaloBarShower.h"
#include "Objects/CRDCaloBarCluster.h"

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

    int getModule() const{
      if(barXCol.size()>0) return barXCol[0].getModule();
      if(barYCol.size()>0) return barYCol[0].getModule();
      return -1;
    }

    int getPart() const{
      if(barXCol.size()>0) return barXCol[0].getPart();
      if(barYCol.size()>0) return barYCol[0].getPart();
      return -1;
    }

    int getStave() const{
      if(barXCol.size()>0) return barXCol[0].getStave();
      if(barYCol.size()>0) return barYCol[0].getStave();
      return -1;
    }

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
