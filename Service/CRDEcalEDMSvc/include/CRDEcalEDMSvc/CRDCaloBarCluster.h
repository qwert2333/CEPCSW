#ifndef _CRD_CALOBAR_CLUSTER_
#define _CRD_CALOBAR_CLUSTER_
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>

#include "CRDEcalEDMSvc/CRDCaloBar.h"

namespace CRDEcalEDM{
  class CRDCaloBarCluster{
  public: 

    CRDCaloBarCluster() {}; 
    CRDCaloBarCluster( std::vector<CRDEcalEDM::CRDCaloBar> _bars, std::vector<CRDEcalEDM::CRDCaloBar> _seeds)
    : Bars(_bars), Seeds(_seeds) {};

    void Clear();
    bool isNeighbor(CRDEcalEDM::CRDCaloBar iBar); 
    bool inCluster(CRDEcalEDM::CRDCaloBar iBar);
    void sortByPos() { std::sort(Bars.begin(), Bars.end()); }

    double getE(); 
    dd4hep::Position getPos();
    double getScndMoment();
    std::vector<CRDEcalEDM::CRDCaloBar> getBars()  { return Bars;  }
    std::vector<CRDEcalEDM::CRDCaloBar> getSeeds() { return Seeds; }

    void PrintBars();
    void PrintSeeds();

    void setBars( std::vector<CRDEcalEDM::CRDCaloBar> _bars ) { Bars = _bars; }
    void setSeeds( std::vector<CRDEcalEDM::CRDCaloBar> _seeds) { Seeds = _seeds; }

  private: 
    std::vector<CRDEcalEDM::CRDCaloBar> Bars;
    std::vector<CRDEcalEDM::CRDCaloBar> Seeds;
    double Energy;
    dd4hep::Position pos;
    int Nseeds;
    double ScndMoment;  //Second moment of this cluster. 
  };

}
#endif
