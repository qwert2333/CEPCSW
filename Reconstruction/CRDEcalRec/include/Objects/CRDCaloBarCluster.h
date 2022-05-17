#ifndef _CRD_CALOBAR_CLUSTER_
#define _CRD_CALOBAR_CLUSTER_
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>

#include "Objects/CRDCaloBar.h"

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

    double getE() const; 
    TVector3 getPos() const;
    double getScndMoment();
    int getNseeds() const { return Seeds.size(); }
    std::vector<CRDEcalEDM::CRDCaloBar> getBars()  const { return Bars;  }
    std::vector<CRDEcalEDM::CRDCaloBar> getSeeds() const { return Seeds; }
    bool getGlobalRange( double& xmin,  double& ymin, double& zmin, double& xmax, double& ymax, double& zmax ) const;
    int  getLeftEdge();
    int  getRightEdge();

    void PrintBars() const;
    void PrintSeeds() const;

    void addBar( CRDEcalEDM::CRDCaloBar _bar ) { Bars.push_back(_bar); }
    void addSeed( CRDEcalEDM::CRDCaloBar _seed ) { Seeds.push_back(_seed); }
    void setScndMoment(double _scndM ) { ScndMoment = _scndM; }
    void setBars( std::vector<CRDEcalEDM::CRDCaloBar> _bars ) { Bars = _bars; }
    void setSeeds( std::vector<CRDEcalEDM::CRDCaloBar> _seeds) { Seeds = _seeds;  Nseeds = _seeds.size(); }

  private: 
    std::vector<CRDEcalEDM::CRDCaloBar> Bars;
    std::vector<CRDEcalEDM::CRDCaloBar> Seeds;
    double Energy;
    TVector3 pos;
    int Nseeds;
    double ScndMoment;  //Second moment of this cluster. 
  };

}
#endif
