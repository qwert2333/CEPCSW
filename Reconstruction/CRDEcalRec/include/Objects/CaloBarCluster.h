#ifndef CALOBAR_CLUSTER_H
#define CALOBAR_CLUSTER_H
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>

#include "Objects/CaloUnit.h"

namespace PandoraPlus{
  class CaloBarCluster{
  public: 

    CaloBarCluster() {}; 
    CaloBarCluster( std::vector<const PandoraPlus::CaloUnit*> _bars, std::vector<const PandoraPlus::CaloUnit*> _seeds)
    : Bars(_bars), Seeds(_seeds) {};
    ~CaloBarCluster() { Clear(); }

    void Clear();
    void Clean();
    void Check();

    bool isNeighbor(const PandoraPlus::CaloUnit* iBar) const; 
    bool inCluster(const PandoraPlus::CaloUnit* iBar) const;
    void sortByPos() { std::sort(Bars.begin(), Bars.end()); }

    double getE() const; 
    TVector3 getPos() const;
    double getScndMoment() const;
    int getNseeds() const { return Seeds.size(); }
    std::vector<const PandoraPlus::CaloUnit*> getBars()  const { return Bars;  }
    std::vector<const PandoraPlus::CaloUnit*> getSeeds() const { return Seeds; }
    bool getGlobalRange( double& xmin,  double& ymin, double& zmin, double& xmax, double& ymax, double& zmax ) const;
    int  getLeftEdge();
    int  getRightEdge();

    void PrintBars() const;
    void PrintSeeds() const;

    void addBar(const PandoraPlus::CaloUnit* _bar ) { Bars.push_back(_bar); }
    void addSeed(const PandoraPlus::CaloUnit* _seed ) { Seeds.push_back(_seed); }
    void setBars( std::vector<const PandoraPlus::CaloUnit*> _bars ) { Bars = _bars; }
    void setSeeds( std::vector<const PandoraPlus::CaloUnit*> _seeds) { Seeds = _seeds; }

  private: 
    std::vector<const PandoraPlus::CaloUnit*> Bars;
    std::vector<const PandoraPlus::CaloUnit*> Seeds;
    double Energy;
    TVector3 pos;
  };

}
#endif
