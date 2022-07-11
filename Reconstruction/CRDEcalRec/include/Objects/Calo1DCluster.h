#ifndef CALO_1DCLUSTER_H
#define CALO_1DCLUSTER_H
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>

#include "Objects/CaloBar.h"

namespace PandoraPlus{
  class Calo1DCluster{
  public: 

    Calo1DCluster() {}; 
    Calo1DCluster( std::vector<const PandoraPlus::CaloBar*> _bars, std::vector<const PandoraPlus::CaloBar*> _seeds)
    : Bars(_bars), Seeds(_seeds) {};
    ~Calo1DCluster() { Clear(); }

    void Clear();
    void Clean();
    void Check();

    inline bool operator == (const Calo1DCluster &x) const{
      return ( Bars == x.getBars()  );
    }

    bool isNeighbor(const PandoraPlus::CaloBar* m_bar) const; 
    bool inCluster(const PandoraPlus::CaloBar* iBar) const;
    void sortByPos() { std::sort(Bars.begin(), Bars.end()); }

    double getEnergy() const; 
    TVector3 getPos() const;
    double getScndMoment() const;
    int getNseeds() const { return Seeds.size(); }
    std::vector<const PandoraPlus::CaloBar*> getBars()  const { return Bars;  }
	  std::vector<const PandoraPlus::CaloBar*> getCluster()  const { return Bars;  }
    std::vector<const PandoraPlus::CaloBar*> getSeeds() const { return Seeds; }
    bool getGlobalRange( double& xmin,  double& ymin, double& zmin, double& xmax, double& ymax, double& zmax ) const;
    int  getLeftEdge();
    int  getRightEdge();

    void PrintBars() const;
    void PrintSeeds() const;
	
	  void addCluster(const PandoraPlus::CaloBar* _bar );
    void addBar(const PandoraPlus::CaloBar* _bar ) { Bars.push_back(_bar); }
    void addSeed(const PandoraPlus::CaloBar* _seed ) { Seeds.push_back(_seed); }
    void setBars( std::vector<const PandoraPlus::CaloBar*> _bars ) { Bars = _bars; }
    void setSeeds( std::vector<const PandoraPlus::CaloBar*> _seeds) { Seeds = _seeds; }

    std::vector<int> getModules() const { return m_modules; }
    std::vector<int> getParts() const { return m_parts; }
    std::vector<int> getStaves() const { return m_staves; }
    
  private: 
    std::vector<const PandoraPlus::CaloBar*> Bars;
    std::vector<const PandoraPlus::CaloBar*> Seeds;
    double Energy;
    TVector3 pos;
    std::vector<int> m_modules;
    std::vector<int> m_parts;
    std::vector<int> m_staves;

    static const int m_module = 7;
    static const int m_modulestart = 0;
    static const int m_part = 4;
    static const int m_stave = 11;
    static const int m_superlayer = 14;
    static const int m_startnumber = 1;
    static const int m_phibarnumber = 60;
    static const int m_zbarnumber = 47;
  };
}
#endif