#ifndef CALO_1DCLUSTER_H
#define CALO_1DCLUSTER_H
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>

#include "Objects/CaloUnit.h"

namespace PandoraPlus{
  class Calo1DCluster{
  public: 

    Calo1DCluster() {}; 
    Calo1DCluster( std::vector<const PandoraPlus::CaloUnit*> _bars, std::vector<const PandoraPlus::CaloUnit*> _seeds)
    : Bars(_bars), Seeds(_seeds) {};
    ~Calo1DCluster() { Clear(); }

    void Clear();
    void Clean();
    void Check();
    std::shared_ptr<PandoraPlus::Calo1DCluster> Clone() const; 

    inline bool operator == (const Calo1DCluster &x) const{
      return ( Bars == x.getBars()  );
    }

    bool isNeighbor(const PandoraPlus::CaloUnit* m_bar) const; 
    bool inCluster(const PandoraPlus::CaloUnit* iBar) const;
    void sortByPos() { std::sort(Bars.begin(), Bars.end()); }

    double getEnergy() const; 
    TVector3 getPos() const;
    double getT1() const;
    double getT2() const;
    double getWidth() const;
    double getScndMoment() const;

    int getNseeds() const { return Seeds.size(); }
    std::vector<const PandoraPlus::CaloUnit*> getBars()  const { return Bars;  }
	  std::vector<const PandoraPlus::CaloUnit*> getCluster()  const { return Bars;  }
    std::vector<const PandoraPlus::CaloUnit*> getSeeds() const { return Seeds; }
    std::vector< const PandoraPlus::Calo1DCluster* > getCousinClusters() const { return CousinClusters; }
    std::vector< const PandoraPlus::Calo1DCluster* > getChildClusters() const { return ChildClusters; }
    bool getGlobalRange( double& xmin,  double& ymin, double& zmin, double& xmax, double& ymax, double& zmax ) const;
    int  getLeftEdge();
    int  getRightEdge();

	  void addUnit(const PandoraPlus::CaloUnit* _bar );
    void addSeed(const PandoraPlus::CaloUnit* _seed ) { Seeds.push_back(_seed); }
    void setBars( std::vector<const PandoraPlus::CaloUnit*> _bars ) { Bars = _bars; }
    void setSeeds( std::vector<const PandoraPlus::CaloUnit*> _seeds) { Seeds = _seeds; }
    void addCousinCluster( const PandoraPlus::Calo1DCluster* clus ) { CousinClusters.push_back(clus); }
    void addChildCluster( const PandoraPlus::Calo1DCluster* clus ) { ChildClusters.push_back(clus); }
    void deleteCousinCluster( const PandoraPlus::Calo1DCluster* _cl );
    void setSeed();  //Set the most energitic unit as seed, Eseed>5 MeV (hardcoded). 
    void setIDInfo(); 

    int getDlayer() const { if(Bars.size()>0) return Bars[0]->getDlayer(); return -99;  }
    int getSlayer() const { if(Bars.size()>0) return Bars[0]->getSlayer(); return -99;  }
    std::vector< std::vector<int> > getTowerID() const { return towerID; }
    
  private: 
    std::vector<const PandoraPlus::CaloUnit*> Bars;
    std::vector<const PandoraPlus::CaloUnit*> Seeds;
    double Energy;
    TVector3 pos;

    std::vector< std::vector<int> > towerID; //[module, part, stave]

    std::vector< const PandoraPlus::Calo1DCluster* > CousinClusters;
    std::vector< const PandoraPlus::Calo1DCluster* > ChildClusters;

  };
}
#endif
