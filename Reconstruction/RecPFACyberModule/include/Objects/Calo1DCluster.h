#ifndef CALO_1DCLUSTER_H
#define CALO_1DCLUSTER_H
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>

#include "Objects/CaloUnit.h"

namespace Cyber{
  class Calo1DCluster{
  public: 

    Calo1DCluster() {}; 
    Calo1DCluster( std::vector<const Cyber::CaloUnit*> _bars, std::vector<const Cyber::CaloUnit*> _seeds)
    : Bars(_bars), Seeds(_seeds) {};
    ~Calo1DCluster() { Clear(); }

    void Clear();
    void Clean();
    void Check();
    std::shared_ptr<Cyber::Calo1DCluster> Clone() const; 

    inline bool operator == (const Calo1DCluster &x) const{
      return ( Bars == x.getBars()  );
    }

    bool isNeighbor(const Cyber::CaloUnit* m_bar) const; 
    bool inCluster(const Cyber::CaloUnit* iBar) const;
    void sortByPos() { std::sort(Bars.begin(), Bars.end(), compPos); }

    double getEnergy() const; 
    TVector3 getPos() const;
    double getT1() const;
    double getT2() const;
    double getWidth() const;
    double getScndMoment() const;

    int getNseeds() const { return Seeds.size(); }
    std::vector<const Cyber::CaloUnit*> getBars()  const { return Bars;  }
	  std::vector<const Cyber::CaloUnit*> getCluster()  const { return Bars;  }
    std::vector<const Cyber::CaloUnit*> getSeeds() const { return Seeds; }
    std::vector< const Cyber::Calo1DCluster* > getCousinClusters() const { return CousinClusters; }
    std::vector< const Cyber::Calo1DCluster* > getChildClusters() const { return ChildClusters; }
    bool getGlobalRange( double& xmin,  double& ymin, double& zmin, double& xmax, double& ymax, double& zmax ) const;
    int  getLeftEdge();
    int  getRightEdge();
    std::vector< std::pair<edm4hep::MCParticle, float> > getLinkedMCP() const { return MCParticleWeight; }
    std::vector< std::pair<edm4hep::MCParticle, float> > getLinkedMCPfromUnit();
    edm4hep::MCParticle getLeadingMCP() const;
    float getLeadingMCPweight() const;

	  void addUnit(const Cyber::CaloUnit* _bar );
    void addSeed(const Cyber::CaloUnit* _seed ) { Seeds.push_back(_seed); }
    void setBars( std::vector<const Cyber::CaloUnit*> _bars ) { Bars = _bars; }
    void setSeeds( std::vector<const Cyber::CaloUnit*> _seeds) { Seeds = _seeds; }
    void addCousinCluster( const Cyber::Calo1DCluster* clus ) { CousinClusters.push_back(clus); }
    void addChildCluster( const Cyber::Calo1DCluster* clus ) { ChildClusters.push_back(clus); }
    void deleteCousinCluster( const Cyber::Calo1DCluster* _cl );
    void addLinkedMCP( std::pair<edm4hep::MCParticle, float> _pair ) { MCParticleWeight.push_back(_pair); }
    void setLinkedMCP( std::vector<std::pair<edm4hep::MCParticle, float>> _pairVec ) { MCParticleWeight.clear(); MCParticleWeight = _pairVec; }
    void setSeed();  //Set the most energitic unit as seed, Eseed>5 MeV (hardcoded). 
    void setIDInfo(); 

    int getDlayer() const { if(Bars.size()>0) return Bars[0]->getDlayer(); return -99;  }
    int getSlayer() const { if(Bars.size()>0) return Bars[0]->getSlayer(); return -99;  }
    int getSystem() const { if(Bars.size()>0) return Bars[0]->getSystem(); return -99;  }
    std::vector< std::vector<int> > getTowerID() const { return towerID; }
    
  private: 
    std::vector<const Cyber::CaloUnit*> Bars;
    std::vector<const Cyber::CaloUnit*> Seeds;
    double Energy;
    TVector3 pos;

    std::vector< std::vector<int> > towerID; //[system, module, stave, part]

    std::vector< const Cyber::Calo1DCluster* > CousinClusters;
    std::vector< const Cyber::Calo1DCluster* > ChildClusters;

    std::vector< std::pair<edm4hep::MCParticle, float> > MCParticleWeight;

    static bool compPos( const Cyber::CaloUnit* hit1, const Cyber::CaloUnit* hit2 )
      { return *hit1 < *hit2; }
  };
}
#endif
