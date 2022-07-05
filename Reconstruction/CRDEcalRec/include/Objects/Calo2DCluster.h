#ifndef CALO_2DCLUSTER_H
#define CALO_2DCLUSTER_H

#include "Objects/CaloBar.h"
#include "Objects/Calo1DCluster.h"
#include "Objects/CaloBarShower.h"

namespace PandoraPlus {

  class Calo2DCluster {
  public: 
    Calo2DCluster() {};
    ~Calo2DCluster() { Clear(); };

    void Clear();
    void ClearShower();
    void Clean();
    void Check(); 

    // bool isNeighbor(const PandoraPlus::Calo1DCluster* i1DCluster) const; 

    // int getModule() const { return module; }
    // int getStave()  const { return stave;  }
    // int getDlayer() const { return dlayer; }
    // int getPart()   const { return part;   }
    std::vector<const CaloBar *> getBarXCol() const { return barXCol; }
    std::vector<const CaloBar *> getBarYCol() const { return barYCol; }
    std::vector<const CaloBarShower*> getShowerXCol() const {return barShowerXCol;}
    std::vector<const CaloBarShower*> getShowerYCol() const {return barShowerYCol;}

    // void setIDInfo( int _m, int _s, int _d, int _p ) { module=_m; stave=_s; dlayer=_d; part=_p; }
    void addBar(const CaloBar* _bar) { if(_bar->getSlayer()==0) barXCol.push_back(_bar); if(_bar->getSlayer()==1) barYCol.push_back(_bar); }
    void setShowerXCol(std::vector<const CaloBarShower*> _sh) { barShowerXCol=_sh; }
    void setShowerYCol(std::vector<const CaloBarShower*> _sh) { barShowerYCol=_sh; }
    void addCluster(const Calo1DCluster* _1dcluster){ if(_1dcluster->getCluster().at(0)->getSlayer()==0) m_v1dclusters.push_back(_1dcluster); if(_1dcluster->getCluster().at(0)->getSlayer()==1) m_u1dclusters.push_back(_1dcluster);}
    std::vector<const Calo1DCluster*> getClusterV() const { return m_v1dclusters; }
    std::vector<const Calo1DCluster*> getClusterU() const { return m_u1dclusters; }
    std::vector<const Calo1DCluster*> getCluster() const; 

    std::vector<int> getTowerID() const;
    std::vector<const CaloBar*> getBars() const;
    double getEnergy() const; 

	
  private:
    // int module;
    // int stave;
    // int dlayer;
    // int part;

    std::vector<const CaloBar *> barXCol;  //slayer == 0.
    std::vector<const CaloBar *> barYCol;  //slayer == 1.
    std::vector<const CaloBarShower *> barShowerXCol; 
    std::vector<const CaloBarShower *> barShowerYCol; 
	
	std::vector<const Calo1DCluster*> m_v1dclusters;
	std::vector<const Calo1DCluster*> m_u1dclusters;

  };

};
#endif