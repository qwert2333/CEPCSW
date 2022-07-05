#ifndef CALO_3DCLUSTER_H
#define CALO_3DCLUSTER_H

#include "Objects/Calo2DCluster.h"
#include "Objects/CaloTower.h"

namespace PandoraPlus {

  class Calo3DCluster {
  public: 
    Calo3DCluster() {};
    ~Calo3DCluster() { Clear(); };

    void Clear();
    void Clean();
    void Check(); 
	
	void addCluster(const Calo2DCluster* _2dcluster){m_2dclusters.push_back(_2dcluster);}
	std::vector<const Calo2DCluster*> getCluster() const { return m_2dclusters; }
	std::vector<int> getTowerID() const;
	std::vector<const CaloBar*> getBars() const;

	double getEnergy() const; 

  private:
	std::vector<const Calo2DCluster*> m_2dclusters;
	std::vector<const CaloTower*> m_towers;
  };

};
#endif