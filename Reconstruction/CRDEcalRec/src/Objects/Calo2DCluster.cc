#ifndef CALO_2DCLUSTER_C
#define CALO_2DCLUSTER_C

#include "Objects/Calo2DCluster.h"
#include <cmath>

namespace PandoraPlus{

  void Calo2DCluster::Clear() {
    m_v1dclusters.clear(); 
	m_u1dclusters.clear(); 
  }

  void Calo2DCluster::ClearShower() {
    barShowerXCol.clear(); 
    barShowerYCol.clear();
  }

  void Calo2DCluster::Check(){
    for(int i=0; i<m_v1dclusters.size(); i++)
      if(!m_v1dclusters[i]) { m_v1dclusters.erase(m_v1dclusters.begin()+i); i--; }
	for(int i=0; i<m_u1dclusters.size(); i++)
      if(!m_u1dclusters[i]) { m_u1dclusters.erase(m_u1dclusters.begin()+i); i--; }
    for(int i=0; i<barXCol.size(); i++)
      if(!barXCol[i]) { barXCol.erase(barXCol.begin()+i); i--; }	  
    for(int i=0; i<barYCol.size(); i++)
      if(!barYCol[i]) { barYCol.erase(barYCol.begin()+i); i--; }
    for(int i=0; i<barShowerXCol.size(); i++)
      if(!barShowerXCol[i]) { barShowerXCol.erase(barShowerXCol.begin()+i); i--; }
    for(int i=0; i<barShowerYCol.size(); i++)
      if(!barShowerYCol[i]) { barShowerYCol.erase(barShowerYCol.begin()+i); i--; }
  	}

  void Calo2DCluster::Clean(){
	for(int i=0; i<m_v1dclusters.size(); i++){ delete m_v1dclusters[i]; m_v1dclusters[i]=NULL; }
	for(int i=0; i<m_u1dclusters.size(); i++){ delete m_u1dclusters[i]; m_u1dclusters[i]=NULL; }
    for(int i=0; i<barXCol.size(); i++) { delete barXCol[i]; barXCol[i]=NULL; }
    for(int i=0; i<barYCol.size(); i++) { delete barYCol[i]; barYCol[i]=NULL; }
    for(int i=0; i<barShowerXCol.size(); i++) { delete barShowerXCol[i]; barShowerXCol[i]=NULL; }
    for(int i=0; i<barShowerYCol.size(); i++) { delete barShowerYCol[i]; barShowerYCol[i]=NULL; }
    Clear();
  	}

	std::vector<const Calo1DCluster*> Calo2DCluster::getCluster() const
	{ 
		std::vector<const Calo1DCluster*> m_1dclusters; 
		m_1dclusters.clear();
		m_1dclusters.insert(m_1dclusters.end(),m_u1dclusters.begin(),m_u1dclusters.end());
		m_1dclusters.insert(m_1dclusters.end(),m_v1dclusters.begin(),m_v1dclusters.end());
		return  m_1dclusters;
	}

	std::vector<int> Calo2DCluster::getTowerID() const
	{
		std::vector<int> towerid;
		towerid.clear();
		int size = 0;
		std::vector<const Calo1DCluster*> m_1dclusters; 
		m_1dclusters.clear();
		m_1dclusters.insert(m_1dclusters.end(),m_u1dclusters.begin(),m_u1dclusters.end());
		m_1dclusters.insert(m_1dclusters.end(),m_v1dclusters.begin(),m_v1dclusters.end());

		for(int i=0; i<m_1dclusters.size(); i++)
		{
			const Calo1DCluster* cluster = m_1dclusters.at(i);
			size = cluster->getTowerID().size();
			for(int j=0; j<size; j++)
			{
				towerid.push_back(cluster->getTowerID().at(j));
			}

		}

		sort(towerid.begin(), towerid.end());
		std::vector<int>::iterator ite = unique(towerid.begin(), towerid.end());
		towerid.erase(ite, towerid.end());

		return towerid;
	}
	
	std::vector<const CaloBar*> Calo2DCluster::getBars() const
	{
		std::vector<const CaloBar*> results;
		results.clear();
		std::vector<const Calo1DCluster*> m_1dclusters; 
		m_1dclusters.clear();
		m_1dclusters.insert(m_1dclusters.end(),m_u1dclusters.begin(),m_u1dclusters.end());
		m_1dclusters.insert(m_1dclusters.end(),m_v1dclusters.begin(),m_v1dclusters.end());
		for(int i=0; i<m_1dclusters.size(); i++)
		{
			for(int j=0; j<m_1dclusters.at(i)->getBars().size(); j++)
			{
				results.push_back(m_1dclusters.at(i)->getBars().at(j));
			}
		}
		return results;
	}
	
	double Calo2DCluster::getEnergy() const 
	{
		double result = 0;
		std::vector<const Calo1DCluster*> m_1dclusters; 
		m_1dclusters.clear();
		m_1dclusters.insert(m_1dclusters.end(),m_u1dclusters.begin(),m_u1dclusters.end());
		m_1dclusters.insert(m_1dclusters.end(),m_v1dclusters.begin(),m_v1dclusters.end());
		for(int m=0; m<m_1dclusters.size(); m++)
		{
			result = result + m_1dclusters.at(m)->getEnergy();
		}
		return result;
	}

};
#endif