#ifndef CALO3D_CLUSTER_C
#define CALO3D_CLUSTER_C

#include "Objects/Calo3DCluster.h"
#include <cmath>

namespace PandoraPlus{

	void Calo3DCluster::Clear() 
	{
		m_2dclusters.clear();
		m_towers.clear();
	}

	void Calo3DCluster::Check()
	{
		for(int i=0; i<m_2dclusters.size(); i++) { delete m_2dclusters[i]; m_2dclusters[i]=NULL; }
		for(int i=0; i<m_towers.size(); i++) { delete m_towers[i]; m_towers[i]=NULL; }
		Clear();
	}

	void Calo3DCluster::Clean()
	{
		for(int i=0; i<m_2dclusters.size(); i++)
		if(!m_2dclusters[i]) { m_2dclusters.erase(m_2dclusters.begin()+i); i--; }
		for(int i=0; i<m_towers.size(); i++)
		if(!m_towers[i]) { m_towers.erase(m_towers.begin()+i); i--; }
	}

	std::vector<int> Calo3DCluster::getTowerID() const
	{
		std::vector<int> towerid;
		towerid.clear();
		int size = 0;

		for(int i=0; i<m_2dclusters.size(); i++)
		{
			const Calo2DCluster* cluster = m_2dclusters.at(i);
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
	
	std::vector<const CaloBar*> Calo3DCluster::getBars() const
	{
		std::vector<const CaloBar*> results;
		results.clear();
		for(int i=0; i<m_2dclusters.size(); i++)
		{
			for(int j=0; j<m_2dclusters.at(i)->getBars().size(); j++)
			{
				results.push_back(m_2dclusters.at(i)->getBars().at(j));
			}
		}
		return results;
	}
	
	double Calo3DCluster::getEnergy() const
	{
		double result = 0;
		for(int m=0; m<m_2dclusters.size(); m++)
		{
			result = result + m_2dclusters.at(m)->getEnergy();
		}
		return result;
	}
};
#endif