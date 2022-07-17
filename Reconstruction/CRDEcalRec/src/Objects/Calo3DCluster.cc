#ifndef CALO_3DCLUSTER_C
#define CALO_3DCLUSTER_C

#include "Objects/Calo3DCluster.h"
#include <cmath>

namespace PandoraPlus{

void Calo3DCluster::Clear() 
{
	m_2dclusters.clear();
	m_towers.clear();
	m_modules.clear();
	m_parts.clear();
	m_staves.clear();
}

void Calo3DCluster::Clean()
{
	for(int i=0; i<m_2dclusters.size(); i++) { delete m_2dclusters[i]; m_2dclusters[i]=NULL; }
	for(int i=0; i<m_towers.size(); i++) { delete m_towers[i]; m_towers[i]=NULL; }
	std::vector<int>().swap(m_modules);
    std::vector<int>().swap(m_parts);
    std::vector<int>().swap(m_staves);
	Clear();
}

void Calo3DCluster::Check()
{
	for(int i=0; i<m_2dclusters.size(); i++)
	if(!m_2dclusters[i]) { m_2dclusters.erase(m_2dclusters.begin()+i); i--; }
	for(int i=0; i<m_towers.size(); i++)
	if(!m_towers[i]) { m_towers.erase(m_towers.begin()+i); i--; }
}

bool Calo3DCluster::isNeighbor(const PandoraPlus::Calo2DCluster* m_2dcluster) const
{
  //Inner module
	for(int i=0; i<m_2dcluster->getModules().size(); i++)
	{
		for(int j=0; j<m_modules.size(); j++)
		{
			if(m_2dcluster->getModules().at(i)==m_modules.at(j)&&m_2dcluster->getParts().at(i)==m_parts.at(j)&&m_2dcluster->getStaves().at(i)==m_staves.at(j))
			{
				return true;
			}
		}
	}
	
  //Over modules
	std::vector<const PandoraPlus::CaloUnit*> bars_2d = m_2dcluster->getBars();
  for(int ib2d=0; ib2d<bars_2d.size(); ib2d++){
    for(int ic=0; ic<m_2dclusters.size(); ic++){
      for(int ib3d=0; ib3d<m_2dclusters[ic]->getBars().size(); ib3d++){
        if(bars_2d[ib2d]->isModuleAdjacent(m_2dclusters[ic]->getBars()[ib3d])) return true;
      }
    }
  }

	return false;
}

void Calo3DCluster::addCluster(const Calo2DCluster* _2dcluster)
{
	m_2dclusters.push_back(_2dcluster);
	std::vector<int> m_2dmodules = _2dcluster->getModules();
	std::vector<int> m_2dparts = _2dcluster->getParts();
	std::vector<int> m_2dstaves = _2dcluster->getStaves();
	m_modules.insert(m_modules.end(),m_2dmodules.begin(),m_2dmodules.end());
	m_parts.insert(m_parts.end(),m_2dparts.begin(),m_2dparts.end());
	m_staves.insert(m_staves.end(),m_2dstaves.begin(),m_2dstaves.end());
}


//

std::vector<const PandoraPlus::CaloUnit*> Calo3DCluster::getBars() const
{
	std::vector<const PandoraPlus::CaloUnit*> results;
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
