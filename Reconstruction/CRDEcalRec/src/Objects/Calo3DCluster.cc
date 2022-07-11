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
	
	std::vector<const PandoraPlus::CaloBar*> bars_2d = m_2dcluster->getBars();
	std::vector<const PandoraPlus::CaloBar*> bars_3d;
	bars_3d.clear();
	for(int i=0; i<m_2dclusters.size(); i++)
	{
		for(int j=0; j<m_2dclusters.at(i)->getBars().size(); j++)
		{
			bars_3d.push_back(m_2dclusters.at(i)->getBars().at(j));
		}
	}

	for(int i=0; i<bars_2d.size(); i++)
	{
		const PandoraPlus::CaloBar* bob = bars_2d.at(i);
		
		for(int j=0; j<bars_3d.size(); j++)
		{
			const PandoraPlus::CaloBar* alice = bars_3d.at(j);
			// cout<<i<<j<<endl;
			if(ifModuleAdjacent(bob,alice))
			{
				// cout<<"ifModuleAdjacent"<<endl;
				return true;
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

bool Calo3DCluster::ifModuleAdjacent(const PandoraPlus::CaloBar* bar_2d, const PandoraPlus::CaloBar* bar_3d) const
{
	PandoraPlus::CaloBar bob; //just first layer
	PandoraPlus::CaloBar alice; //second new
	if(bar_2d->getModule()==m_module && bar_3d->getModule()==m_modulestart)
	{
		bob = *bar_2d;
		alice = *bar_3d;
		return ifAdjacent(bob, alice);	
	}
	else if(bar_2d->getModule()==m_modulestart && bar_3d->getModule()==m_module)
	{
		bob = *bar_3d;
		alice = *bar_2d;
		return ifAdjacent(bob, alice);
	}
	else if(bar_2d->getModule() == bar_3d->getModule())
	{
		return false;
	}
	else if((bar_2d->getModule() < bar_3d->getModule()) && (bar_2d->getModule() - bar_3d->getModule())==-1)
	{
		bob = *bar_2d;
		alice = *bar_3d;
		return ifAdjacent(bob, alice);
	}
	else if((bar_2d->getModule() > bar_3d->getModule()) && (bar_2d->getModule() - bar_3d->getModule())==1)
	{
		bob = *bar_3d;
		alice = *bar_2d;
		return ifAdjacent(bob, alice);
	}	
	else
	{
		return false;
	}	
}

//
bool Calo3DCluster::ifAdjacent(PandoraPlus::CaloBar &bob, PandoraPlus::CaloBar &alice) const
{
	if( (bob.getDlayer()==1 && bob.getSlayer()==0 && bob.getPart()==m_part && alice.getSlayer()==0 && alice.getPart()==1 && bob.getStave()==alice.getStave() && abs(bob.getBar()-alice.getBar())<=1 ) ||
		(bob.getDlayer()==1 && bob.getSlayer()==0 && bob.getPart()==m_part && alice.getSlayer()==1 && alice.getPart()==1 && bob.getStave()==alice.getStave() && alice.getBar()==1) ||
		(bob.getDlayer()==1 && bob.getSlayer()==0 && bob.getPart()==m_part && alice.getSlayer()==0 && alice.getPart()==1 && abs(bob.getStave()-alice.getStave())<=1 && ((bob.getBar()==1&&alice.getBar()==60)||(bob.getBar()==60&&alice.getBar()==1)))
	) 
	{
		return true;
	}
	else
	{
		return false;
	}
}

std::vector<const PandoraPlus::CaloBar*> Calo3DCluster::getBars() const
{
	std::vector<const PandoraPlus::CaloBar*> results;
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