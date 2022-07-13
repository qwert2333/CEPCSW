#ifndef CALO_2DCLUSTER_C
#define CALO_2DCLUSTER_C

#include "Objects/Calo2DCluster.h"
#include <cmath>

namespace PandoraPlus{

  void Calo2DCluster::Clear() {
	m_modules.clear(); 
    m_parts.clear(); 
    m_staves.clear();
    barXCol.clear();
    barYCol.clear();
    barShowerXCol.clear(); 
    barShowerYCol.clear();  
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
	std::vector<int>().swap(m_modules);
    std::vector<int>().swap(m_parts);
    std::vector<int>().swap(m_staves);
    Clear();
  	}

  bool Calo2DCluster::isNeighbor(const PandoraPlus::Calo1DCluster* m_1dcluster) const
  {
	assert(m_1dcluster->getBars().size() > 0 && getCluster().at(0)->getBars().size()>0 );
	if(m_1dcluster->getBars().at(0)->getDlayer() == getCluster().at(0)->getBars().at(0)->getDlayer()  )
	{
		for(int i=0; i<m_1dcluster->getModules().size(); i++)
		{
			for(int j=0; j<m_modules.size(); j++)
			{
				if(m_1dcluster->getModules().at(i)==m_modules.at(j)&&m_1dcluster->getParts().at(i)==m_parts.at(j)&&m_1dcluster->getStaves().at(i)==m_staves.at(j))
				{
					return true;
				}
			}
		}
	}
	return false;
  } 
  void Calo2DCluster::addCluster(const Calo1DCluster* _1dcluster)
  {
	assert(_1dcluster->getCluster().size()>0);
	if(_1dcluster->getCluster().at(0)->getSlayer()==0) m_v1dclusters.push_back(_1dcluster); 
	if(_1dcluster->getCluster().at(0)->getSlayer()==1) m_u1dclusters.push_back(_1dcluster);
	std::vector<int> m_1dmodules = _1dcluster->getModules();
	std::vector<int> m_1dparts = _1dcluster->getParts();
	std::vector<int> m_1dstaves = _1dcluster->getStaves();
	m_modules.insert(m_modules.end(),m_1dmodules.begin(),m_1dmodules.end());
	m_parts.insert(m_parts.end(),m_1dparts.begin(),m_1dparts.end());
	m_staves.insert(m_staves.end(),m_1dstaves.begin(),m_1dstaves.end());
  }

  std::vector<const Calo1DCluster*> Calo2DCluster::getCluster() const
  { 
	std::vector<const Calo1DCluster*> m_1dclusters; 
	m_1dclusters.clear();
	m_1dclusters.insert(m_1dclusters.end(),m_u1dclusters.begin(),m_u1dclusters.end());
	m_1dclusters.insert(m_1dclusters.end(),m_v1dclusters.begin(),m_v1dclusters.end());
	return  m_1dclusters;
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