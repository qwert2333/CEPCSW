#ifndef CALO_HALFCLUSTER_C
#define CALO_HALFCLUSTER_C

#include "Objects/CaloHalfCluster.h"

namespace PandoraPlus{

  void CaloHalfCluster::Clear() {
    m_slayer=99;
    m_1dclusters.clear(); 
  }

  void CaloHalfCluster::Check(){
    for(int i=0; i<m_1dclusters.size(); i++)
      if(!m_1dclusters[i]) { m_1dclusters.erase(m_1dclusters.begin()+i); i--; }
  	}

  void CaloHalfCluster::Clean(){
	for(int i=0; i<m_1dclusters.size(); i++){ delete m_1dclusters[i]; m_1dclusters[i]=NULL; }
    Clear();
  	}

  bool CaloHalfCluster::isNeighbor(const PandoraPlus::Calo1DCluster* m_1dcluster) const{
    assert(m_1dcluster->getBars().size() > 0 && getCluster().at(0)->getBars().size()>0 );
    if(m_1dcluster->getSlayer() != getSlayer()  ) return false; 

    for(int i1d=0; i1d<m_1dclusters.size(); i1d++)
    {
      for(int ibar=0; ibar<m_1dcluster->getBars().size(); ibar++)
      {
        for(int jbar=0; jbar<m_1dclusters.at(i1d)->getBars().size(); jbar++)
        {
          if( m_1dcluster->getBars().at(ibar)->isLongiNeighbor(m_1dclusters.at(i1d)->getBars().at(jbar)) ) return true;
          if( m_1dcluster->getBars().at(ibar)->isLongiModuleAdjacent(m_1dclusters.at(i1d)->getBars().at(jbar)) ) return true;
        }
      }
    }

    // std::vector<const PandoraPlus::CaloUnit*> bars_1d = m_1dcluster->getBars();
    // for(int ib1d=0; ib1d<bars_1d.size(); ib1d++)
    // {
    //   for(int ic=0; ic<m_1dclusters.size(); ic++)
    //   {
    //     for(int ib2d=0; ib2d<m_1dclusters[ic]->getBars().size(); ib2d++)
    //     {
    //       if(bars_1d[ib1d]->isModuleAdjacent(m_1dclusters[ic]->getBars()[ib2d])) return true;
    //     }
    //   } 
    // }    
    
    return false;
  }

 
  void CaloHalfCluster::addUnit(const Calo1DCluster* _1dcluster)
  {
    if(_1dcluster->getSlayer()==0) m_slayer=0;
    if(_1dcluster->getSlayer()==1) m_slayer=1;
    m_1dclusters.push_back(_1dcluster);
  }
  
  std::vector<const CaloUnit*> CaloHalfCluster::getBars() const
  {
    std::vector<const CaloUnit*> results;
    results.clear();
    for(int i=0; i<m_1dclusters.size(); i++)
    {
      for(int j=0; j<m_1dclusters.at(i)->getBars().size(); j++)
      {
        results.push_back(m_1dclusters.at(i)->getBars().at(j));
      }
    }
    return results;
  }

  double CaloHalfCluster::getEnergy() const {
    double sumE = 0;
    for(int i=0; i<m_1dclusters.size(); i++)
    {
      sumE = sumE + m_1dclusters.at(i)->getEnergy();
    }
    return sumE;
  }

  std::vector<const PandoraPlus::Calo1DCluster*> CaloHalfCluster::getLocalMaxCol(std::string name) const{
    std::vector<const PandoraPlus::Calo1DCluster*> emptyCol; emptyCol.clear(); 
    if(map_localMax.find(name)!=map_localMax.end()) emptyCol = map_localMax.at(name);
    return emptyCol;
  }

  std::vector<const PandoraPlus::LongiCluster*> CaloHalfCluster::getLongiClusterCol(std::string name) const{
    std::vector<const PandoraPlus::LongiCluster*> emptyCol; emptyCol.clear(); 
    if(map_longiClusCol.find(name)!=map_longiClusCol.end()) emptyCol = map_longiClusCol.at(name);
    return emptyCol;
  }

};
#endif
