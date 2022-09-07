#ifndef CALO_3DCLUSTER_C
#define CALO_3DCLUSTER_C

#include "Objects/Calo3DCluster.h"
#include <cmath>
using namespace std;
namespace PandoraPlus{

  void Calo3DCluster::Clear() 
  {
	  m_2dclusters.clear();
	  m_towers.clear();
    towerID.clear(); 
	  //m_modules.clear();
	  //m_parts.clear();
	  //m_staves.clear();
  }

  void Calo3DCluster::Clean(){
    for(int i=0; i<m_2dclusters.size(); i++) { delete m_2dclusters[i]; m_2dclusters[i]=NULL; }
    //for(int i=0; i<m_towers.size(); i++) { delete m_towers[i]; m_towers[i]=NULL; }
    //std::vector<int>().swap(m_modules);
    //std::vector<int>().swap(m_parts);
    //std::vector<int>().swap(m_staves);
    Clear();
  }

  void Calo3DCluster::Check()
  {
    for(int i=0; i<m_2dclusters.size(); i++)
    if(!m_2dclusters[i]) { m_2dclusters.erase(m_2dclusters.begin()+i); i--; }
    //for(int i=0; i<m_towers.size(); i++)
    //if(!m_towers[i]) { m_towers.erase(m_towers.begin()+i); i--; }
  }

  bool Calo3DCluster::isNeighbor(const PandoraPlus::Calo2DCluster* m_2dcluster) const
  {
    //Inner module
    for(int i=0; i<m_2dcluster->getTowerID().size(); i++){
    for(int j=0; j<towerID.size(); j++){
      if( m_2dcluster->getTowerID()[i]==towerID[j] ) return true;
    }}

    //for(int i=0; i<m_2dcluster->getModules().size(); i++){
    //  for(int j=0; j<m_modules.size(); j++){
    //    if(m_2dcluster->getModules().at(i)==m_modules.at(j) && m_2dcluster->getParts().at(i)==m_parts.at(j) && m_2dcluster->getStaves().at(i)==m_staves.at(j)){
    //      return true;
    //  }}
    //}

    //Inner module but in adjacent Dlayer
    for(int i=0; i<m_2dclusters.size(); i++){
      if(fabs(m_2dcluster->getDlayer()-m_2dclusters[i]->getDlayer())>1) continue; 

      std::vector<const PandoraPlus::CaloUnit*> m_EdgeBarU_clus1; m_EdgeBarU_clus1.clear(); //Income 2DCluster
      std::vector<const PandoraPlus::CaloUnit*> m_EdgeBarV_clus1; m_EdgeBarV_clus1.clear();
      std::vector<const PandoraPlus::CaloUnit*> m_EdgeBarU_clus2; m_EdgeBarU_clus2.clear(); //2DCluster in present 3DCluster. 
      std::vector<const PandoraPlus::CaloUnit*> m_EdgeBarV_clus2; m_EdgeBarV_clus2.clear();
      
      for(int ib=0; ib<m_2dcluster->getBarUCol().size(); ib++){
        if(m_2dcluster->getBarUCol()[ib]->isAtLowerEdgeZ() || m_2dcluster->getBarUCol()[ib]->isAtUpperEdgeZ()) m_EdgeBarU_clus1.push_back(m_2dcluster->getBarUCol()[ib]);
      }
      for(int ib=0; ib<m_2dcluster->getBarVCol().size(); ib++){
        if(m_2dcluster->getBarVCol()[ib]->isAtLowerEdgePhi() || m_2dcluster->getBarVCol()[ib]->isAtUpperEdgePhi()) m_EdgeBarV_clus1.push_back(m_2dcluster->getBarVCol()[ib]);
      }
      for(int ib=0; ib<m_2dclusters[i]->getBarUCol().size(); ib++){
        if(m_2dclusters[i]->getBarUCol()[ib]->isAtLowerEdgeZ() || m_2dclusters[i]->getBarUCol()[ib]->isAtUpperEdgeZ()) m_EdgeBarU_clus2.push_back(m_2dclusters[i]->getBarUCol()[ib]);
      }
      for(int ib=0; ib<m_2dclusters[i]->getBarVCol().size(); ib++){
        if(m_2dclusters[i]->getBarVCol()[ib]->isAtLowerEdgePhi() || m_2dclusters[i]->getBarVCol()[ib]->isAtUpperEdgePhi()) m_EdgeBarV_clus2.push_back(m_2dclusters[i]->getBarVCol()[ib]);
      }


      for(int ib1=0; ib1<m_EdgeBarU_clus1.size(); ib1++){
        for(int ib2=0; ib2<m_2dclusters[i]->getBarVCol().size(); ib2++){
          if( m_EdgeBarU_clus1[ib1]->getPart()==m_2dclusters[i]->getBarVCol()[ib2]->getPart() && abs(m_EdgeBarU_clus1[ib1]->getStave()-m_2dclusters[i]->getBarVCol()[ib2]->getStave())==1 ) 
            return true; 
      }}
      for(int ib1=0; ib1<m_EdgeBarV_clus1.size(); ib1++){
        for(int ib2=0; ib2<m_2dclusters[i]->getBarUCol().size(); ib2++){
          if( m_EdgeBarV_clus1[ib1]->getStave()==m_2dclusters[i]->getBarUCol()[ib2]->getStave() && abs(m_EdgeBarV_clus1[ib1]->getPart()-m_2dclusters[i]->getBarUCol()[ib2]->getPart())==1 )
            return true;
      }}
      for(int ib1=0; ib1<m_EdgeBarU_clus2.size(); ib1++){
        for(int ib2=0; ib2<m_2dcluster->getBarVCol().size(); ib2++){
          if( m_EdgeBarU_clus2[ib1]->getPart()==m_2dcluster->getBarVCol()[ib2]->getPart() && abs(m_EdgeBarU_clus2[ib1]->getStave()-m_2dcluster->getBarVCol()[ib2]->getStave())==1 )
            return true;
      }}
      for(int ib1=0; ib1<m_EdgeBarV_clus2.size(); ib1++){
        for(int ib2=0; ib2<m_2dcluster->getBarUCol().size(); ib2++){
          if( m_EdgeBarV_clus2[ib1]->getStave()==m_2dcluster->getBarUCol()[ib2]->getStave() && abs(m_EdgeBarV_clus2[ib1]->getPart()-m_2dcluster->getBarUCol()[ib2]->getPart())==1 )
            return true;
      }}

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

  void Calo3DCluster::addUnit(const Calo2DCluster* _2dcluster){

    m_2dclusters.push_back(_2dcluster);
    std::vector< std::vector<int> > id = _2dcluster->getTowerID();
    for(int ii=0; ii<id.size(); ii++)
      if( find(towerID.begin(), towerID.end(), id[ii])==towerID.end() ) towerID.push_back(id[ii]);

    //std::vector<int> m_2dmodules = _2dcluster->getModules();
    //std::vector<int> m_2dparts = _2dcluster->getParts();
    //std::vector<int> m_2dstaves = _2dcluster->getStaves();
    //m_modules.insert(m_modules.end(),m_2dmodules.begin(),m_2dmodules.end());
    //m_parts.insert(m_parts.end(),m_2dparts.begin(),m_2dparts.end());
    //m_staves.insert(m_staves.end(),m_2dstaves.begin(),m_2dstaves.end());
  }


  std::vector<const PandoraPlus::CaloUnit*> Calo3DCluster::getBars() const{
    std::vector<const PandoraPlus::CaloUnit*> results; results.clear();
    for(int i=0; i<m_2dclusters.size(); i++){
      for(int j=0; j<m_2dclusters.at(i)->getBars().size(); j++){
        results.push_back(m_2dclusters.at(i)->getBars().at(j));
      }
    }
    return results;
  }

  double Calo3DCluster::getEnergy() const{
    double result = 0;
    for(int m=0; m<m_2dclusters.size(); m++)
      result += m_2dclusters[m]->getEnergy();
	
    return result;
  }

  std::vector<const PandoraPlus::LongiCluster*> Calo3DCluster::getLongiClusterUCol(std::string name) const {
    std::vector<const LongiCluster*> emptyCol; emptyCol.clear(); 
    if(map_longiClusUCol.find(name)!=map_longiClusUCol.end()) emptyCol = map_longiClusUCol.at(name);
    return emptyCol;
  }

  std::vector<const PandoraPlus::LongiCluster*> Calo3DCluster::getLongiClusterVCol(std::string name) const {
    std::vector<const LongiCluster*> emptyCol; emptyCol.clear(); 
    if(map_longiClusVCol.find(name)!=map_longiClusVCol.end()) emptyCol = map_longiClusVCol.at(name);
    return emptyCol;
  }

  std::vector<const PandoraPlus::Calo1DCluster*> Calo3DCluster::getLocalMaxUCol(std::string name) const{
    std::vector<const Calo1DCluster*> emptyCol; emptyCol.clear(); 
    if(map_localMaxU.find(name)!=map_localMaxU.end()) emptyCol = map_localMaxU.at(name);
    return emptyCol;
  }

  std::vector<const PandoraPlus::Calo1DCluster*> Calo3DCluster::getLocalMaxVCol(std::string name) const{
    std::vector<const Calo1DCluster*> emptyCol; emptyCol.clear(); 
    if(map_localMaxV.find(name)!=map_localMaxV.end()) emptyCol = map_localMaxV.at(name);
    return emptyCol; 
  }


};
#endif
