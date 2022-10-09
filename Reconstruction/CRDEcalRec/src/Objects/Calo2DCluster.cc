#ifndef CALO_2DCLUSTER_C
#define CALO_2DCLUSTER_C

#include "Objects/Calo2DCluster.h"
#include <cmath>
using namespace std;
namespace PandoraPlus{

  void Calo2DCluster::Clear() {
    towerID.clear();
    //m_modules.clear(); 
    //m_parts.clear(); 
    //m_staves.clear();
    barUCol.clear();
    barVCol.clear();
    barShowerUCol.clear(); 
    barShowerVCol.clear();  
    barClusterUCol.clear(); 
    barClusterVCol.clear(); 
  }

  void Calo2DCluster::ClearShower() {
    barShowerUCol.clear(); 
    barShowerVCol.clear();
  }

  void Calo2DCluster::Check(){
    for(int i=0; i<barClusterVCol.size(); i++)
      if(!barClusterVCol[i]) { barClusterVCol.erase(barClusterVCol.begin()+i); i--; }
	for(int i=0; i<barClusterUCol.size(); i++)
      if(!barClusterUCol[i]) { barClusterUCol.erase(barClusterUCol.begin()+i); i--; }
    for(int i=0; i<barUCol.size(); i++)
      if(!barUCol[i]) { barUCol.erase(barUCol.begin()+i); i--; }	  
    for(int i=0; i<barVCol.size(); i++)
      if(!barVCol[i]) { barVCol.erase(barVCol.begin()+i); i--; }
    for(int i=0; i<barShowerUCol.size(); i++)
      if(!barShowerUCol[i]) { barShowerUCol.erase(barShowerUCol.begin()+i); i--; }
    for(int i=0; i<barShowerVCol.size(); i++)
      if(!barShowerVCol[i]) { barShowerVCol.erase(barShowerVCol.begin()+i); i--; }
  	}

  void Calo2DCluster::Clean(){
	for(int i=0; i<barClusterVCol.size(); i++){ delete barClusterVCol[i]; barClusterVCol[i]=NULL; }
	for(int i=0; i<barClusterUCol.size(); i++){ delete barClusterUCol[i]; barClusterUCol[i]=NULL; }
    for(int i=0; i<barUCol.size(); i++) { delete barUCol[i]; barUCol[i]=NULL; }
    for(int i=0; i<barVCol.size(); i++) { delete barVCol[i]; barVCol[i]=NULL; }
    for(int i=0; i<barShowerUCol.size(); i++) { delete barShowerUCol[i]; barShowerUCol[i]=NULL; }
    for(int i=0; i<barShowerVCol.size(); i++) { delete barShowerVCol[i]; barShowerVCol[i]=NULL; }
    //std::vector<int>().swap(m_modules);
    //std::vector<int>().swap(m_parts);
    //std::vector<int>().swap(m_staves);
    Clear();
  	}

  bool Calo2DCluster::isNeighbor(const PandoraPlus::Calo1DCluster* m_1dcluster) const{
    assert(m_1dcluster->getBars().size() > 0 && getCluster().at(0)->getBars().size()>0 );
    if(m_1dcluster->getDlayer() != getDlayer()  ) return false; 

    for(int i=0; i<m_1dcluster->getTowerID().size(); i++){
		for(int j=0; j<towerID.size(); j++){
      if( m_1dcluster->getTowerID()[i]==towerID[j] ) return true;
    }}
    
    return false;
  }

 
  void Calo2DCluster::addUnit(const Calo1DCluster* _1dcluster)
  {
    if(_1dcluster->getSlayer()==0) barClusterUCol.push_back(_1dcluster);
    if(_1dcluster->getSlayer()==1) barClusterVCol.push_back(_1dcluster); 
    for(int ib=0; ib<_1dcluster->getBars().size(); ib++) addBar(_1dcluster->getBars()[ib]);

    std::vector< std::vector<int> > id = _1dcluster->getTowerID();
    for(int ii=0; ii<id.size(); ii++) 
      if( find(towerID.begin(), towerID.end(), id[ii])==towerID.end() ) towerID.push_back(id[ii]);

    //std::vector<int> m_1dmodules = _1dcluster->getModules();
    //std::vector<int> m_1dparts = _1dcluster->getParts();
    //std::vector<int> m_1dstaves = _1dcluster->getStaves();
    //m_modules.insert(m_modules.end(),m_1dmodules.begin(),m_1dmodules.end());
    //m_parts.insert(m_parts.end(),m_1dparts.begin(),m_1dparts.end());
    //m_staves.insert(m_staves.end(),m_1dstaves.begin(),m_1dstaves.end());
  }

  std::vector<const Calo1DCluster*> Calo2DCluster::getCluster() const{ 
    std::vector<const Calo1DCluster*> m_1dclusters; 
    m_1dclusters.clear();
    m_1dclusters.insert(m_1dclusters.end(),barClusterUCol.begin(),barClusterUCol.end());
    m_1dclusters.insert(m_1dclusters.end(),barClusterVCol.begin(),barClusterVCol.end());
    return  m_1dclusters;
  }
	
  std::vector<const CaloUnit*> Calo2DCluster::getBars() const
  {
	std::vector<const CaloUnit*> results;
	results.clear();
	std::vector<const Calo1DCluster*> m_1dclusters; 
	m_1dclusters.clear();
	m_1dclusters.insert(m_1dclusters.end(),barClusterUCol.begin(),barClusterUCol.end());
	m_1dclusters.insert(m_1dclusters.end(),barClusterVCol.begin(),barClusterVCol.end());
	for(int i=0; i<m_1dclusters.size(); i++)
	{
		for(int j=0; j<m_1dclusters.at(i)->getBars().size(); j++)
		{
			results.push_back(m_1dclusters.at(i)->getBars().at(j));
		}
	}
	return results;
  }
	
  double Calo2DCluster::getEnergy() const {
    double sumE = 0;

    //Priority order: barShower -> barCluster. 
    if(barShowerUCol.size()==0 && barShowerVCol.size()==0){
      for(int m=0; m<barClusterUCol.size(); m++) sumE += barClusterUCol[m]->getEnergy();
      for(int m=0; m<barClusterVCol.size(); m++) sumE += barClusterVCol[m]->getEnergy();
    }
    else{
      for(int m=0; m<barShowerUCol.size(); m++) sumE += barShowerUCol[m]->getEnergy();
      for(int m=0; m<barShowerVCol.size(); m++) sumE += barShowerVCol[m]->getEnergy();
    }
    return sumE;
  }


  TVector3 Calo2DCluster::getPos() const{
    TVector3 m_pos(0, 0, 0); 
    if(towerID.size()==0) return m_pos; 

    float rotAngle = towerID[0][0]*TMath::Pi()/4.;
    TVector3 m_vecX(0., 0., 0.);  
    TVector3 m_vecY(0., 0., 0.);
    if(barShowerUCol.size()==0 && barShowerVCol.size()==0){
      for(int m=0; m<barClusterUCol.size(); m++) m_vecX += barClusterUCol[m]->getPos();
      m_vecX *= (1./barClusterUCol.size());
      for(int m=0; m<barClusterVCol.size(); m++) m_vecY += barClusterVCol[m]->getPos();
      m_vecY *= (1./barClusterVCol.size());
    }
    else{
      for(int m=0; m<barShowerUCol.size(); m++) m_vecX += barShowerUCol[m]->getPos();
      m_vecX *= (1./barShowerUCol.size());
      for(int m=0; m<barShowerVCol.size(); m++) m_vecY += barShowerVCol[m]->getPos();
      m_vecY *= (1./barShowerVCol.size());
    }
    m_vecX.RotateZ(rotAngle);
    m_vecY.RotateZ(rotAngle);
    m_pos.SetXYZ( m_vecY.x(), (m_vecX.y()+m_vecY.y())/2 , m_vecX.z() );
    m_pos.RotateZ(-rotAngle);
    return m_pos;
  }


};
#endif
