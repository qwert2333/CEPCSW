#ifndef CALOBAR_CLUSTER_C
#define CALOBAR_CLUSTER_C

#include "Objects/CaloUnit.h"
#include "Objects/CaloBarCluster.h"

namespace PandoraPlus{

  void CaloBarCluster::Clear(){
    Bars.clear();
    Seeds.clear();
    Energy=0.;
    pos.SetXYZ(0.,0.,0.);
  }

  void CaloBarCluster::Clean() {
    for(int i=0; i<Bars.size(); i++) { delete Bars[i]; Bars[i]=NULL; }
    for(int i=0; i<Seeds.size(); i++) { delete Seeds[i]; Seeds[i]=NULL; }
    Clear();
  }

  void CaloBarCluster::Check() {
    for(int i=0; i<Bars.size(); i++)
      if(!Bars[i]) { Bars.erase(Bars.begin()+i); i--; }
    for(int i=0; i<Seeds.size(); i++)
      if(!Seeds[i]) { Seeds.erase(Seeds.begin()+i); i--; }
  }

  bool CaloBarCluster::isNeighbor(const PandoraPlus::CaloUnit* iBar) const{
    for(int i=0;i<Bars.size();i++){
      if(!inCluster(iBar) && (iBar->getBar()==Bars[i]->getBar()+1 || iBar->getBar()==Bars[i]->getBar()-1) ) return true;
    }
    return false;
  }

  bool CaloBarCluster::inCluster(const PandoraPlus::CaloUnit* iBar) const{
    return (find(Bars.begin(), Bars.end(), iBar)!=Bars.end() );
  }

  double CaloBarCluster::getE() const{
    double E=0;
    for(int i=0;i<Bars.size();i++) E+=Bars[i]->getEnergy();
    return E;
  }

  TVector3 CaloBarCluster::getPos() const{
    TVector3 pos(0,0,0);
    double Etot=getE();
    for(int i=0;i<Bars.size();i++) pos += Bars[i]->getPosition() * (Bars[i]->getEnergy()/Etot);
    return pos;
  }

  double CaloBarCluster::getScndMoment() const{
    if(Bars.size()<=1) return 0.;
    TVector3 pos = getPos();
    double Etot = getE();
    double scndM = 0;
    for(int i=0;i<Bars.size();i++) scndM += (Bars[i]->getEnergy() * (pos-Bars[i]->getPosition()).Mag2()) / Etot;
    return scndM;
  }

  bool CaloBarCluster::getGlobalRange( double& xmin,  double& ymin, double& zmin, double& xmax, double& ymax, double& zmax ) const{
    if(Bars.size()==0) return false; 

    std::vector<double> posx; posx.clear();
    std::vector<double> posy; posx.clear();
    std::vector<double> posz; posx.clear();
    for(int i=0; i<Bars.size(); i++){  
      posx.push_back(Bars[i]->getPosition().x());  
      posy.push_back(Bars[i]->getPosition().y());
      posz.push_back(Bars[i]->getPosition().z());
    }
    auto xminptr = std::min_element( std::begin(posx), std::end(posx) );
    auto xmaxptr = std::max_element( std::begin(posx), std::end(posx) );
    auto yminptr = std::min_element( std::begin(posy), std::end(posy) );
    auto ymaxptr = std::max_element( std::begin(posy), std::end(posy) );
    auto zminptr = std::min_element( std::begin(posz), std::end(posz) );
    auto zmaxptr = std::max_element( std::begin(posz), std::end(posz) );

    int slayer = Bars[0]->getSlayer(); 
    if(slayer==0){
      xmin = -9999.;    xmax = 9999.; 
      ymin = -9999.;    ymax = 9999.; 
      zmin = *zminptr-6.;  zmax = *zmaxptr+6.;
    }
    else if(slayer==1){
      xmin = *xminptr-6;  xmax = *xmaxptr+6;
      ymin = *yminptr-6;  ymax = *ymaxptr+6;
      zmin = -9999;     zmax = 9999; 
    }
    else return false; 

    return true;
  }

  int CaloBarCluster::getLeftEdge(){
    std::sort(Bars.begin(), Bars.end());
    if(Bars.size()==0) return -99;
    return Bars[0]->getBar();
  }

  int CaloBarCluster::getRightEdge(){
    std::sort(Bars.begin(), Bars.end());
    if(Bars.size()==0) return -99;
    return Bars[Bars.size()-1]->getBar();
  }  


  void CaloBarCluster::PrintBars() const{
    if(Bars.size()!=0){
      std::cout<<"BarCluster::PrintBars"<<std::endl;
      printf("#bar \t x \t y \t z \t E \t cellID \n");
      for(int i=0;i<Bars.size();i++) printf("%d \t %f \t %f \t %f \t %f \t (%d, %d, %d, %d, %d) \n",i, Bars[i]->getPosition().x(), Bars[i]->getPosition().y(), Bars[i]->getPosition().z(), Bars[i]->getEnergy(), Bars[i]->getModule(), Bars[i]->getStave(), Bars[i]->getDlayer(), Bars[i]->getPart(),  Bars[i]->getSlayer() );
    }
  }

  void CaloBarCluster::PrintSeeds() const{
    if(Seeds.size()>0){
      std::cout<<"BarCluster::PrintSeeds"<<std::endl;
      printf("#Seed \t x \t y \t z \t E \n");
      for(int i=0;i<Seeds.size();i++) printf("%d \t %f \t %f \t %f \t %f \n",i, Seeds[i]->getPosition().x(), Seeds[i]->getPosition().y(), Seeds[i]->getPosition().z(), Seeds[i]->getEnergy() );
    }
  }


};
#endif
