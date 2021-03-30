#ifndef _CRD_CALOBAR_CLUSTER_C
#define _CRD_CALOBAR_CLUSTER_C

#include "CRDEcalEDMSvc/CRDCaloBar.h"
#include "CRDEcalEDMSvc/CRDCaloBarCluster.h"

namespace CRDEcalEDM{

  void CRDCaloBarCluster::Clear(){
    Bars.clear();
    Seeds.clear();
    Energy=0.;
    pos.SetXYZ(0.,0.,0.);
    Nseeds=0;
    ScndMoment=0;
  }

  bool CRDCaloBarCluster::isNeighbor(CRDEcalEDM::CRDCaloBar iBar){
    for(int i=0;i<Bars.size();i++){
      if(!inCluster(iBar) && (iBar.getBar()==Bars[i].getBar()+1 || iBar.getBar()==Bars[i].getBar()-1) ) return true;
    }
    return false;
  }

  bool CRDCaloBarCluster::inCluster(CRDEcalEDM::CRDCaloBar iBar){
    for(int i=0;i<Bars.size();i++){
      if(iBar.getcellID() == Bars[i].getcellID()) return true;
    }
    return false;
  }

  double CRDCaloBarCluster::getE() const{
    double E=0;
    for(int i=0;i<Bars.size();i++) E+=Bars[i].getEnergy();
    return E;
  }

  dd4hep::Position CRDCaloBarCluster::getPos() const{
    dd4hep::Position pos(0,0,0);
    double Etot=getE();
    for(int i=0;i<Bars.size();i++) pos += (Bars[i].getPosition() * Bars[i].getEnergy())/Etot;
    return pos;
  }

  double CRDCaloBarCluster::getScndMoment(){
    dd4hep::Position pos = getPos();
    double Etot = getE();
    double scndM = 0;
    for(int i=0;i<Bars.size();i++) scndM += (Bars[i].getEnergy() * (pos-Bars[i].getPosition()).Mag2()) / Etot;
    ScndMoment = scndM;
    return ScndMoment;
  }

  bool CRDCaloBarCluster::getGlobalRange( double& xmin,  double& ymin, double& zmin, double& xmax, double& ymax, double& zmax ) const{
    if(Bars.size()==0) return false; 

    std::vector<double> posx; posx.clear();
    std::vector<double> posy; posx.clear();
    std::vector<double> posz; posx.clear();
    for(int i=0; i<Bars.size(); i++){  
      posx.push_back(Bars[i].getPosition().x());  
      posy.push_back(Bars[i].getPosition().y());
      posz.push_back(Bars[i].getPosition().z());
    }
    auto xminptr = std::min_element( std::begin(posx), std::end(posx) );
    auto xmaxptr = std::max_element( std::begin(posx), std::end(posx) );
    auto yminptr = std::min_element( std::begin(posy), std::end(posy) );
    auto ymaxptr = std::max_element( std::begin(posy), std::end(posy) );
    auto zminptr = std::min_element( std::begin(posz), std::end(posz) );
    auto zmaxptr = std::max_element( std::begin(posz), std::end(posz) );

    int slayer = Bars[0].getSlayer(); 
    if(slayer==0){
      xmin = -9999.;    xmax = 9999.; 
      ymin = -9999.;    ymax = 9999.; 
      zmin = *zminptr;  zmax = *zmaxptr;
    }
    else if(slayer==1){
      xmin = *xminptr-6;  xmax = *xmaxptr+6;
      ymin = *yminptr-6;  ymax = *ymaxptr+6;
      zmin = -9999;     zmax = 9999; 
    }
    else return false; 

    return true;
  }


  void CRDCaloBarCluster::PrintBars(){
    if(Bars.size()!=0){
      std::cout<<"BarCluster::PrintBars"<<std::endl;
      printf("#bar \t x \t y \t z \t E \n");
      for(int i=0;i<Bars.size();i++) printf("%d \t %f \t %f \t %f \t %f \n",i, Bars[i].getPosition().x(), Bars[i].getPosition().y(), Bars[i].getPosition().z(), Bars[i].getEnergy() );
    }
  }

  void CRDCaloBarCluster::PrintSeeds(){
    if(Seeds.size()!=0){
      std::cout<<"BarCluster::PrintSeeds"<<std::endl;
      printf("#Seed \t x \t y \t z \t E \n");
      for(int i=0;i<Seeds.size();i++) printf("%d \t %f \t %f \t %f \t %f \n",i, Seeds[i].getPosition().x(), Seeds[i].getPosition().y(), Seeds[i].getPosition().z(), Seeds[i].getEnergy() );
    }
  }


};
#endif
