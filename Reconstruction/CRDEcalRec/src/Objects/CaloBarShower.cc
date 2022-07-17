#ifndef CALOBARSHOWER_C
#define CALOBARSHOWER_C

#include "Objects/CaloBarShower.h"
#include <cmath>
#include <vector>

namespace PandoraPlus {

  void CaloBarShower::Clear(){
    Bars.clear();
  }

 void CaloBarShower::Clean(){
    delete Seed; Seed = NULL;
    for(int i=0; i<Bars.size(); i++) { delete Bars[i]; Bars[i] = NULL; }
    Clear();
  }

  void CaloBarShower::Check(){
    for(int i=0; i<Bars.size(); i++)
      if(!Bars[i]) { Bars.erase(Bars.begin()+i); i--; }
  }

/*    bool CaloBarShower::isNeighbor(PandoraPlus::CaloUnit iBar){
      for(int i=0;i<Bars.size();i++){
        if(!inShower(iBar) && (iBar.getBar()==Bars[i].getBar()+1 || iBar.getBar()==Bars[i].getBar()-1) ) return true;
      }
      return false;
    }

    bool CaloBarShower::inShower(PandoraPlus::CaloUnit iBar){
      for(int i=0;i<Bars.size();i++){
        if(iBar.getcellID() == Bars[i].getcellID() ) return true;
      }
      return false;
    }
*/
    double CaloBarShower::getE() const{
      double E=0;
      for(int i=0;i<Bars.size();i++) E+=Bars[i]->getEnergy();
      return E;
  }

  TVector3 CaloBarShower::getPos() const{
    TVector3 pos(0,0,0);
    double Etot=getE();
    for(int i=0;i<Bars.size();i++) pos += Bars[i]->getPosition() * (Bars[i]->getEnergy()/Etot);
    return pos;
  }

  double CaloBarShower::getT1() const{
    double T1=0;
    double Etot = getE();
    for(int i=0;i<Bars.size();i++) T1 += (Bars[i]->getT1() * Bars[i]->getEnergy())/Etot;
    return T1;
  }

  double CaloBarShower::getT2() const{
    double T2=0;
    double Etot = getE();
    for(int i=0;i<Bars.size();i++) T2 += (Bars[i]->getT2() * Bars[i]->getEnergy())/Etot;
    return T2;
  }

  double CaloBarShower::getWidth() const{
    TVector3 centPos = getPos(); 
    double Etot = getE(); 
    double sigmax=0; 
    double sigmay=0; 
    double sigmaz=0; 
    for(int i=0; i<Bars.size(); i++){
      double wi = Bars[i]->getEnergy()/Etot;
      sigmax += wi*(Bars[i]->getPosition().x()-centPos.x())*(Bars[i]->getPosition().x()-centPos.x());
      sigmay += wi*(Bars[i]->getPosition().y()-centPos.y())*(Bars[i]->getPosition().y()-centPos.y());
      sigmaz += wi*(Bars[i]->getPosition().z()-centPos.z())*(Bars[i]->getPosition().z()-centPos.z());
    }

    if(sigmaz!=0) return sigmaz;  //sLayer=1, bars along z-axis. 
    else if(sigmax==0 && sigmaz==0) return sigmay; //Module 2, 6
    else if(sigmay==0 && sigmaz==0) return sigmax; //Module 0, 4
    else if(sigmax!=0 && sigmay!=0 && sigmaz==0) return sqrt(sigmax*sigmax+sigmay*sigmay); //Module 1, 3, 5, 7; 
    else return 0.; 

  }

  void CaloBarShower::setIDInfo(){
    if(Bars.size()==0) { std::cout<<"WARNING: Empty shower does not have ID Info! "<<std::endl;  return; }
    module = Bars[0]->getModule(); 
    stave = Bars[0]->getStave();
    part = Bars[0]->getPart(); 
    dlayer = Bars[0]->getDlayer(); 
    slayer = Bars[0]->getSlayer(); 
  }

};
#endif
