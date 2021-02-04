#ifndef _CRD_CALOBARSHOWER_C
#define _CRD_CALOBARSHOWER_C

#include "CRDEcalEDMSvc/CRDCaloBarShower.h"
#include <cmath>
#include <vector>

namespace CRDEcalEDM {

    bool CRDCaloBarShower::isNeighbor(CRDEcalEDM::CRDCaloBar iBar){
      for(int i=0;i<Bars.size();i++){
        if(!inCluster(iBar) && (iBar.getBar()==Bars[i].getBar()+1 || iBar.getBar()==Bars[i].getBar()-1) ) return true;
      }
      return false;
    }

    bool CRDCaloBarShower::inCluster(CRDEcalEDM::CRDCaloBar iBar){
      for(int i=0;i<Bars.size();i++){
        if(iBar.getcellID() == Bars[i].getcellID() ) return true;
      }
      return false;
    }

    double CRDCaloBarShower::getE(){
      double E=0;
      for(int i=0;i<Bars.size();i++) E+=Bars[i].getEnergy();
      return E;
    }

    dd4hep::Position CRDCaloBarShower::getPos(){
      dd4hep::Position pos(0,0,0);
      double Etot=getE();
      for(int i=0;i<Bars.size();i++) pos += (Bars[i].getPosition() * Bars[i].getEnergy())/Etot;
      return pos;
    }

    double CRDCaloBarShower::getT1(){
      double T1=0;
      double Etot = getE();
      for(int i=0;i<Bars.size();i++) T1 += (Bars[i].getT1() * Bars[i].getEnergy())/Etot;
      return T1;
    }

    double CRDCaloBarShower::getT2(){
      double T2=0;
      double Etot = getE();
      for(int i=0;i<Bars.size();i++) T2 += (Bars[i].getT2() * Bars[i].getEnergy())/Etot;
      return T2;
    }


};
#endif
