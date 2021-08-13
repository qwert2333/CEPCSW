#ifndef _CRD_CALOBARSHOWER_C
#define _CRD_CALOBARSHOWER_C

#include "Objects/CRDCaloBarShower.h"
#include <cmath>
#include <vector>

namespace CRDEcalEDM {

    bool CRDCaloBarShower::isNeighbor(CRDEcalEDM::CRDCaloBar iBar){
      for(int i=0;i<Bars.size();i++){
        if(!inShower(iBar) && (iBar.getBar()==Bars[i].getBar()+1 || iBar.getBar()==Bars[i].getBar()-1) ) return true;
      }
      return false;
    }

    bool CRDCaloBarShower::inShower(CRDEcalEDM::CRDCaloBar iBar){
      for(int i=0;i<Bars.size();i++){
        if(iBar.getcellID() == Bars[i].getcellID() ) return true;
      }
      return false;
    }

    double CRDCaloBarShower::getE() const{
      double E=0;
      for(int i=0;i<Bars.size();i++) E+=Bars[i].getEnergy();
      return E;
    }

    dd4hep::Position CRDCaloBarShower::getPos() const{
      dd4hep::Position pos(0,0,0);
      double Etot=getE();
      for(int i=0;i<Bars.size();i++) pos += (Bars[i].getPosition() * Bars[i].getEnergy())/Etot;
      return pos;
    }
   
    double CRDCaloBarShower::getT1() const{
      double T1=0;
      double Etot = getE();
      for(int i=0;i<Bars.size();i++) T1 += (Bars[i].getT1() * Bars[i].getEnergy())/Etot;
      return T1;
    }

    double CRDCaloBarShower::getT2() const{
      double T2=0;
      double Etot = getE();
      for(int i=0;i<Bars.size();i++) T2 += (Bars[i].getT2() * Bars[i].getEnergy())/Etot;
      return T2;
    }

    double CRDCaloBarShower::getWidth() const{
      dd4hep::Position centPos = getPos(); 
      double Etot = getE(); 
      double sigmax=0; 
      double sigmay=0; 
      double sigmaz=0; 
      for(int i=0; i<Bars.size(); i++){
        double wi = Bars[i].getEnergy()/Etot;
        sigmax += wi*(Bars[i].getPosition().x()-centPos.x())*(Bars[i].getPosition().x()-centPos.x());
        sigmay += wi*(Bars[i].getPosition().y()-centPos.y())*(Bars[i].getPosition().y()-centPos.y());
        sigmaz += wi*(Bars[i].getPosition().z()-centPos.z())*(Bars[i].getPosition().z()-centPos.z());
      }

      if(sigmaz!=0) return sigmaz;  //sLayer=1, bars along z-axis. 
      else if(sigmax==0 && sigmaz==0) return sigmay; //Module 2, 6
      else if(sigmay==0 && sigmaz==0) return sigmax; //Module 0, 4
      else if(sigmax!=0 && sigmay!=0 && sigmaz==0) return sqrt(sigmax*sigmax+sigmay*sigmay); //Module 1, 3, 5, 7; 
      else return 0.; 

    }


    int CRDCaloBarShower::getDlayer() const{
      if(Bars.size()>0) return Bars[0].getDlayer();
      else return -1;
    }

    std::vector<CRDEcalEDM::CRDShadowCluster> CRDCaloBarShower::getAllCandiCol() const{
      std::vector<CRDEcalEDM::CRDShadowCluster> m_allcandi; m_allcandi.clear();
      m_allcandi = NeuCandidateCol;
      m_allcandi.insert(m_allcandi.end(), TrkCandidateCol.begin(), TrkCandidateCol.end());
      return m_allcandi;
    }



};
#endif
