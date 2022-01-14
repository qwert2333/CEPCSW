#ifndef _CRD_CALOHIT_LONGICLUS_C
#define _CRD_CALOHIT_LONGICLUS_C

#include "Objects/CRDCaloHitLongiCluster.h"
#include <cmath>

namespace CRDEcalEDM{

  int CRDCaloHitLongiCluster::getSlayer() const{
    if(barShowerCol.size()>0) return barShowerCol[0].getSlayer();
    else return -1;
  }

  TVector3 CRDCaloHitLongiCluster::getPos() const{
    TVector3 pos(0, 0, 0);
    for(int i=0; i<barShowerCol.size(); i++){
      TVector3 m_pos(barShowerCol[i].getPos().x(), barShowerCol[i].getPos().y(), barShowerCol[i].getPos().z());
      pos += m_pos;
    }
    return pos; 
  }

  int CRDCaloHitLongiCluster::getBeginningDlayer() const{
    int Lstart = 99; 
    for(int i=0; i<barShowerCol.size(); i++)
      if(barShowerCol[i].getDlayer()<Lstart) Lstart = barShowerCol[i].getDlayer(); 
    if(Lstart==99) return -99;
    else return Lstart;  
  }


  int CRDCaloHitLongiCluster::getEndDlayer() const{
    int Lend = -99;
    for(int i=0; i<barShowerCol.size(); i++)
      if(barShowerCol[i].getDlayer()>Lend) Lend = barShowerCol[i].getDlayer(); 
    return Lend;
  }



  TVector3 CRDCaloHitLongiCluster::getExpPos(int& dlayer) const{
    double maxE = -99;
    int dlayer_maxE = 0;
    TVector3 pos_maxE(0,0,0);

    if(dlayer<getBeginningDlayer() || dlayer>getEndDlayer()) return pos_maxE;

    for(int i=0; i<barShowerCol.size(); i++){
      if(barShowerCol[i].getE()>maxE){
        maxE = barShowerCol[i].getE(); 
        dlayer_maxE = barShowerCol[i].getDlayer();
        pos_maxE.SetXYZ( barShowerCol[i].getPos().x(), barShowerCol[i].getPos().y(), barShowerCol[i].getPos().z() );
    }}

    return pos_maxE + axis*20*(dlayer-dlayer_maxE);
  }

  bool CRDCaloHitLongiCluster::isContinue() const{
    int ly_max = -1; 
    int ly_min = 99;
    std::vector<int> vec_layers; vec_layers.clear();
    for(int il=0; il<barShowerCol.size(); il++){
      vec_layers.push_back(barShowerCol[il].getDlayer());
      if(barShowerCol[il].getDlayer()>ly_max) ly_max = barShowerCol[il].getDlayer();
      if(barShowerCol[il].getDlayer()<ly_min) ly_min = barShowerCol[il].getDlayer();
    }

    bool flag = true;
    std::sort(vec_layers.begin(), vec_layers.end());
    for(int il=0; il<vec_layers.size() && vec_layers[il]!=ly_max; il++)
      if( find(vec_layers.begin(), vec_layers.end(), vec_layers[il]+1) == vec_layers.end() ){ flag=false; break; }

    return flag;
  }

  bool CRDCaloHitLongiCluster::isContinueN(int n) const{
    if(n<=0) return true; 

    int ly_max = -1;
    int ly_min = 99;
    std::vector<int> vec_layers; vec_layers.clear();
    for(int il=0; il<barShowerCol.size(); il++){
      vec_layers.push_back(barShowerCol[il].getDlayer());
      if(barShowerCol[il].getDlayer()>ly_max) ly_max = barShowerCol[il].getDlayer();
      if(barShowerCol[il].getDlayer()<ly_min) ly_min = barShowerCol[il].getDlayer();
    }
    if(n>(ly_max-ly_min+1)) return false; 

    bool flag = false;
    std::sort(vec_layers.begin(), vec_layers.end());
    for(int il=0; il<vec_layers.size(); il++){
      bool fl_continueN = true;
      for(int in=1; in<n; in++) 
        if( find(vec_layers.begin(), vec_layers.end(), vec_layers[il]+in)!=vec_layers.end() ) { fl_continueN=false; break; }
      if(fl_continueN) {flag = true; break;}
    }
    return flag;
  }

  bool CRDCaloHitLongiCluster::isSubset(CRDCaloHitLongiCluster& clus) const{
    bool flag = true; 
    for(int is=0; is<clus.getBarShowers().size(); is++)
      if( find(barShowerCol.begin(), barShowerCol.end(), clus.getBarShowers()[is])==barShowerCol.end() ) {flag=false; break;}
    return flag; 
  }

  void CRDCaloHitLongiCluster::FitAxis(){
    if(barShowerCol.size()==0){ axis.SetXYZ(0,0,0); return; }
    else if(barShowerCol.size()==1){ 
      axis.SetXYZ(barShowerCol[0].getPos().x(), barShowerCol[0].getPos().y(), barShowerCol[0].getPos().z());
      axis = axis.Unit(); 
      return;
    }
    else if(barShowerCol.size()==2){
      dd4hep::Position rpos = barShowerCol.back().getPos() - barShowerCol.front().getPos();
      axis.SetXYZ( rpos.x(), rpos.y(), rpos.z() );
      axis = axis.Unit();
      return; 
    }
    else{
      track->clear(); 
      double barAngle = (barShowerCol[0].getModule()+2)*TMath::Pi()/4.;
      double posErr = 10./sqrt(12);
      if(barAngle>=TMath::TwoPi()) barAngle = barAngle-TMath::TwoPi();
      track->setBarAngle(barAngle);
      for(int is=0; is<barShowerCol.size(); is++){
        dd4hep::Position b_pos = barShowerCol[is].getPos(); 
        track->setGlobalPoint(0, b_pos.x(), posErr, b_pos.y(), posErr, b_pos.z(), posErr);
        track->setGlobalPoint(1, b_pos.x(), posErr, b_pos.y(), posErr, b_pos.z(), posErr);
      }

      track->fitTrack();
      double fitPhi = track->getTrkPar(2);
      double fitTheta = track->getTrkPar(3);

      axis.SetMag(1.); 
      axis.SetPhi(fitPhi);
      axis.SetTheta(fitTheta);
    }
  }





};
#endif
