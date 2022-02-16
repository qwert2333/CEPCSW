#ifndef _CALOCLUSTER_
#define _CALOCLUSTER_
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "CaloBar.h"
#include <algorithm>

using namespace std;
class CaloCluster{

public: 
  CaloCluster() {};

  void Clear() { Bars.clear(); pos.SetXYZ(0.,0.,0.); phi=0; z=0; chi2=-99; alpha=-99; beta=-99;}
  vector<CaloBar> getBars() { return Bars; }

  void FitProfile(){

    //TCanvas *c1 = new TCanvas();
    double initPhi = getPhiStart(); 
    double X0 = 11.2;   //Radiation length, unit: mm
    TH1D *h_hitz = new TH1D("h_hitz", "h_hitz", 30, 0, 30);
    for(int i=0; i<Bars.size(); i++){
      double relR = (Bars[i].getCrystal()*180./571 - initPhi)*1900*3.1415/(180*X0);
      h_hitz->Fill(relR, Bars[i].getEnergy() );
    }    
    h_hitz->Draw();

    TF1 *func = new TF1("fc1", "[2]*[1]*pow([1]*x, ([0]-1))*exp(-[1]*x)/TMath::Gamma([0])", 0, 30);
    func->SetParName(0, "alpha");
    func->SetParName(1, "beta");
    func->SetParName(2, "E0");
    func->SetParameter(0, 4.8);
    func->SetParameter(1, 0.5);
    func->SetParameter(2, 10);
    h_hitz->Fit("fc1","Q","",0,30);
    //c1->SaveAs("/cefs/higgs/guofy/CEPCSW_FanEcal/run/LongiFit.C");

    chi2 = func->GetChisquare()/func->GetNDF();
    alpha = func->GetParameter(0);
    beta = func->GetParameter(1);

    delete h_hitz;
    //delete c1;
  }

  double getPhiStart(){
    vector<int> phiID; phiID.clear();
    double totE = getEnergy();
    for(int i=0; i<Bars.size(); i++)
      if(Bars[i].getEnergy()>totE*0.01) phiID.push_back(Bars[i].getCrystal());
    sort(phiID.begin(), phiID.end());
    if(phiID.size()==0) return -99;
    return (double)phiID[0]*180./(571.) ; 
  }

  double getPhiStartFit(){
    FitProfile();
    double totE = getEnergy();
    double avePhi = 0;
    for(int i=0; i<Bars.size(); i++) avePhi += Bars[i].getEnergy()*Bars[i].getCrystal()*(180/571.)/totE; 

    double X0 = 11.2;   //Radiation length, unit: mm
    double depth_Phi = (alpha/beta)*X0/1900. * (180./3.14159);  // 
    return avePhi - depth_Phi;
  }

  double getZ(){
    double aveZ=0;
    double totE = getEnergy();
    for(int i=0; i<Bars.size(); i++) aveZ += Bars[i].getEnergy()*Bars[i].getPosition().z()/totE; 
    return aveZ; 
  }

  double getEnergy(){
    double totE = 0;
    for(int i=0; i<Bars.size(); i++) totE += Bars[i].getEnergy();
    return totE; 
  }

  double getFitAlpha() {return alpha;}
  double getFitBeta() {return beta;}
  double getFitChi2() {return chi2;}

  bool inCluster( CaloBar& _ibar ){
    for(int i=0; i<Bars.size(); i++) 
      if(_ibar==Bars[i]) return true; 
    return false; 
  }

  bool isNeighbor( CaloBar& _ibar ){
    if(inCluster(_ibar)) return false; 

    for(int i=0; i<Bars.size(); i++)
      if( Bars[i].isNeighbor(_ibar) ) return true; 
    
    return false; 
  }

  void AddBar( CaloBar& _ibar) { Bars.push_back(_ibar); }


private: 
  vector<CaloBar> Bars;
  TVector3 pos;
  double phi;
  double z; 

  double chi2;
  double alpha;
  double beta;

};

#endif
