#ifndef _CALOCLUSTER_
#define _CALOCLUSTER_
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "CaloBar.h"
#include <stdlib.h>

#include <algorithm>

using namespace std;
class CaloCluster{

public: 
  CaloCluster() {getEcalpar(); vector<CaloBar> getNeighbor(CaloBar& _ibar);};

  void Clear() { Bars.clear(); pos.SetXYZ(0.,0.,0.); phi=0; z=0; chi2=-99; alpha=-99; beta=-99;}
  vector<CaloBar> getBars() { return Bars; }

  void getEcalpar(){

   // char *command;
  //  strcpy(command, "python ../script/readxml.py");
  //  system(command);    

    ifstream file("/afs/ihep.ac.cn/users/z/zhaoxiao/cms/CEPCSW/Digitisers/FanEcalDigi/script/Ecal_Rotated_parameter.txt");
    int i;
    string s;
    string name;
    double value;

    for(i=0;getline(file,s);++i){
        istringstream sin(s);
        sin>>name>>value;
        if(name == "Ecal_crystal_y_width"){
           Ecal_crystal_y_width = value;
        }
        if(name == "Ecal_crystal_rotate_angle"){
           Ecal_alpha = value;
        }
        if(name == "Ecal_crystal_envelope_length"){
           Ecal_crystal_envelope_length = value;
        }
        if(name == "Ecal_barrel_outer_radius_redef"){
           Ecal_rmax = value;
        }
        if(name == "Ecal_barrel_inner_radius"){
           Ecal_rmin = value;
        }
        if(name == "Ecal_barrel_thickness"){
           Ecal_barrel_thickness = value;
        }
        if(name == "Ecal_barrel_outer_radius"){
           Ecal_barrel_outer_radius = value;
        }
        if(name == "Ecal_barrel_half_length"){
           Ecal_barrel_half_length = value;
        }
	if(name == "nphi"){
           Ecal_nphi = int(value);
	}
        if(name == "numberZ"){
           Ecal_nz = value;
        }
        if(name == "zhalf"){
           Ecal_zhalf = value;
        }
      }
   } 

  void FitProfile(){

  //  TCanvas *c1 = new TCanvas();
    double initPhi = getPhiStart();  
    double X0 = 11.2;   //Radiation length, unit: mm
    TH1D *h_hitz = new TH1D("h_hitz", "h_hitz", 30, 0, 30);
    for(int i=0; i<Bars.size(); i++){
//      double relR = (Bars[i].getCrystal()*180./571 - initPhi)*1900*3.1415/(180*X0);
  	double relR = (Bars[i].getCrystal()*ang_pi/int(Ecal_nphi/2) - initPhi)*Ecal_rmin*pi/(ang_pi*X0);
        h_hitz->Fill(relR, Bars[i].getEnergy());
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
  //  c1->SaveAs("/afs/ihep.ac.cn/users/z/zhaoxiao/cms/CEPCSW/Detector/DetCRD/scripts/data/h_hitz.png");

    chi2 = func->GetChisquare()/func->GetNDF();
    alpha = func->GetParameter(0);
    beta = func->GetParameter(1);

    delete h_hitz;
  //  delete c1;
  }

  void Fit2Profile(){

    TCanvas *c1 = new TCanvas();
    double initPhi = getPhiStart();
   
    double X0 = 11.2;   //Radiation length, unit: mm
    TH1D *h_hitz = new TH1D("h_hitz", "h_hitz", 30, 0, 30);

    for(int i=0; i<Bars.size(); i++){
	double relR = (Bars[i].getCrystal()*ang_pi/(Ecal_nphi/2) - initPhi)*Ecal_rmin*pi/(ang_pi*X0);
        h_hitz->Fill(relR, Bars[i].getEnergy());
    }
    h_hitz->Draw();

    TF1 *func = new TF1("fc1", "([2]*[1]*pow([1]*x, ([0]-1))*exp(-[1]*x)/TMath::Gamma([0])) + ([5]*[4]*pow([4]*(x+[6]), ([3]-1))*exp(-[4]*(x+[6]))/TMath::Gamma([3]))", 0, 30);
    func->SetParName(0, "alpha1");
    func->SetParName(1, "beta1");
    func->SetParName(2, "E1");
    func->SetParName(3, "alpha2");
    func->SetParName(4, "beta2");
    func->SetParName(5, "E2");
    func->SetParName(6, "t0");
    func->SetParameter(0, 4.8);
    func->SetParameter(1, 0.5);
    func->SetParameter(2, 10);
    func->SetParameter(3, 4.8);
    func->SetParameter(4, 0.5);
    func->SetParameter(5, 10);
    func->SetParameter(6, 0.1);
    h_hitz->Fit("fc1","Q","",0,30);
    c1->SaveAs("/afs/ihep.ac.cn/users/z/zhaoxiao/cms/CEPCSW/Detector/DetCRD/scripts/data/2h_hitz.png");

    chi2 = func->GetChisquare()/func->GetNDF();
    alpha = func->GetParameter(0);
    beta = func->GetParameter(1);

    delete h_hitz;
    }

  double getPhiStart(){
    vector<int> phiID; phiID.clear();
    double totE = getEnergy();
    for(int i=0; i<Bars.size(); i++)
      if(Bars[i].getEnergy()>totE*0.01) phiID.push_back(Bars[i].getCrystal());
    sort(phiID.begin(), phiID.end());
    if(phiID.size()==0) return -99;
    return (double)phiID[0]*180./int(Ecal_nphi/2) ; //571.int(Ecal_nphi/2.)
  }

  double getPhiStartFit(){
    FitProfile();
    double totE = getEnergy();
    double avePhi = 0;
    for(int i=0; i<Bars.size(); i++) avePhi += Bars[i].getEnergy()*Bars[i].getCrystal()*(180./int(Ecal_nphi/2))/totE;//571. 

    double X0 = 11.2;   //Radiation length, unit: mm
    double depth_Phi = (alpha/beta)*X0/Ecal_rmin * (ang_pi/pi);  // 
    return avePhi - depth_Phi;
  }

  double getZ(){
    double aveZ=0;
    double totE = getEnergy();
    for(int i=0; i<Bars.size(); i++) aveZ += Bars[i].getEnergy()*Bars[i].getPosition().z()/totE; 
    return aveZ; 
  }


  double getAvePhi(){
    double avephi=0;
    double subTotE = getSubEnergy();
    for(int i=0; i<Bars.size(); i++) {
	if (Bars[i].getEnergy() > 0.05) {
	  //  subTotE += Bars[i].getEnergy();
	    avephi += (Bars[i].getCrystal()*(180./int(Ecal_nphi/2))*Bars[i].getEnergy())/subTotE;
	}
//	avephi = avephi/subTotE;
    }
    return avephi;
  }

  double getEnergy(){
    double totE = 0;
    for(int i=0; i<Bars.size(); i++) totE += Bars[i].getEnergy();
    return totE; 
  }
  
  double getSubEnergy(){
    double subTotE = 0;
    for(int i=0; i<Bars.size(); i++) {
	if(Bars[i].getEnergy() > 0.05) subTotE += Bars[i].getEnergy();
    }
    return subTotE;
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

  vector<CaloBar> getNeighbor(CaloBar& _ibar){
	n_Bars.clear();
    for(int i=0; i<Bars.size(); i++){
      if(Bars[i].isNeighbor(_ibar)){
	n_Bars.push_back(Bars[i]);
      }
    }
    return n_Bars;
  }
private: 
  vector<CaloBar> Bars;
  vector<CaloBar> n_Bars;
  TVector3 pos;
  double phi;
  double z; 
  int id;
  
  double chi2;
  double alpha;
  double beta;

  double Ecal_crystal_y_width;
  double Ecal_barrel_inner_radius;//mm from ~/cms/CEPCSW/Detector/DetCRD/compact/Standalone/Dimensions_v01_01.xml
  double Ecal_barrel_outer_radius;
  double Ecal_barrel_thickness;
  double Ecal_barrel_half_length;
  double Ecal_crystal_envelope_length; 
  double Ecal_zhalf;//Ecal_barrel_half_length_correct from ~/cms/CEPCSW/Detector/DetCRD/compact/CRD_common_v01/Ecal_Rotated_Crystal_v01_01.xml
  double Ecal_alpha;//Ecal_crystal_rotate_angle
  int Ecal_nphi;
  double Ecal_rmin;//Ecal_barrel_inner_radius
  double Ecal_rmax;//Ecal_barrel_outer_radius_redef
  int Ecal_nz;//numberZ

  const double pi = 3.1415926;
  double ang_pi = 180.;
};

#endif
