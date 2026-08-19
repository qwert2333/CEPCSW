#ifndef _CALOCLUSTER_
#define _CALOCLUSTER_
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "CaloBar.h"
#include "CaloCrossBar.h"
#include "CaloCrossCluster.h"
#include <stdlib.h>

#include <algorithm>

using namespace std;
class CaloCluster{

public: 
  CaloCluster() {getEcalpar(); vector<CaloBar> getNeighbor(CaloBar& _ibar); vector<CaloCrossBar> getCrossPoint( CaloBar& _ibar ); };

  void Clear() { Bars.clear(); pos.SetXYZ(0.,0.,0.); phi=0; z=0; chi2=-99; alpha=-99; beta=-99;}
  vector<CaloBar> getBars() { return Bars; }
 
  int getClusterID() {return Bars[0].getClusterID(); }

  void getEcalpar(){
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
 
  double getPhiStart(){
    vector<int> phiID; phiID.clear();
    double totE = getEnergy();
    for(int i=0; i<Bars.size(); i++)
      if(Bars[i].getEnergy()>totE*0.002) phiID.push_back(Bars[i].getCrystal());
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

  double getClusModule(){
    double module = 0;
    double totE = getEnergy();
    for(int i=0; i<Bars.size(); i++) module += Bars[i].getEnergy()*((Bars[i].getModule()-336.)*10.)/totE;
    return module;
  }

  double getETHAvePhi(double E_threshold = 0){
    double thr_phi = 0;
    double subTotE = getSubEnergy(E_threshold);
    for(int i=0; i<Bars.size(); i++){
       if(Bars[i].getEnergy() > E_threshold) {
           thr_phi += (Bars[i].getCrystal()*(180./int(Ecal_nphi/2))*Bars[i].getEnergy())/subTotE;
       }
    }
    return thr_phi;
  }

  double getAvePhi(bool useEnergy = true){
    double avephi=0;
    double E = getEnergy();
    for(int i=0; i<Bars.size(); i++) {
       if(useEnergy){
         avephi += (Bars[i].getCrystal()*(180./int(Ecal_nphi/2))*Bars[i].getEnergy())/E;
       }else{
         avephi += Bars[i].getCrystal()*(180./int(Ecal_nphi/2));
      }
    }
    if(!useEnergy){avephi = avephi/Bars.size();}
    return avephi;
  }

  double getMaxPhi(){
    double maxPhi = 0;
    double maxE = -1;
    for(int i=0; i<Bars.size(); i++){
       if(Bars[i].getEnergy() >= maxE){
	  maxE = Bars[i].getEnergy();
	  maxPhi = Bars[i].getCrystal()*(180./int(Ecal_nphi/2));
       }
    } 
    return maxPhi;
  }

  double getEnergy(){
    double totE = 0;
    for(int i=0; i<Bars.size(); i++) totE += Bars[i].getEnergy();
    return totE; 
  }
  
  double getSubEnergy(double E_th = 0){
    double subTotE = 0;
    for(int i=0; i<Bars.size(); i++) {
	if(Bars[i].getEnergy() > E_th) subTotE += Bars[i].getEnergy();
    }
    return subTotE;
  }
//=================================21/21==================================
/* 
  vector<int> getClusModule(){
     vector<int> module_ID; module_ID.clear();
     for(int i=0; i<Bars.size(); i++){
	module_ID.push_back(Bars[i].getModule());
     }
     sort(module_ID.begin(),module_ID.end());
     module_ID.erase(unique(module_ID.begin(),module_ID.end()),module_ID.end());
     return module_ID;
  }
*/
  double ToZ(int module){   
     return (module-336.)*10.;
  }  
 
 vector<CaloCrossBar> getCrossPoint( CaloBar& _ibar ){
    dd4hep::Position c_pos;
    CaloCrossBar c_bar;
    vector<CaloCrossBar> bars; bars.clear();
    CaloCrossCluster clus_bar; clus_bar.Clear();
 
    double O_phi=0;
    double M_phi=0;
    double delta_phi = 0;
    double bar_R=0;
    double real_R = 0;
    double bar_phi=0;
    double bar_z=0;
    double bar_E = 0;
    double D1=0;

    for(int i=0; i<Bars.size(); i++) {  
      if(Bars[i].isNextModule(_ibar)){
         O_phi = _ibar.getCrystal()/1276.*360.;
         M_phi = Bars[i].getCrystal()/1276.*360.;
         delta_phi = abs(O_phi-M_phi);

        if(_ibar.getModule()%2 == 0 && (_ibar.getCrystal()-Bars[i].getCrystal()<20&&_ibar.getCrystal()-Bars[i].getCrystal()>=0)){
           bar_z = (ToZ(_ibar.getModule())+ToZ(Bars[i].getModule()))/2.;
           bar_phi = O_phi-(delta_phi/2.);
           D1 = Ecal_rmin*sin((delta_phi/2.)*toRad)/sin((20.-delta_phi/2.)*toRad);
           bar_R = sqrt(Ecal_rmin*Ecal_rmin + D1*D1 + 2*D1*Ecal_rmin*cos(20.*toRad));
           real_R = sqrt(bar_R*bar_R + bar_z*bar_z);
           bar_E = _ibar.getEnergy()+Bars[i].getEnergy();
           c_pos.SetXYZ(real_R,bar_phi,bar_z);
           c_bar.setPosition(c_pos);
           c_bar.setEnergy(bar_E);
           bars.push_back(c_bar);
           clus_bar.AddBar(c_bar);
           if(_ibar.getModule()>335&&_ibar.getModule()<337){
          cout<<"real R: "<<real_R<<endl;
          cout<<"bar_R: "<<c_bar.getPositionR()<<" bar_phi: "<<c_bar.getPositionPhi()<<" bar_z: "<<bar_z<<" bar_E: "<<bar_E<<endl;
           cout<<"c1:"<<_ibar.getCrystal()<<" m1: "<<_ibar.getModule()<<" E1:"<<_ibar.getEnergy()<<" c2:"<<Bars[i].getCrystal()<<" m2:"<<Bars[i].getModule()<<" E2:"<<Bars[i].getEnergy()<<endl; 
           }
        }
        if(_ibar.getModule()%2 != 0 && (Bars[i].getCrystal()-_ibar.getCrystal()<20&&Bars[i].getCrystal()-_ibar.getCrystal()>=0)){
           D1 = Ecal_rmin*sin((delta_phi/2.)*toRad)/sin((20.-delta_phi/2.)*toRad);
           bar_R = sqrt(Ecal_rmin*Ecal_rmin + D1*D1 + 2*D1*Ecal_rmin*cos(20.*toRad));
           bar_phi = O_phi+(delta_phi/2.);
           bar_z = (ToZ(_ibar.getModule())+ToZ(Bars[i].getModule()))/2.;
           real_R = sqrt(bar_R*bar_R + bar_z*bar_z);
           bar_E = _ibar.getEnergy()+Bars[i].getEnergy();
           c_pos.SetXYZ(real_R,bar_phi,bar_z);
           c_bar.setPosition(c_pos);
           c_bar.setEnergy(bar_E);
           bars.push_back(c_bar);
           clus_bar.AddBar(c_bar);
        }
      }
    }
    return bars;
 }

  double getCrossR(vector<CaloCrossBar>& m_cross_bars){
    CaloCrossCluster m_cross_clus;  m_cross_clus.Clear();
    for(int i=0;i<m_cross_bars.size();i++){
	m_cross_clus.AddBar(m_cross_bars[i]);
    }
    cout<<"????????????????? R: "<<m_cross_clus.getCrossR(true)<<endl;
    return m_cross_clus.getCrossR(true);
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
  dd4hep::Position c_pos;
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

  double pi = 3.14;
  double ang_pi = 180.;
  double toAng = ang_pi/pi;
  double toRad = pi/ang_pi;
  double nlayer = 10;
};

#endif
