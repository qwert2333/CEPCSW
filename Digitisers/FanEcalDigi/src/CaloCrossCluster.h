#ifndef _CALOCROSSCLUSTER_
#define _CALOCROSSCLUSTER_
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "CaloBar.h"
#include "CaloCrossBar.h"
#include <stdlib.h>

#include <algorithm>

using namespace std;
class CaloCrossCluster{

public: 
  CaloCrossCluster() {getEcalpar(); };

  void Clear() { cBars.clear(); pos.SetXYZ(0.,0.,0.); phi=0; z=0;}
  vector<CaloCrossBar> getcBars() { return cBars; }

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
 
  double getCrossRStart(){
    vector<double> RID; RID.clear();
    double totE = getEnergy();
    for(int i=0; i<cBars.size(); i++) RID.push_back(cBars[i].getPositionR());
    sort(RID.begin(), RID.end());
    if(RID.size()==0) return -99;
    return (double)RID[0];
  }

  double getCrossR(bool useEnergy = true){
    double aveR=0;
    double num = 0;
//    double subTotE = getSubEnergy();
    double totE = getEnergy();
    for(int i=0; i<cBars.size(); i++){
       if(useEnergy){
       aveR += cBars[i].getPositionR()*cBars[i].getEnergy()/totE;
       }else{
        aveR += cBars[i].getPositionR();
       }
        num+=1;
    }
    if(!useEnergy){aveR = aveR/num;}
    return aveR;
  }
  
  double getCrossZ(){
    double aveZ=0;
    double totE = getEnergy();
    for(int i=0; i<cBars.size(); i++) aveZ += cBars[i].getEnergy()*cBars[i].getPositionZ()/totE; 
    return aveZ; 
  }

  double getAvePhi(bool useEnergy = true){
    double avephi=0;
    double num = 0;
    double subTotE = getSubEnergy();
    for(int i=0; i<cBars.size(); i++) {
	if (cBars[i].getEnergy() > 0.001) {
	    if(useEnergy){
	    avephi += cBars[i].getPositionPhi()*cBars[i].getEnergy()/subTotE;
	    }else{
	     avephi += cBars[i].getPositionPhi();
 	    }
	num+=1;
	}	
    }
    if(!useEnergy){avephi = avephi/num;}
    return avephi;
  }

  double getSubEnergy(){
    double subtotE = 0;
    for(int i=0; i<cBars.size(); i++){
	if(cBars[i].getEnergy() > 0.001) subtotE += cBars[i].getEnergy();
    }
    return subtotE;
  }

  double getEnergy(){
    double totE = 0;
    for(int i=0; i<cBars.size(); i++) totE += cBars[i].getEnergy();
    return totE; 
  }
  
 CaloCrossBar getMaxEBar(){
   
    double max_bar_E = 0; 
    int idx = 0;
    CaloCrossBar bar;
    for(int i=0; i<cBars.size(); i++) {
       if(cBars[i].getEnergy()>max_bar_E) {
  //       cout<<"E: "<<cBars[i].getEnergy()<<" ; R:"<<cBars[i].getPositionR()<<" Phi: "<<cBars[i].getPositionPhi()<<" Z:"<<cBars[i].getPositionZ()<<endl;
         max_bar_E = cBars[i].getEnergy();
         idx = i;
       }
    }
    bar = cBars[idx];
    return bar;
  }

  void SortBarR(vector<CaloCrossBar>& bars){
   for (int i=0;i<bars.size();i++){
     int k = i;
     for (int j = k+1;j<bars.size();j++){
        if (bars[j].getPositionR()<bars[k].getPositionR()){
           k = j;
        }
     }
     CaloCrossBar tempBar = bars[i];
     bars[i] = bars[k];
     bars[k] = tempBar;
   }
 }

  vector<CaloCrossBar> removeDuplicates(vector<CaloCrossBar>& bars){
    if(bars.size()==0){
       cout<<"vector size = 0"<<endl;
    }
    int idx = 0;
    for(int i=1; i<bars.size(); i++){
       if(bars[idx].getPositionR() != bars[i].getPositionR()){
          idx++;
          bars[idx] = bars[i];
       }
    } 
    return bars; 
  }

  double getShowerStartR(){
    CaloCrossBar max_bar = getMaxEBar();
    double max_bar_R = max_bar.getPositionR();
    double max_bar_phi = max_bar.getPositionPhi();
    double max_bar_z = max_bar.getPositionZ();
    double max_bar_E = max_bar.getEnergy();
  //  cout<<"maxE:"<<max_bar_E<<" ; phi:"<<max_bar_phi<<" Z: "<<max_bar_z<<endl;
    vector<CaloCrossBar> ShowerBars; ShowerBars.clear();
    vector<CaloCrossBar> lshowerBars; lshowerBars.clear();
    vector<CaloCrossBar> rshowerBars; rshowerBars.clear();

    for(int i=0; i<cBars.size(); i++) {
       if(max_bar_z<0){
    //     cout<<"xiaoyu 0: "<<cBars[i].getPositionZ()<<" ; "<<cBars[i].getPositionPhi()<<endl;
           if((cBars[i].getPositionZ()<(max_bar_z+1)&&cBars[i].getPositionZ()>(max_bar_z-1))/* && cBars[i].getPositionPhi()==max_bar_phi*/) {
               ShowerBars.push_back(cBars[i]);
    //      cout<<"XE: "<<cBars[i].getEnergy()<<" ; R:"<<cBars[i].getPositionR()<<" Phi: "<<cBars[i].getPositionPhi()<<" Z:"<<cBars[i].getPositionZ()<<endl;
       }
           if((cBars[i].getPositionZ()<(max_bar_z+11)&&cBars[i].getPositionZ()>(max_bar_z+9))/* &&cBars[i].getPositionPhi()==max_bar_phi*/) {
             lshowerBars.push_back(cBars[i]);
      //       cout<<"YE: "<<cBars[i].getEnergy()<<" ; R:"<<cBars[i].getPositionR()<<" Phi: "<<cBars[i].getPositionPhi()<<" Z:"<<cBars[i].getPositionZ()<<endl;
           }
           if((cBars[i].getPositionZ()<(max_bar_z-9)&&cBars[i].getPositionZ()>(max_bar_z-11))/* &&cBars[i].getPositionPhi()==max_bar_phi*/) {
             rshowerBars.push_back(cBars[i]);
           } 
       } 
       if(max_bar_z>0){
      //     cout<<"dayu 0: "<<cBars[i].getPositionZ()<<" ; "<<cBars[i].getPositionPhi()<<endl;
             if((cBars[i].getPositionZ()<(max_bar_z+1)&&cBars[i].getPositionZ()>(max_bar_z-1))/* && cBars[i].getPositionPhi()==max_bar_phi*/) {
          ShowerBars.push_back(cBars[i]);
      //    cout<<"XE: "<<cBars[i].getEnergy()<<" ; R:"<<cBars[i].getPositionR()<<" Phi: "<<cBars[i].getPositionPhi()<<" Z:"<<cBars[i].getPositionZ()<<endl;
             } 
           if((cBars[i].getPositionZ()<(max_bar_z-9)&&cBars[i].getPositionZ()>(max_bar_z-11))/* &&cBars[i].getPositionPhi()==max_bar_phi*/) {
             lshowerBars.push_back(cBars[i]);
        //     cout<<"YE: "<<cBars[i].getEnergy()<<" ; R:"<<cBars[i].getPositionR()<<" Phi: "<<cBars[i].getPositionPhi()<<" Z:"<<cBars[i].getPositionZ()<<endl;
            }
            if((cBars[i].getPositionZ()<(max_bar_z+11)&&cBars[i].getPositionZ()>(max_bar_z+9))/* &&cBars[i].getPositionPhi()==max_bar_phi*/) {
               rshowerBars.push_back(cBars[i]);
            }
       }
    }
    SortBarR(ShowerBars);
    SortBarR(lshowerBars);
    SortBarR(rshowerBars);
    removeDuplicates(ShowerBars);  
    removeDuplicates(lshowerBars);
    removeDuplicates(rshowerBars);
  /*  for(int m=0;m<ShowerBars.size();m++){
      cout<<"aaaa: "<<ShowerBars[m].getPositionR()<<"; z: "<<ShowerBars[m].getPositionZ()<<"; phi: "<<ShowerBars[m].getPositionPhi()<<"; E: "<<ShowerBars[m].getEnergy()<<endl;
    }
*/
    SortBarR(ShowerBars);
    SortBarR(lshowerBars);
    SortBarR(rshowerBars);

    double minR = ShowerBars[0].getPositionR();
    double minE = ShowerBars[0].getEnergy();
    double minR1 = lshowerBars[0].getPositionR();
    double minE1 = lshowerBars[0].getEnergy();
    double minR2 = rshowerBars[0].getPositionR();
    double minE2 = rshowerBars[0].getEnergy();
    for(int j=0;j<ShowerBars.size();j++){
     if(ShowerBars[j].getEnergy()<0.1){
       if(ShowerBars[j].getEnergy()>0.001 && ShowerBars[j+1].getEnergy()>3*ShowerBars[j].getEnergy() && ShowerBars[j+2].getEnergy()>3*ShowerBars[j].getEnergy()){
        cout<<"1: "<<ShowerBars[j].getPositionR()<<"; z: "<<ShowerBars[j].getPositionZ()<<"; phi: "<<ShowerBars[j].getPositionPhi()<<"; E: "<<ShowerBars[j].getEnergy()<<endl;
        minR = ShowerBars[j].getPositionR();
        minE = ShowerBars[j].getEnergy();
        break;
        }
      } 
      if(ShowerBars[j].getEnergy()>0.1){
        if(ShowerBars[j+1].getEnergy()>ShowerBars[j].getEnergy()&&ShowerBars[j+2].getEnergy()>ShowerBars[j].getEnergy()){
           minR = ShowerBars[j].getPositionR();
           minE = ShowerBars[j].getEnergy();
           cout<<">0.1     "<<ShowerBars[j].getPositionR()<<"; z: "<<ShowerBars[j].getPositionZ()<<"; phi: "<<ShowerBars[j].getPositionPhi()<<"; E: "<<ShowerBars[j].getEnergy()<<endl; 
           break;  
        }
      }
   }
   for(int k=0;k<lshowerBars.size();k++){  
      if(lshowerBars[k].getEnergy()<0.1){
       if(lshowerBars[k].getEnergy()>0.001 && lshowerBars[k+1].getEnergy()>3*lshowerBars[k].getEnergy() && lshowerBars[k+2].getEnergy()>3*lshowerBars[k].getEnergy()){
        cout<<"2: "<<lshowerBars[k].getPositionR()<<"; z: "<<lshowerBars[k].getPositionZ()<<"; phi: "<<lshowerBars[k].getPositionPhi()<<"; E: "<<lshowerBars[k].getEnergy()<<endl;
        minR1 = lshowerBars[k].getPositionR();
        minE1 = lshowerBars[k].getEnergy();
        break;
        }
      }
      if(lshowerBars[k].getEnergy()>0.1){
        if(lshowerBars[k+1].getEnergy()>lshowerBars[k].getEnergy()&&lshowerBars[k+2].getEnergy()>lshowerBars[k].getEnergy()){
           minR1 = lshowerBars[k].getPositionR();
           minE1 = lshowerBars[k].getEnergy();
           cout<<"2: >0.1    "<<lshowerBars[k].getPositionR()<<"; z: "<<lshowerBars[k].getPositionZ()<<"; phi: "<<lshowerBars[k].getPositionPhi()<<"; E: "<<lshowerBars[k].getEnergy()<<endl;
           break;
        }
      }
    }
     for(int k=0;k<rshowerBars.size();k++){
      if(rshowerBars[k].getEnergy()<0.1){
       if(rshowerBars[k].getEnergy()>0.001 && rshowerBars[k+1].getEnergy()>3*rshowerBars[k].getEnergy() && rshowerBars[k+2].getEnergy()>3*rshowerBars[k].getEnergy()){
        cout<<"2: "<<rshowerBars[k].getPositionR()<<"; z: "<<rshowerBars[k].getPositionZ()<<"; phi: "<<rshowerBars[k].getPositionPhi()<<"; E: "<<rshowerBars[k].getEnergy()<<endl;
        minR2 = rshowerBars[k].getPositionR();
        minE2 = rshowerBars[k].getEnergy();
        break;
        }
      }
      if(rshowerBars[k].getEnergy()>0.1){
        if(rshowerBars[k+1].getEnergy()>rshowerBars[k].getEnergy()&&rshowerBars[k+2].getEnergy()>rshowerBars[k].getEnergy()){
           minR2 = rshowerBars[k].getPositionR();
           minE2 = rshowerBars[k].getEnergy();
           cout<<"2: >0.1    "<<rshowerBars[k].getPositionR()<<"; z: "<<rshowerBars[k].getPositionZ()<<"; phi: "<<rshowerBars[k].getPositionPhi()<<"; E: "<<rshowerBars[k].getEnergy()<<endl;
           break;
        }
      }
    }

   double tempE = 0; 
    if(minE1>minE && minE1>minE2){
       minR = minR1;
    }
    if(minE2>minE && minE2>minE1){
       minR = minR2;
    }

    cout<<"========minR========"<<minR<<endl;
    return minR;
  }

/*  double getScaleEnergy(double f_E){
    double totE = 0;
    totE = f_E*getEnergy();
    return totE;
  }
*/
  void AddBar( CaloCrossBar& _ibar) { cBars.push_back(_ibar); }

private: 
  vector<CaloCrossBar> cBars;
  TVector3 pos;
  dd4hep::Position c_pos;
  double phi;
  double z; 
  int id;

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
  double toAng = ang_pi/pi;
  double toRad = pi/ang_pi;
  double nlayer = 10;
};

#endif
