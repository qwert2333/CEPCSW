#ifndef _CALOBAR_
#define _CALOBAR_

#include <DD4hep/Objects.h>
#include "TVector3.h"
#include "TRandom.h"

using namespace dd4hep;
using namespace std;
using namespace dd4hep::detail;

class CaloBar{

public:
  CaloBar(unsigned long long _cellID, int _system, int _module, int _crystal, dd4hep::Position _pos, double _Q1, double _Q2, double _T1, double _T2, int _clus_id)
  : cellID(_cellID), system(_system), module(_module), crystal(_crystal), position(_pos), Q1(_Q1), Q2(_Q2), T1(_T1), T2(_T2) ,clus_id(_clus_id) {}; 
  CaloBar() {};

  inline bool operator == (const CaloBar &x) const{
    return cellID == x.cellID;
  }
  unsigned long long getcellID() const { return cellID; }
  int getSystem()  const { return system; }
  int getModule()  const { return module; }
  int getCrystal() const { return crystal;}
  double getQ1()   const { return Q1;     }
  double getQ2()   const { return Q2;     }
  double getT1()   const { return T1;     }
  double getT2()   const { return T2;     }
  int getClusterID() const { return clus_id; }

  dd4hep::Position getPosition() const { return position; }
  double getDigitizedEnergy() const { return Digi_frac*(Q1+Q2)/2.; } //E after digitized
  double getEnergy() const { return (Q1+Q2)/2.; }

//=======================================SiPM readout====================================
 double SiPMDigi(double edepCry){
	int sPix = 0;
	sPix = gRandom->Poisson(edepCry/fMIPEnergy*fCryMeanLY);
	if (sPix<5)
	    return 0;
	double sAdcMean = sPix * fADCMean;
	double sAdcSigma = sqrt(sPix) * fADCSigma;
	double sAdc = -1;
	while(sAdc<0)
	   sAdc = gRandom->Gaus(sAdcMean, sAdcSigma);
	double sMIP = sAdc/fADCMean/fCryMeanLY;
	if (sMIP < MIP_Thre)
	   return 0;
	return sMIP * fMIPEnergy;
 }
 
  double getSiPMEnergy() { return SiPMDigi(1000.*(Q1+Q2)/2.)/(Cab_frac*1000.); }

  bool isNeighbor(CaloBar &x){
    if( (fabs(x.getModule()-module)==1 || fabs(x.getCrystal()-crystal)==1) ) return true;
    else return false; 
  }

  bool isNextModule(CaloBar &x){
    if( fabs(x.getModule()-module)==1 ) return true;
   //  if( x.getModule()-module==1 ) return true;
     else return false;
  }
  
  bool isCrossNeighbor(CaloBar &x){
    if( abs(x.getModule()-module)==1 && abs(x.getCrystal()-crystal)<=10) return true;
    else return false;
  }

  void setcellID(unsigned long long _cellid) { cellID = _cellid; }
  void setcellID(int _system, int _module, int _crystal ) { system=_system; module=_module; crystal=_crystal; }
  void setPosition( dd4hep::Position pos) { position = pos; }
  void setPosition( TVector3 posv3) { position.SetXYZ( posv3.x(), posv3.y(), posv3.z() ); }
  void setQ(double _q1, double _q2) { Q1=_q1; Q2=_q2; }
  void setT(double _t1, double _t2) { T1=_t1; T2=_t2; }
  void setClusterID(int _clus_id) { clus_id=_clus_id; }

private:
	unsigned long long cellID;
	int system;
	int module;
	int crystal;
	double E;
	dd4hep::Position position;
	double Q1;      // Q in left readout
	double Q2;      // Q in right readout;
	double T1;    // T in left readout;
	double T2;    // T in right readout;
	int clus_id;     //cluster id of bar;
 
        double fCryMeanLY = 850; // 850 p.e./MIP, 3x3mm SiPM
        double fMIPEnergy = 9; // 9 MeV / MIP in 1 cm BGO
        double fADCMean = 15;
        double fADCSigma = 3;
        double MIP_Thre = 0.3; // 0.3 MIP     
        double Cab_frac = 0.952; // CALIBRATION
        double Digi_frac = 1.041;
};
  
#endif
