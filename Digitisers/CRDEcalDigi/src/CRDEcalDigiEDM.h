#ifndef _CRDECAL_DIGI_EDM_
#define _CRDECAL_DIGI_EDM_

#include "k4FWCore/DataHandle.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include "DetInterface/IGeomSvc.h"
#include "TVector3.h"
#include "TString.h"

using namespace std;
class CRDEcalDigiEDM {
public: 

  class DigiBar {
  public:
		unsigned long long cellID = 0;
		int system;
		int module;
		int dlayer;
		int part;
		int block;
		int slayer;
		int bar;
		dd4hep::Position position;
		double Q1=0;      // Q in left readout
		double Q2=0;      // Q in right readout;
		double T1=999;    // T in left readout;
		double T2=999;    // T in right readout;
		inline bool operator < (const DigiBar &x) const {
			return bar<x.bar ;
		}
		inline bool operator == (const DigiBar &x) const{
			return cellID == x.cellID;
		}


		double getEnergy(){ return (Q1+Q2)/2.; }
      bool isNeighbor(DigiBar &x){
			if(x.cellID!=cellID && 
				x.system==system && 
				x.module==module &&
				x.dlayer==dlayer && 
				x.part==part && 
				x.block==block &&
				x.slayer==slayer &&
				(x.bar==bar+1 || x.bar==bar-1) 
			) return true;
         return false;
      }

  };


  class StepDigiOut {
  public: 
    double Q;
    double T;
    inline bool operator < (const StepDigiOut &x) const { 
      return T<x.T ;
    } 
  };


  class BarCluster{
	public: 
		bool isNeighbor(CRDEcalDigiEDM::DigiBar iBar){
			for(int i=0;i<Bars.size();i++){
				if(!inCluster(iBar) && (iBar.bar==Bars[i].bar+1 || iBar.bar==Bars[i].bar-1) ) return true;
			}
			return false;
		}
		bool inCluster(CRDEcalDigiEDM::DigiBar iBar){
			for(int i=0;i<Bars.size();i++){
				if(iBar.cellID == Bars[i].cellID) return true;
			}
			return false;
		}
		void sortByPos(){
			std::sort(Bars.begin(), Bars.end());
		}
		//void sortByE(){
		//	std::sort(Bars.begin(), Bars.end(), compE);
		//}
		double getE(){
			double E=0;
			for(int i=0;i<Bars.size();i++) E+=Bars[i].getEnergy();
			return E;
		}
		dd4hep::Position getPos(){
			dd4hep::Position pos(0,0,0);
			double Etot=getE();
			for(int i=0;i<Bars.size();i++) pos += (Bars[i].position * Bars[i].getEnergy())/Etot;
			return pos;
		}
		double getScndMoment(){
			dd4hep::Position pos = getPos();
			double Etot = getE();
			double scndM = 0;
			for(int i=0;i<Bars.size();i++) scndM += (Bars[i].getEnergy() * (pos-Bars[i].position).Mag2()) / Etot;
			return scndM;
		}

		//members in this class
		std::vector<CRDEcalDigiEDM::DigiBar> Bars;
		std::vector<CRDEcalDigiEDM::DigiBar> Seeds;
		double Energy;
		dd4hep::Position pos;
		int Nseeds=0;
		double ScndMoment=0.;  //Second moment of this cluster. 
  };


  class BarCollection{
	public: 
      bool isNeighbor(CRDEcalDigiEDM::DigiBar iBar){
         for(int i=0;i<Bars.size();i++){
            if(!inCluster(iBar) && (iBar.bar==Bars[i].bar+1 || iBar.bar==Bars[i].bar-1) ) return true;
         }
         return false;
      }
      bool inCluster(CRDEcalDigiEDM::DigiBar iBar){
         for(int i=0;i<Bars.size();i++){
            if(iBar.cellID == Bars[i].cellID) return true;
         }
         return false;
      }
      double getE(){
         double E=0;
         for(int i=0;i<Bars.size();i++) E+=Bars[i].getEnergy();
         return E;
      }
		dd4hep::Position getPos(){
			dd4hep::Position pos(0,0,0);
			double Etot=getE();
			for(int i=0;i<Bars.size();i++) pos += (Bars[i].position * Bars[i].getEnergy())/Etot;
			return pos;
		}

		double getT1(){
			double T1;
			double Etot = getE();
			for(int i=0;i<Bars.size();i++) T1 += (Bars[i].T1 * Bars[i].getEnergy())/Etot;
			return T1;
		}
		double getT2(){
			double T2;
			double Etot = getE();
			for(int i=0;i<Bars.size();i++) T2 += (Bars[i].T2 * Bars[i].getEnergy())/Etot;
			return T2;
		}

		double Energy;
		dd4hep::Position pos;
		std::vector<CRDEcalDigiEDM::DigiBar> Bars;
		CRDEcalDigiEDM::DigiBar Seed;
  };

	//bool compPos(DigiBar &a, DigiBar &b){ return (a.bar < b.bar); }
	//bool compE(DigiBar &a, DigiBar &b){ return (a.getEnergy() < b.getEnergy()); }


};
#endif
