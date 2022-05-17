#ifndef _CRD_CALOHIT_3DSHOWER_H
#define _CRD_CALOHIT_3DSHOWER_H

#include "k4FWCore/DataHandle.h"
//#include "edm4hep/CalorimeterHit.h"
//#include "edm4hep/CalorimeterHitCollection.h"
//#include "edm4hep/SimCalorimeterHitCollection.h"
#include "Objects/CRDCaloHitTransShower.h"
#include "Objects/TrackFitInEcal.h"

#include "TMath.h"
#include "TVector3.h"
#include "TString.h"

using namespace std;
namespace CRDEcalEDM{

  class CRDCaloHit3DCluster{
  public:
    CRDCaloHit3DCluster(){};

    //edm4hep::ConstCalorimeterHit getClusterInitialHit() const; 
    //TVector3 getHitCenter() const;
    TVector3 getShowerCenter() const;
    //double getHitsE() const;
    double getShowerE() const;

    double getExpEnergy(int& dlayer) const;
    TVector3 getExpPos(int& dlayer) const;

    //void FitProfile();
    //void FitAxis(); 
    void AddShower(CRDEcalEDM::CRDCaloHitTransShower& shower);
    void MergeCluster(CRDEcalEDM::CRDCaloHit3DCluster& clus);
    void IdentifyCluster(); 

    inline bool operator == ( const CRDEcalEDM::CRDCaloHit3DCluster &x) const {
      //TVector3 vec = this->getHitCenter();
      //double En = this->getHitsE();
      //return ( vec==x.getHitCenter() && En==x.getHitsE());
      return true;
    }

    std::vector<CRDEcalEDM::CRDCaloHitTransShower> get2DShowers() const { return ShowerinLayer; }
    //std::vector<edm4hep::ConstCalorimeterHit> getCaloHits() const { return CaloHits; }
    double getFitAlpha() const { return alpha; }
    double getFitBeta() const { return beta; }
    double getFitExpEn() const {return ExpEn; }
    double getShowerMax() const {return showerMax; }
    double getChi2() const { return chi2; }
    TVector3 getAxis() const { return axis; }
    void getFitPars(double& _alpha, double& _beta ) const { _alpha = alpha; _beta = beta; }
    int getBeginningDlayer() const; 
    int getEndDlayer() const;

    std::vector<double> getEnInLayer() const; 
    int getMaxELayer() const; 
    double getAveE() const; 
    double getStdDevE() const; 

    //std::vector<double> getClusterWidth() const;
    double getMaxWidth() const; 
    int getType() const { return type; }


    void setShowerVec( std::vector<CRDEcalEDM::CRDCaloHitTransShower>& _svec ) { ShowerinLayer = _svec;}
    //void setCaloHits( std::vector<edm4hep::ConstCalorimeterHit>& _hits ) { CaloHits = _hits; }
    void setType( int _t ) { type = _t; }
    void Clear() { ShowerinLayer.clear(); /*CaloHits.clear()*/; showerEnergy=0; hitsEnergy=0; showerMax=0; chi2=0; alpha=0; beta=0; ExpEn=0; axis.SetXYZ(0,0,0); type=-1; }

  private: 
    std::vector<CRDEcalEDM::CRDCaloHitTransShower> ShowerinLayer;
    //std::vector<edm4hep::ConstCalorimeterHit> CaloHits;
    TVector3 axis;
    double showerEnergy;
    double hitsEnergy;
    double showerMax;
    double chi2;
    double alpha;
    double beta;
    int    type;  //0: MIP shower.  1: EM shower.  2: hadronic shower
    double ExpEn; 

    TrackFitInEcal* track = new TrackFitInEcal(); 

  };
};
#endif
