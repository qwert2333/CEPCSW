#ifndef _CRD_CALOHIT_3DSHOWER_H
#define _CRD_CALOHIT_3DSHOWER_H

#include "k4FWCore/DataHandle.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "CRDEcalEDMSvc/CRDCaloHit2DShower.h"
#include "CRDEcalEDMSvc/TrackFitInEcal.h"

#include "TMath.h"
#include "TVector3.h"
#include "TString.h"

using namespace std;
namespace CRDEcalEDM{

  class CRDCaloHit3DShower{
  public:
    CRDCaloHit3DShower(){};

    edm4hep::ConstCalorimeterHit getClusterInitialHit() const; 
    TVector3 getHitCenter() const;
    TVector3 getShowerCenter() const;
    double getHitsE() const;
    double getShowerE() const;
    bool   isEMShowerPre();

    double getExpEnergy(int& dlayer) const;
    TVector3 getExpPos(int& dlayer) const;

    void FitProfile();
    void FitAxis(); 
    void AddShower(CRDEcalEDM::CRDCaloHit2DShower& shower);
    void MergeCluster(CRDEcalEDM::CRDCaloHit3DShower& clus);

    inline bool operator == ( const CRDEcalEDM::CRDCaloHit3DShower &x) const {
      TVector3 vec = this->getHitCenter();
      double En = this->getHitsE();
      return ( vec==x.getHitCenter() && En==x.getHitsE());
    }

    std::vector<CRDEcalEDM::CRDCaloHit2DShower> get2DShowers() const { return ShowerinLayer; }
    std::vector<edm4hep::ConstCalorimeterHit> getCaloHits() const { return CaloHits; }
    double getShowerMax() const {return showerMax; }
    double getChi2() const { return chi2; }
    TVector3 getAxis() const { return axis; }
    void getFitPars(double& _alpha, double& _beta ) const { _alpha = alpha; _beta = beta; }
    int getBeginningDlayer() const; 
    int getEndDlayer() const;     


    void setShowerVec( std::vector<CRDEcalEDM::CRDCaloHit2DShower>& _svec ) { ShowerinLayer = _svec;}
    void setCaloHits( std::vector<edm4hep::ConstCalorimeterHit>& _hits ) { CaloHits = _hits; }
    void Clear() { ShowerinLayer.clear(); CaloHits.clear(); showerEnergy=0; hitsEnergy=0; showerMax=0; chi2=0; alpha=0; beta=0; axis.SetXYZ(0,0,0); }

  private: 
    std::vector<CRDEcalEDM::CRDCaloHit2DShower> ShowerinLayer;
    std::vector<edm4hep::ConstCalorimeterHit> CaloHits;
    TVector3 axis;
    double showerEnergy;
    double hitsEnergy;
    double showerMax;
    double chi2;
    double alpha;
    double beta;

    TrackFitInEcal* track = new TrackFitInEcal(); 

  };
};
#endif
