#ifndef _CRD_CALOHIT_3DSHOWER_H
#define _CRD_CALOHIT_3DSHOWER_H

#include "k4FWCore/DataHandle.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "Objects/CRDCaloHit2DShower.h"
#include "Objects/TrackFitInEcal.h"

#include "TMath.h"
#include "TVector3.h"
#include "TString.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace std;

namespace TMVA {
  class Reader;
}

namespace CRDEcalEDM{

  class CRDCaloHit3DCluster{
  public:
    CRDCaloHit3DCluster();
    //~CRDCaloHit3DCluster(){ delete track;}

    int getclusterIDlayer() const;  //swz_1 
    double getlateral() const;

    edm4hep::ConstCalorimeterHit getClusterInitialHit() const; 
    TVector3 getHitCenter() const;
    TVector3 getShowerCenter() const;
    double getHitsE() const;
    double getShowerE() const;

    double getExpEnergy(int& dlayer) const;
    TVector3 getExpPos(int& dlayer) const;

    void FitProfile();
    void FitAxis(); 
    void AddShower(CRDEcalEDM::CRDCaloHit2DShower& shower);
    void MergeCluster(CRDEcalEDM::CRDCaloHit3DCluster& clus);
    void IdentifyCluster(); 

    inline bool operator == ( const CRDEcalEDM::CRDCaloHit3DCluster &x) const {
      TVector3 vec = this->getHitCenter();
      double En = this->getHitsE();
      return ( vec==x.getHitCenter() && En==x.getHitsE());
    }

    std::vector<CRDEcalEDM::CRDCaloHit2DShower> get2DShowers() const { return ShowerinLayer; }
    std::vector<edm4hep::ConstCalorimeterHit> getCaloHits() const { return CaloHits; }
    float getFitAlpha() const { return alpha; }
    float getFitBeta() const { return beta; }
    double getFitExpEn() const {return ExpEn; }
    float getChi2() const { return chi2; }
    TVector3 getAxis() const { return axis; }
    void getFitPars(float& _alpha, float& _beta ) const { _alpha = alpha; _beta = beta; }
    int getBeginningDlayer() const; 
    int getEndDlayer() const;

    std::vector<double> getEnInLayer() const; 
    int getMaxELayer() const; 
    double getAveE() const; 
    double getStdDevE() const; 

    std::vector<double> getClusterWidth() const;
    double getMaxWidth() const; 
    int getType() const { return type; }
    double getBDTValue() const { return ID_BDT; }

    void setShowerVec( std::vector<CRDEcalEDM::CRDCaloHit2DShower>& _svec ) { ShowerinLayer = _svec;}
    void setCaloHits( std::vector<edm4hep::ConstCalorimeterHit>& _hits ) { CaloHits = _hits; }
    void setType( int _t ) { type = _t; }
    void Clear() { ShowerinLayer.clear(); CaloHits.clear(); showerEnergy=0; hitsEnergy=0; chi2=0; alpha=0; beta=0; ExpEn=0; axis.SetXYZ(0,0,0); type=-1; ID_BDT=-999; }

  private: 
    std::vector<CRDEcalEDM::CRDCaloHit2DShower> ShowerinLayer;
    std::vector<edm4hep::ConstCalorimeterHit> CaloHits;
    TVector3 axis;
    double showerEnergy;
    double hitsEnergy;
    float chi2;
    float alpha;
    float beta;
    int    type;  //0: MIP shower.  1: EM shower.  2: hadronic shower
    double ID_BDT;
    double ExpEn; 

    TrackFitInEcal* track = new TrackFitInEcal(); 

    float Lend; 
    float LmaxWidth; 
    float LmaxE;
    float maxEnergy; 
    float maxWidth; 
    //TMVA::Reader *reader = new TMVA::Reader("Silent");  
  };
};
#endif
