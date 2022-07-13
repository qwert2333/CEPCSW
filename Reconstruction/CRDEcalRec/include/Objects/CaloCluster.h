#ifndef CALO_CLUSTER_H
#define CALO_CLUSTER_H

#include "k4FWCore/DataHandle.h"
//#include "edm4hep/CalorimeterHit.h"
//#include "edm4hep/CalorimeterHitCollection.h"
//#include "edm4hep/SimCalorimeterHitCollection.h"
#include "Objects/TransShower.h"
#include "Objects/CaloHit.h"
#include "Objects/Track.h"
#include "Tools/TrackFitInEcal.h"

#include "TMath.h"
#include "TVector3.h"
#include "TString.h"

using namespace std;
namespace PandoraPlus{

  class CaloCluster{
  public:
    CaloCluster(){};
    ~CaloCluster() { Clear(); }

    void Clear();
    void Clean();
    void Check();

    const PandoraPlus::CaloHit* getClusterInitialHit() const; 
    TVector3 getHitCenter() const;
    TVector3 getShowerCenter() const;
    double getHitsE() const;
    double getShowerE() const;

    //void FitProfile();
    void FitAxis(); 
    void FitAxisHit(); 
    void addShower(const PandoraPlus::TransShower* _shower) { showers.push_back(_shower); };
    void addHit(const PandoraPlus::CaloHit* _hit) { hits.push_back(_hit); };
    void MergeCluster(const PandoraPlus::CaloCluster* clus);
    void IdentifyCluster(); 

    inline bool operator == ( const PandoraPlus::CaloCluster &x) const {
      //TVector3 vec = this->getHitCenter();
      //double En = this->getHitsE();
      //return ( vec==x.getHitCenter() && En==x.getHitsE());
      return (hits == x.hits);
    }

    std::vector<const PandoraPlus::TransShower*> getShowers() const { return showers; }
    std::vector<const PandoraPlus::CaloHit*> getCaloHits() const { return hits; }
    double getFitAlpha() const { return alpha; }
    double getFitBeta() const { return beta; }
    double getFitExpEn() const {return ExpEn; }
    double getShowerMax() const {return showerMax; }
    double getChi2() const { return chi2; }
    TVector3 getAxis(); 
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


    void setShowers( std::vector<const PandoraPlus::TransShower*>& _svec ) { showers = _svec;}
    void setCaloHits( std::vector<const PandoraPlus::CaloHit*> _hits ) { hits = _hits; }
    void setType( int _t ) { type = _t; }
    void Print() const; 

  private: 
    std::vector<const PandoraPlus::CaloHit*> hits;
    std::vector<const PandoraPlus::Track*> tracks;
    std::vector<const PandoraPlus::CaloCluster*> daughter_clusters;
    std::vector<const PandoraPlus::TransShower*> showers; //Specific for Bar Ecal.
    TVector3 axis;
    double showerMax;
    double chi2;
    double alpha;
    double beta;
    double ExpEn;
    int    type;  //0: MIP shower.  1: EM shower.  2: hadronic shower

    TrackFitInEcal trackFitter;

  };
};
#endif
