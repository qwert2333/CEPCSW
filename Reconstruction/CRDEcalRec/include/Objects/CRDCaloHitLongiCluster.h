#ifndef _CRD_CALOHIT_LONGICLUS_H
#define _CRD_CALOHIT_LONGICLUS_H
#include <vector>
#include <map>
#include "Objects/CRDCaloBarShower.h"
#include "Objects/TrackFitInEcal.h"
#include "TVector3.h"

namespace CRDEcalEDM{

  class CRDCaloHitLongiCluster{
  public: 	
    CRDCaloHitLongiCluster(){};
    ~CRDCaloHitLongiCluster(){};

  void Clear() { axis.SetXYZ(0.,0.,0.); barShowerCol.clear(); }
  inline bool operator == (const CRDCaloHitLongiCluster &x) const{
    return barShowerCol == x.barShowerCol;
  }  

  TVector3 getPos() const; 
  TVector3 getAxis() const {return axis; }
  int getSlayer() const; 
  std::vector<CRDEcalEDM::CRDCaloBarShower> getBarShowers() const { return barShowerCol; }

  int getBeginningDlayer() const; 
  int getEndDlayer() const; 
  TVector3 getExpPos(int& dlayer) const; 
  bool isContinue() const; 
  bool isContinueN(int n) const; 
  bool isSubset( CRDCaloHitLongiCluster& clus) const; 

  void FitAxis(); 
  void AddBarShower( CRDEcalEDM::CRDCaloBarShower _shower ) { barShowerCol.push_back(_shower); FitAxis(); }
  void SetBarShowers( std::vector<CRDEcalEDM::CRDCaloBarShower> _barshwoers ) { barShowerCol = _barshwoers; }

  private:
    TVector3 axis; 
    std::vector<CRDEcalEDM::CRDCaloBarShower> barShowerCol;
    TrackFitInEcal* track = new TrackFitInEcal();

  };


};
#endif
