#ifndef CALO_LONGICLUS_H
#define CALO_LONGICLUS_H
#include <vector>
#include <map>
#include "Objects/CaloBarShower.h"
//#include "Objects/TrackFitInEcal.h"
#include "TVector3.h"

namespace PandoraPlus{

  class LongiCluster{
  public: 	
    LongiCluster(){};
    ~LongiCluster(){};

  void Clear() { axis.SetXYZ(0.,0.,0.); barShowerCol.clear(); }
  void Clean();
  void Check();

  inline bool operator == (const LongiCluster &x) const{
    return barShowerCol == x.barShowerCol;
  }  

  TVector3 getPos() const; 
  //TVector3 getAxis() { FitAxis(); return axis; }
  double getE() const; 
  int getSlayer() const; 
  std::vector<const PandoraPlus::CaloBarShower*> getBarShowers() const { return barShowerCol; }
  std::vector<const PandoraPlus::CaloBarShower*> getBarShowersInLayer(int _layer) const; 
  double getHoughAlpha() const { return Hough_alpha; }
  double getHoughRho() const { return Hough_rho; }
  double getHoughIntercept() const { return Hough_intercept; }
  //double getFitTrkDr() { FitAxis(); return trk_dr; }
  //double getFitTrkDz() { FitAxis(); return trk_dz; }

  int getBeginningDlayer() const; 
  int getEndDlayer() const; 
  bool isContinue() const; 
  bool isContinueN(int n) const; 
  bool isSubset(const LongiCluster* clus) const; 
  double OverlapRatioE( const LongiCluster* clus ) const;

  //void FitAxis(); 
  void addBarShower( const PandoraPlus::CaloBarShower* _shower ) { barShowerCol.push_back(_shower); /*FitAxis();*/ }
  void SortBarShowersByLayer() { std::sort(barShowerCol.begin(), barShowerCol.end(), compLayer); }
  void setBarShowers( std::vector<const PandoraPlus::CaloBarShower*> _barshwoers ) { barShowerCol = _barshwoers; }
  void setHoughPars(double _a, double _r) { Hough_alpha=_a; Hough_rho=_r; }
  void setIntercept(double _in) { Hough_intercept=_in; }
  void MergeCluster( const LongiCluster* clus ); 
  void RemoveShowers( std::vector<const PandoraPlus::CaloBarShower*>& _showers ); 

  private:
    TVector3 axis; 
    double trk_dr; 
    double trk_dz;
    double Hough_alpha;
    double Hough_rho;
    double Hough_intercept;
    std::vector<const PandoraPlus::CaloBarShower*> barShowerCol;
    //TrackFitInEcal* track = new TrackFitInEcal();

    static bool compLayer( const PandoraPlus::CaloBarShower* hit1, const PandoraPlus::CaloBarShower* hit2 ) 
      { return hit1->getDlayer() < hit2->getDlayer(); }

  };


};
#endif
