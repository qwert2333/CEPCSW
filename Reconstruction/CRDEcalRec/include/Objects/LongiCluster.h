#ifndef CALO_LONGICLUS_H
#define CALO_LONGICLUS_H
#include <vector>
#include <map>
#include "Objects/Calo1DCluster.h"
#include "Tools/TrackFitInEcal.h"
#include "TVector3.h"

namespace PandoraPlus{

  class LongiCluster{
  public: 	
    LongiCluster(){};
    ~LongiCluster(){ Clear(); };

  void Clear() { axis.SetXYZ(0.,0.,0.); barShowerCol.clear(); }
  void Clean();
  void Check();

  inline bool operator == (const LongiCluster &x) const{
    return barShowerCol == x.barShowerCol;
  }  

  TVector3 getPos() const; 
  TVector3 getAxis() { FitAxis(); return axis; }
  double getEnergy() const; 
  int getSlayer() const; 
  std::vector<const PandoraPlus::Calo1DCluster*> getBarShowers() const { return barShowerCol; }
  std::vector<const PandoraPlus::Calo1DCluster*> getBarShowersInLayer(int _layer) const; 
  std::vector<const PandoraPlus::LongiCluster*>  getCousinClusters() const { return CousinClusters; }
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

  void FitAxis(); 
  void addBarShower( const PandoraPlus::Calo1DCluster* _shower, int option ); //option=0: do not merge the BarShowers in the same layer. 
                                                                              //option=1: merge the BarShowers in the same layer.  
  void addCousinCluster(const PandoraPlus::LongiCluster* _clus) { CousinClusters.push_back(_clus); }
  void SortBarShowersByLayer() { std::sort(barShowerCol.begin(), barShowerCol.end(), compLayer); }
  void setBarShowers( std::vector<const PandoraPlus::Calo1DCluster*> _barshwoers ) { barShowerCol = _barshwoers; }
  void setHoughPars(double _a, double _r) { Hough_alpha=_a; Hough_rho=_r; }
  void setIntercept(double _in) { Hough_intercept=_in; }
  void MergeCluster( const LongiCluster* clus ); 
  void RemoveShowers( std::vector<const PandoraPlus::Calo1DCluster*>& _showers ); 

  private:
    TVector3 axis; 
    double trk_dr; 
    double trk_dz;
    double Hough_alpha;
    double Hough_rho;
    double Hough_intercept;
    std::vector<const PandoraPlus::Calo1DCluster*> barShowerCol;
    TrackFitInEcal* track = new TrackFitInEcal();

    //const PandoraPlus::LongiCluster* ParentCluster; 
    std::vector<const PandoraPlus::LongiCluster*> CousinClusters; 


    static bool compLayer( const PandoraPlus::Calo1DCluster* hit1, const PandoraPlus::Calo1DCluster* hit2 ) 
      { return hit1->getDlayer() < hit2->getDlayer(); }

  };


};
#endif
