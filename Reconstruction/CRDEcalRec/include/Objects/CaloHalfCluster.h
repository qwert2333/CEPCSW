#ifndef CALO_HALFCLUSTER_H
#define CALO_HALFCLUSTER_H

#include "Objects/CaloUnit.h"
#include "Objects/Calo1DCluster.h"
#include "Objects/Track.h"
#include "Tools/TrackFitInEcal.h"

namespace PandoraPlus {

  class CaloHalfCluster {
  public: 
    CaloHalfCluster() {};
    ~CaloHalfCluster() { Clear(); };

    void Clear();
    void Clean();
    void Check(); 

    inline bool operator == (const CaloHalfCluster &x) const{
      return m_1dclusters==x.getCluster();
    }
    CaloHalfCluster* Clone() const; 

    bool isNeighbor(const PandoraPlus::Calo1DCluster* m_1dcluster) const; 

    double getEnergy() const; 
    TVector3 getPos() const; 
    TVector3 getAxis() const { return axis; }
    TVector3 getEnergyCenter() const;
    int getSlayer() const { return m_slayer; }
    std::vector< std::vector<int> > getTowerID() const { return towerID; }
    double getHoughAlpha() const { return Hough_alpha; }
    double getHoughRho() const { return Hough_rho; }
    double getHoughIntercept() const { return Hough_intercept; }

    std::vector<const CaloUnit*> getBars() const;
    std::vector<const Calo1DCluster*> getCluster() const { return m_1dclusters;};
    std::vector<const Calo1DCluster*> getLocalMaxCol(std::string name) const;
    std::vector<const Calo1DCluster*> getAllLocalMaxCol() const;
    std::vector<const Calo1DCluster*> getClusterInLayer(int _layer) const;
    std::vector<const CaloHalfCluster*> getHalfClusterCol(std::string name) const;
    std::vector<const CaloHalfCluster*> getAllHalfClusterCol() const;
    std::map<std::string, std::vector<const PandoraPlus::CaloHalfCluster*> > getHalfClusterMap() const {return map_halfClusCol; }
    std::vector<const PandoraPlus::Track*> getAssotiatedTracks() const { return m_TrackCol; }

    int getBeginningDlayer() const;
    int getEndDlayer() const;
    bool isContinue() const;
    bool isContinueN(int n) const;
    bool isSubset(const CaloHalfCluster* clus) const;
    double OverlapRatioE( const CaloHalfCluster* clus ) const;

    void fitAxis( std::string name );
    void sortBarShowersByLayer() { std::sort(m_1dclusters.begin(), m_1dclusters.end(), compLayer); }
    void addUnit(const Calo1DCluster* _1dcluster);
    void setLocalMax( std::string name, std::vector<const Calo1DCluster*> _col) { map_localMax[name]=_col; }
    void setHalfClusters( std::string name, std::vector<const PandoraPlus::CaloHalfCluster*>& _cl) { map_halfClusCol[name]=_cl; }
    void addCousinCluster( const PandoraPlus::CaloHalfCluster* _cl ) { map_halfClusCol["CousinCluster"].push_back(_cl); }
    void setHoughPars(double _a, double _r) { Hough_alpha=_a; Hough_rho=_r; }
    void setIntercept(double _in) { Hough_intercept=_in; }
    void mergeHalfCluster( const CaloHalfCluster* clus );
    void addTowerID(int _m, int _p, int _s) { std::vector<int> id(3); id[0] = _m; id[1] = _p; id[2] = _s; towerID.push_back(id); }
    void addTowerID(std::vector<int> id) { towerID.push_back(id); }
    void addAssociatedTrack(const PandoraPlus::Track* _track){ m_TrackCol.push_back(_track); }

  private:
    std::vector< std::vector<int> > towerID; //[module, part, stave]
    int m_slayer;
    TVector3 axis;
    double trk_dr;
    double trk_dz;
    double Hough_alpha;
    double Hough_rho;
    double Hough_intercept;
    std::vector<const Calo1DCluster*> m_1dclusters; 
    std::map<std::string, std::vector<const PandoraPlus::Calo1DCluster*> > map_localMax;
    std::map<std::string, std::vector<const PandoraPlus::CaloHalfCluster*> > map_halfClusCol;  
    std::vector<const PandoraPlus::Track*> m_TrackCol;

    TrackFitInEcal* track = new TrackFitInEcal();

    static bool compLayer( const PandoraPlus::Calo1DCluster* hit1, const PandoraPlus::Calo1DCluster* hit2 )
      { return hit1->getDlayer() < hit2->getDlayer(); }

  };

};
#endif
