#ifndef CALO_3DCLUSTER_H
#define CALO_3DCLUSTER_H

#include "Objects/Calo2DCluster.h"
#include "Objects/Track.h"
#include "Objects/CaloHalfCluster.h"
//#include "Objects/CaloTower.h"
#include "Tools/TrackFitInEcal.h"

namespace PandoraPlus {

  class Calo3DCluster {
  public: 
    Calo3DCluster() {};
    ~Calo3DCluster() { Clear(); };

    void Clear();
    void Clean();
    void Check(); 
    void Clear2DClusters() { m_2dclusters.clear(); };
    void CleanLongiClusters();
    std::shared_ptr<PandoraPlus::Calo3DCluster> Clone() const;

    inline bool operator == (const Calo3DCluster &x) const{
      return ( m_2dclusters == x.getCluster()  );
    }

    std::vector< std::vector<int> > getTowerID() const { return towerID; }
    bool isNeighbor(const PandoraPlus::Calo2DCluster* m_2dcluster) const; 
    std::vector<const PandoraPlus::CaloHit*> getCaloHits() const { return hits; }
    std::vector<const Calo2DCluster*> getCluster() const { return m_2dclusters; }
    std::vector<const Calo1DCluster*> getLocalMaxUCol(std::string name) const;
    std::vector<const Calo1DCluster*> getLocalMaxVCol(std::string name) const;
    std::map<std::string, std::vector<const PandoraPlus::Calo1DCluster*> >  getLocalMaxUMap() const { return map_localMaxU; }
    std::map<std::string, std::vector<const PandoraPlus::Calo1DCluster*> >  getLocalMaxVMap() const { return map_localMaxV; }
    std::vector<const CaloHalfCluster*> getHalfClusterUCol(std::string name) const;
    std::vector<const CaloHalfCluster*> getHalfClusterVCol(std::string name) const; 
    std::map<std::string, std::vector<const PandoraPlus::CaloHalfCluster*> > getHalfClusterUMap() const { return map_halfClusUCol; }
    std::map<std::string, std::vector<const PandoraPlus::CaloHalfCluster*> > getHalfClusterVMap() const { return map_halfClusVCol; }
    std::vector<const PandoraPlus::Track*> getAssociatedTracks() const { return m_TrackCol; }
    std::vector< std::pair<edm4hep::MCParticle, float> > getLinkedMCP() const { return MCParticleWeight; }
    std::vector< std::pair<edm4hep::MCParticle, float> > getLinkedMCPfromHFCluster(std::string name);
    std::vector< std::pair<edm4hep::MCParticle, float> > getLinkedMCPfromHit();
    edm4hep::MCParticle getLeadingMCP() const;
    float getLeadingMCPweight() const;

    std::vector<const Calo3DCluster*> getTowers() const {return m_towers; }
    std::vector<const CaloUnit*> getBars() const;
    double getEnergy() const; 
    double getHitsE() const;
    double getLongiE() const;
    TVector3 getHitCenter() const;
    TVector3 getShowerCenter() const; 
    TVector3 getAxis() const { return axis; }
    int getBeginningDlayer() const;
    int getEndDlayer() const;

    void setCaloHits( std::vector<const PandoraPlus::CaloHit*> _hits ) { hits = _hits; }
    void setTowers(std::vector<const Calo3DCluster*> _t) { m_towers = _t; }
    void setHalfClusters( std::string name1, std::vector<const CaloHalfCluster*>& _clU, 
                           std::string name2, std::vector<const CaloHalfCluster*>& _clV )
    { map_halfClusUCol[name1]=_clU; map_halfClusVCol[name2]=_clV;  }

    void setLocalMax( std::string name1, std::vector<const Calo1DCluster*>& _colU,
                      std::string name2, std::vector<const Calo1DCluster*>& _colV )
    { map_localMaxU[name1]=_colU; map_localMaxV[name2]=_colV; }
    void setClusters(std::vector<const Calo2DCluster*> _2dcol) { m_2dclusters = _2dcol; }
    void addTowerID(int _m, int _p, int _s) { std::vector<int> id(3); id[0] = _m; id[1] = _p; id[2] = _s; towerID.push_back(id); }
    void addTowerID( std::vector<int> id ) { towerID.push_back(id); }

    void addUnit(const Calo2DCluster* _2dcluster);
    void addHit(const PandoraPlus::CaloHit* _hit) { hits.push_back(_hit); };
    void addTower( const Calo3DCluster* _tower ) { m_towers.push_back(_tower); }
    void addHalfClusterU( std::string name, const CaloHalfCluster* _clU ) { map_halfClusUCol[name].push_back(_clU); }
    void addHalfClusterV( std::string name, const CaloHalfCluster* _clV ) { map_halfClusVCol[name].push_back(_clV); }
    void addLocalMaxU( std::string name, const Calo1DCluster* _shU ) { map_localMaxU[name].push_back(_shU); }
    void addLocalMaxV( std::string name, const Calo1DCluster* _shV ) { map_localMaxV[name].push_back(_shV); }
    void addAssociatedTrack(const PandoraPlus::Track* _track){ m_TrackCol.push_back(_track); }
    void addLinkedMCP( std::pair<edm4hep::MCParticle, float> _pair ) { MCParticleWeight.push_back(_pair); }
    void setLinkedMCP( std::vector<std::pair<edm4hep::MCParticle, float>> _pairVec ) { MCParticleWeight.clear(); MCParticleWeight = _pairVec; }
    void mergeCluster( const PandoraPlus::Calo3DCluster* _clus );

    //void FitProfile();
    void FitAxis();
    //void FitAxisHit();

    bool isHCALNeighbor(const PandoraPlus::CaloHit* m_hit) const;

  private:
    std::vector<const PandoraPlus::CaloHit*> hits; 
    std::vector<const PandoraPlus::Calo2DCluster*> m_2dclusters;
    std::vector<const Calo3DCluster*> m_towers;
    std::map<std::string, std::vector<const PandoraPlus::Calo1DCluster*> > map_localMaxU;
    std::map<std::string, std::vector<const PandoraPlus::Calo1DCluster*> > map_localMaxV;
    std::map<std::string, std::vector<const PandoraPlus::CaloHalfCluster*> > map_halfClusUCol;
    std::map<std::string, std::vector<const PandoraPlus::CaloHalfCluster*> > map_halfClusVCol;
    std::vector<const PandoraPlus::Track*> m_TrackCol;
    std::vector< std::pair<edm4hep::MCParticle, float> > MCParticleWeight;

    std::vector< std::vector<int> > towerID; //[module, part, stave]
    TVector3 axis;
    double chi2;
    double alpha;
    double beta;
    int    type;  //0: MIP shower.  1: EM shower.  2: hadronic shower
    TrackFitInEcal trackFitter;

  };

};
#endif
