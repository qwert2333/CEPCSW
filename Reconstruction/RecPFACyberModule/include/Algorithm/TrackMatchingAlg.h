#ifndef TRACKMATCHING_H
#define TRACKMATCHING_H

#include "CyberDataCol.h"
#include "Tools/Algorithm.h"
#include "TVector3.h"
using namespace Cyber;

class TrackMatchingAlg: public Cyber::Algorithm{

public:
  TrackMatchingAlg () {};
  ~TrackMatchingAlg() {};

  class Factory : public Cyber::AlgorithmFactory
  {
  public:
    Cyber::Algorithm* CreateAlgorithm() const{ return new TrackMatchingAlg(); }
  };

  StatusCode ReadSettings(Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode GetGoodTracks(std::vector<std::shared_ptr<Cyber::Track>>& rawTracks, std::vector<std::shared_ptr<Cyber::Track>>& goodTracks);
  StatusCode GetExtrpoECALPoints(const Cyber::Track* track, std::vector<TVector3>& extrapo_points);
  StatusCode CreateTrackAxis(vector<TVector3>& extrapo_points, std::vector<const Cyber::Calo1DCluster*>& localMaxVCol,
                             Cyber::CaloHalfCluster* t_track_axis);
  StatusCode GetNearby(const std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>* p_HalfClusterV, 
                    const std::vector<TVector3>& extrapo_points, 
                    std::vector<Cyber::CaloHalfCluster*>& t_nearbyHalfClusters, 
                    std::vector<const Cyber::Calo1DCluster*>& t_nearbyLocalMax);
  StatusCode LongiConeLinking(const std::vector<TVector3>& extrapo_points, 
                            std::vector<const Cyber::Calo1DCluster*>& nearbyLocalMax, 
                            std::vector<const Cyber::Calo1DCluster*>& cone_axis);
  bool isStopLinking(const std::vector<TVector3>& extrapo_points, 
                    const Cyber::Calo1DCluster* final_cone_hit);
  TVector2 GetProjectedRelR(const Cyber::Calo1DCluster* m_shower1, const Cyber::Calo1DCluster* m_shower2 );
  TVector2 GetProjectedAxis(const std::vector<TVector3>& extrapo_points, const Cyber::Calo1DCluster* m_shower);
  StatusCode CreatConeAxis(CyberDataCol& m_datacol, Cyber::Track* track, std::vector<Cyber::CaloHalfCluster*>& nearbyHalfClusters, 
                          std::vector<const Cyber::Calo1DCluster*>& cone_axis);

private:
  std::vector<Cyber::Track*> m_TrackCol;
  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>* p_HalfClusterV;
  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>* p_HalfClusterU;

  // std::vector<const Cyber::CaloHalfCluster*> m_trackAxisVCol;
  // std::vector<const Cyber::CaloHalfCluster*> m_trackAxisUCol;
  static bool compTrkP( std::shared_ptr<Cyber::Track> trk1, std::shared_ptr<Cyber::Track> trk2 ){
    return trk1->getMomentum() > trk2->getMomentum();
  }
};
#endif
