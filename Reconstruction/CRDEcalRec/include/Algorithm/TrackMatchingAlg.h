#ifndef TRACKMATCHING_H
#define TRACKMATCHING_H

#include "PandoraPlusDataCol.h"
#include "Tools/Algorithm.h"
#include "TVector3.h"
using namespace PandoraPlus;

class TrackMatchingAlg: public PandoraPlus::Algorithm{

public:
  TrackMatchingAlg () {};
  ~TrackMatchingAlg() {};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public:
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new TrackMatchingAlg(); }
  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();


  StatusCode GetExtrpoECALPoints(const PandoraPlus::Track* track, std::vector<TVector3>& extrapo_points);
  StatusCode CreateTrackAxis(vector<TVector3>& extrapo_points, std::vector<const PandoraPlus::Calo1DCluster*>& localMaxVCol,
                             PandoraPlus::CaloHalfCluster* t_track_axis);
  StatusCode GetNearby(const std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>* p_HalfClusterV, 
                    const std::vector<TVector3>& extrapo_points, 
                    std::vector<PandoraPlus::CaloHalfCluster*>& t_nearbyHalfClusters, 
                    std::vector<const PandoraPlus::Calo1DCluster*>& t_nearbyLocalMax);
  StatusCode LongiConeLinking(const std::vector<TVector3>& extrapo_points, 
                            std::vector<const PandoraPlus::Calo1DCluster*>& nearbyLocalMax, 
                            std::vector<const PandoraPlus::Calo1DCluster*>& cone_axis);
  bool isStopLinking(const std::vector<TVector3>& extrapo_points, 
                    const PandoraPlus::Calo1DCluster* final_cone_hit);
  TVector2 GetProjectedRelR(const PandoraPlus::Calo1DCluster* m_shower1, const PandoraPlus::Calo1DCluster* m_shower2 );
  TVector2 GetProjectedAxis(const std::vector<TVector3>& extrapo_points, const PandoraPlus::Calo1DCluster* m_shower);
  StatusCode CreatConeAxis(PandoraPlusDataCol& m_datacol, PandoraPlus::Track* track, std::vector<PandoraPlus::CaloHalfCluster*>& nearbyHalfClusters, 
                          std::vector<const PandoraPlus::Calo1DCluster*>& cone_axis);

private:
  std::vector<PandoraPlus::Track*> m_TrackCol;
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>* p_HalfClusterV;
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>* p_HalfClusterU;

  // std::vector<const PandoraPlus::CaloHalfCluster*> m_trackAxisVCol;
  // std::vector<const PandoraPlus::CaloHalfCluster*> m_trackAxisUCol;


};
#endif
