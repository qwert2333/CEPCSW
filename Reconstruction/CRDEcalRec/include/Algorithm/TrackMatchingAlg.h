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


  StatusCode GetExtrpoPoints(const PandoraPlus::Track* track, std::vector<TVector3>& extrapo_points);
  StatusCode CreateTrackAxis(vector<TVector3>& extrapo_points, std::vector<const PandoraPlus::Calo1DCluster*>& localMaxVCol,
                             PandoraPlus::CaloHalfCluster* t_track_axis);
  float GetBarHalfLength(const PandoraPlus::Calo1DCluster* local_max);

private:
  std::vector<PandoraPlus::Track*> m_TrackCol;
  std::vector<PandoraPlus::CaloHalfCluster*> m_HalfClusterV;
  std::vector<PandoraPlus::CaloHalfCluster*> m_HalfClusterU;

  std::vector<const PandoraPlus::CaloHalfCluster*> m_trackAxisVCol;
  std::vector<const PandoraPlus::CaloHalfCluster*> m_trackAxisUCol;


};
#endif