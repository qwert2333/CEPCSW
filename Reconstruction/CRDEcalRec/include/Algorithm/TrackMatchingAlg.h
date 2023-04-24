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

private:
  std::vector<PandoraPlus::Track*> m_TrackCol;
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>* p_HalfClusterV;
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>* p_HalfClusterU;

  std::vector<const PandoraPlus::CaloHalfCluster*> m_trackAxisVCol;
  std::vector<const PandoraPlus::CaloHalfCluster*> m_trackAxisUCol;


};
#endif
