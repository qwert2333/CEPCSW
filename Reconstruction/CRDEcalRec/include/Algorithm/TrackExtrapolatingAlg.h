#ifndef _TRACKEXTRAPOLATING_ALG_H
#define _TRACKEXTRAPOLATING_ALG_H

#include "PandoraPlusDataCol.h"
#include "Tools/Algorithm.h"
#include "Objects/Track.h"

using namespace PandoraPlus;
class TrackExtrapolatingAlg: public PandoraPlus::Algorithm{
public: 

  TrackExtrapolatingAlg(){};
  ~TrackExtrapolatingAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new TrackExtrapolatingAlg(); } 

  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize();
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  //Self defined algorithms
  // StatusCode SelfAlg1(); 

  // Get normal vectors of planes in each module
  StatusCode GetPlaneNormalVector(std::vector<TVector2> & normal_vectors);
  // Get points in each plane of layer
  StatusCode GetLayerPoints(const std::vector<TVector2> & normal_vectors, 
                          std::vector<std::vector<TVector2>> & layer_points);
  // If the track reach barrel ECAL
  bool IsReachECAL(PandoraPlus::Track * track);
  // Get track state at calorimeter
  StatusCode GetECALTrackState(PandoraPlus::Track * track, 
                             PandoraPlus::TrackState & ECAL_trk_state);
  // get extrapolated points
  StatusCode ExtrapolateByLayer(const std::vector<TVector2> & normal_vectors,
                              const std::vector<std::vector<TVector2>> & layer_points, 
                              const PandoraPlus::TrackState & ECAL_trk_state, 
                              PandoraPlus::Track* p_track);
  // Get the radius rho 
  float GetRho(const PandoraPlus::TrackState & trk_state);
  // Get coordinates of the center of the circle
  TVector2 GetCenterOfCircle(const PandoraPlus::TrackState & trk_state, const float & rho);
  // phase from center to reference point
  float GetRefAlpha0(const PandoraPlus::TrackState & trk_state, const TVector2 & center);
  // If the charged particle return back 
  bool IsReturn(float rho, TVector2 & center);

private: 

};
#endif
