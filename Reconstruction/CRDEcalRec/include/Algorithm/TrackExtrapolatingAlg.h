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
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
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
  StatusCode GetTrackPoints(const PandoraPlus::TrackState & ECAL_trk_state, PandoraPlus::Track* p_track);

    // Get phi0 of extrapolated points. Note that this phi0 is not same as the definition of the phi0 in TrackState, but will be stored in TrackState
  float GetExtrapolatedPhi0(float Kappa, float ECAL_phi0, TVector2 center, TVector3 ext_point);
  // To sort the extrapolatedpoints, define the following comparison function
  static bool SortByPhi0(const PandoraPlus::TrackState& trk_state1, const PandoraPlus::TrackState& trk_state2 ) 
  { return TMath::Abs(trk_state1.phi0) < TMath::Abs(trk_state2.phi0); }
  
private: 

};
#endif
