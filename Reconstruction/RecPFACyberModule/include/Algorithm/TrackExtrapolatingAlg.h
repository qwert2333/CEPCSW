#ifndef _TRACKEXTRAPOLATING_ALG_H
#define _TRACKEXTRAPOLATING_ALG_H

#include "CyberDataCol.h"
#include "Tools/Algorithm.h"
#include "Objects/Track.h"

using namespace Cyber;
class TrackExtrapolatingAlg: public Cyber::Algorithm{
public: 

  TrackExtrapolatingAlg(){};
  ~TrackExtrapolatingAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
  public: 
    Cyber::Algorithm* CreateAlgorithm() const{ return new TrackExtrapolatingAlg(); } 

  };

  StatusCode ReadSettings(Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();


  // Get points in each plane of layer
  StatusCode GetLayerRadius(std::vector<float> & ECAL_layer_radius, std::vector<float> & HCAL_layer_radius);
  StatusCode GetLayerZ(std::vector<float> & ECAL_layer_z, std::vector<float> & HCAL_layer_z);
  // If the track reach barrel ECAL
  bool IsReachECAL(Cyber::Track * track);
  // Get track state at calorimeter
  StatusCode GetTrackStateAtCalo(Cyber::Track * track, 
                           Cyber::TrackState & trk_state_at_calo);
  // get extrapolated points
  StatusCode Extrapolate( const std::vector<float> & ECAL_layer_radius, const std::vector<float> & ECAL_layer_z, 
                          const std::vector<float> & HCAL_layer_radius, const std::vector<float> & HCAL_layer_z,
                          const Cyber::TrackState & CALO_trk_state, Cyber::Track* p_track);
  // Get the radius rho 
  float GetRho(const Cyber::TrackState & trk_state);
  // Get coordinates of the center of the circle
  TVector2 GetCenterOfCircle(const Cyber::TrackState & trk_state, const float & rho);
  // phase from center to reference point
  float GetRefAlpha0(const Cyber::TrackState & trk_state, const TVector2 & center);
  // If the charged particle return back 
  bool IsReturn(float rho, TVector2 & center);

  StatusCode GetDeltaPhi( float rho, TVector2 center, float alpha0, 
                          const std::vector<float> & layer_radius, const std::vector<float> & layer_z, 
                          const Cyber::TrackState & CALO_trk_state,
                          vector<float> & barrel_delta_phi, vector<float> & endcap_delta_phi);
  std::vector<TVector3> GetExtrapoPoints(std::string calo_name, 
                                         float rho, TVector2 center, float alpha0, 
                                         const Cyber::TrackState & CALO_trk_state,
                                         const std::vector<float>& barrel_delta_phi, const std::vector<float>& endcap_delta_phi);

    // Get phi0 of extrapolated points. Note that this phi0 is not same as the definition of the phi0 in TrackState, but will be stored in TrackState
  float GetExtrapolatedPhi0(float Kappa, float ECAL_phi0, TVector2 center, TVector3 ext_point);
  // To sort the extrapolatedpoints, define the following comparison function
  static bool SortByPhi0(const Cyber::TrackState& trk_state1, const Cyber::TrackState& trk_state2 ) 
  { return TMath::Abs(trk_state1.phi0) < TMath::Abs(trk_state2.phi0); }
  
private: 

};
#endif
