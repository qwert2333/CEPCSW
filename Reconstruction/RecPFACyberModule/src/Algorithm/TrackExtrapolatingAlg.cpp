#ifndef _TRACKEXTRAPOLATING_ALG_C
#define _TRACKEXTRAPOLATING_ALG_C

#include "TVector2.h"

#include "Algorithm/TrackExtrapolatingAlg.h"
#include "Objects/Track.h"
#include "Objects/TrackState.h"

using namespace TMath;
using namespace std;

StatusCode TrackExtrapolatingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  // ECAL parameters
  if(settings.map_floatPars.find("ECAL_innermost_distance")==settings.map_floatPars.end()) 
    settings.map_floatPars["ECAL_innermost_distance"] = 1830;
  if(settings.map_floatPars.find("ECAL_outermost_distance")==settings.map_floatPars.end()) 
    settings.map_floatPars["ECAL_outermost_distance"] = 1830+300;
  if(settings.map_intPars.find("ECAL_Nlayers")==settings.map_intPars.end()) 
    settings.map_intPars["ECAL_Nlayers"] = 28;
  if(settings.map_floatPars.find("ECAL_layer_width")==settings.map_floatPars.end())
    settings.map_floatPars["ECAL_layer_width"] = Cyber::CaloUnit::barsize;
  if(settings.map_floatPars.find("ECAL_half_length")==settings.map_floatPars.end())
    settings.map_floatPars["ECAL_half_length"] = 2900;
  if(settings.map_floatPars.find("ECAL_endcap_zmin")==settings.map_floatPars.end())
    settings.map_floatPars["ECAL_endcap_zmin"] = 2930;
  if(settings.map_floatPars.find("ECAL_endcap_zmax")==settings.map_floatPars.end())
    settings.map_floatPars["ECAL_endcap_zmax"] = 3230;  // 2930+300

  // HCAL parameters
  if(settings.map_floatPars.find("HCAL_innermost_distance")==settings.map_floatPars.end()) 
    settings.map_floatPars["HCAL_innermost_distance"] = 2140;
  if(settings.map_floatPars.find("HCAL_outermost_distance")==settings.map_floatPars.end()) 
    settings.map_floatPars["HCAL_outermost_distance"] = 3455;
  if(settings.map_intPars.find("HCAL_Nlayers")==settings.map_intPars.end()) 
    settings.map_intPars["HCAL_Nlayers"] = 48;
  if(settings.map_floatPars.find("HCAL_layer_width")==settings.map_floatPars.end())
    settings.map_floatPars["HCAL_layer_width"] = 26.61;
  if(settings.map_floatPars.find("HCAL_sensitive_distance")==settings.map_floatPars.end())
    settings.map_floatPars["HCAL_sensitive_distance"] = 22.81;  // distance between sensitive material and front face of each layer
  if(settings.map_floatPars.find("HCAL_half_length")==settings.map_floatPars.end())
    settings.map_floatPars["HCAL_half_length"] = 3230;
  if(settings.map_floatPars.find("HCAL_endcap_zmin")==settings.map_floatPars.end())
    settings.map_floatPars["HCAL_endcap_zmin"] = 3260;
  if(settings.map_floatPars.find("HCAL_endcap_zmax")==settings.map_floatPars.end())
    settings.map_floatPars["HCAL_endcap_zmax"] = 4575;  // 3260 + (3455-2140)

  if(settings.map_intPars.find("Nmodule")==settings.map_intPars.end()) 
    settings.map_intPars["Nmodule"] = Cyber::CaloUnit::Nmodule;
  if(settings.map_floatPars.find("B_field")==settings.map_floatPars.end())
    settings.map_floatPars["B_field"] = 3.0;

  if(settings.map_intPars.find("Input_track")==settings.map_intPars.end())
    settings.map_intPars["Input_track"] = 0;    // 0: reconstructed tracks.  1: MC particle track

  return StatusCode::SUCCESS;
};


StatusCode TrackExtrapolatingAlg::Initialize( CyberDataCol& m_datacol ){
  std::cout<<"Initialize TrackExtrapolatingAlg"<<std::endl;

  return StatusCode::SUCCESS;
};


StatusCode TrackExtrapolatingAlg::RunAlgorithm( CyberDataCol& m_datacol ){
  std::vector<std::shared_ptr<Cyber::Track>>* p_tracks = &(m_datacol.TrackCol);
  // (Barrel) Layer radius from ECAL(HCAL) innermost distance to outermost distance; spacing: ECAL_layer_width/2 
  std::vector<float> ECAL_layer_radius;  
  std::vector<float> HCAL_layer_radius;
  GetLayerRadius(ECAL_layer_radius, HCAL_layer_radius);
  // (Endcap) Layer z from ECAL innermost distance to outermost distance; spacing: ECAL_layer_width/2
  std::vector<float> ECAL_layer_z;  
  std::vector<float> HCAL_layer_z;
  GetLayerZ(ECAL_layer_z, HCAL_layer_z);

  for(int itrk=0; itrk<p_tracks->size(); itrk++){
    // Only tracks that reach ECAL should be processed.
    if(!IsReachECAL( p_tracks->at(itrk).get() )) continue;
    // get track state at calorimeter
    Cyber::TrackState CALO_trk_state;
    GetTrackStateAtCalo(p_tracks->at(itrk).get(), CALO_trk_state);
    Extrapolate(ECAL_layer_radius, ECAL_layer_z, HCAL_layer_radius, HCAL_layer_z, CALO_trk_state, p_tracks->at(itrk).get());
  } // end loop tracks


  p_tracks = nullptr;
  return StatusCode::SUCCESS;
}; // RunAlgorithm end


StatusCode TrackExtrapolatingAlg::ClearAlgorithm(){
  std::cout<<"End run TrackExtrapolatingAlg. Clean it."<<std::endl;

  return StatusCode::SUCCESS;
};


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
StatusCode TrackExtrapolatingAlg::GetLayerRadius(std::vector<float> & ECAL_layer_radius, std::vector<float> & HCAL_layer_radius){
  // ECAL
  float tmp_ecal_radius = settings.map_floatPars["ECAL_innermost_distance"] + 0.5*settings.map_floatPars["ECAL_layer_width"];
  while(tmp_ecal_radius < settings.map_floatPars["ECAL_outermost_distance"]){
    ECAL_layer_radius.push_back(tmp_ecal_radius);
    tmp_ecal_radius += 0.5 * settings.map_floatPars["ECAL_layer_width"];
  }
  
  // HCAL
  float tmp_hcal_radius = settings.map_floatPars["HCAL_innermost_distance"] + 0.5*settings.map_floatPars["HCAL_layer_width"];
  while(tmp_hcal_radius < settings.map_floatPars["HCAL_outermost_distance"]){
    HCAL_layer_radius.push_back(tmp_hcal_radius);
    tmp_hcal_radius += 0.5 * settings.map_floatPars["HCAL_layer_width"];
  }

  return StatusCode::SUCCESS;
}

StatusCode TrackExtrapolatingAlg::GetLayerZ(std::vector<float> & ECAL_layer_z, std::vector<float> & HCAL_layer_z){
  // ECAL
  float tmp_ecal_z = settings.map_floatPars["ECAL_endcap_zmin"] + 0.5*settings.map_floatPars["ECAL_layer_width"];
  while(tmp_ecal_z < settings.map_floatPars["ECAL_endcap_zmax"]){
    ECAL_layer_z.push_back(tmp_ecal_z);
    tmp_ecal_z += 0.5 * settings.map_floatPars["ECAL_layer_width"];
  }

  // HCAL
  float tmp_hcal_z = settings.map_floatPars["HCAL_endcap_zmin"] + 0.5*settings.map_floatPars["HCAL_layer_width"];
  while(tmp_hcal_z < settings.map_floatPars["HCAL_endcap_zmax"]){
    HCAL_layer_z.push_back(tmp_hcal_z);
    tmp_hcal_z += 0.5 * settings.map_floatPars["HCAL_layer_width"];
  }

  return StatusCode::SUCCESS;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
bool TrackExtrapolatingAlg::IsReachECAL(Cyber::Track * track){
  if(settings.map_intPars["Input_track"] == 0){
    // The track is reconstructed in ECAL. If the track reach ECAL, it should have track state at calorimeter
    std::vector<TrackState> input_trackstates = track->getTrackStates("Input");
    int count=0;
    TVector3 t_vec;
    for(int i=0; i<input_trackstates.size(); i++){
      if(input_trackstates[i].location==Cyber::TrackState::AtCalorimeter){
        count++;
        t_vec = input_trackstates[i].referencePoint;
        break;
      }
    }
    if(count==0){
      // The track has no track state at calorimeter
      return false;
    } 
    // if( Abs(Abs(t_vec.Z())-settings.map_floatPars["ECAL_half_length"]) < 1.0 ){
    //   // The track escape from endcap
    //   return false;
    // } 
    
    return true;
  }
  else if(settings.map_intPars["Input_track"] == 1){
    // The track is from MC particle as ideal helix. 
    // The pT should large enough to reach ECAL. The pz should not be so large that it escape from endcap
    std::vector<TrackState> input_trackstates = track->getTrackStates("Input");
    if(input_trackstates.size()==0){
      std::cout << "Error! No track state!" << std::endl;
      return false;
    }

    TrackState IP_trk_state;
    for(int i=0; i<input_trackstates.size(); i++){
      if(input_trackstates[i].location==Cyber::TrackState::AtIP)
        IP_trk_state = input_trackstates[i];
      break;
    }

    TVector3 ref_point = IP_trk_state.referencePoint;
    double rho = GetRho(IP_trk_state);
    double r_max = TMath::Sqrt(ref_point.X()*ref_point.X() + ref_point.Y()*ref_point.Y()) + rho*2;

    // if(r_max<settings.map_floatPars["ECAL_innermost_distance"]){ return false; }

    return true;
  }
  else{
    std::cout << "Error, wrong source of input tracks for TrackExtrapolatingAlg!" << std:: endl;
    return false;
  }
  
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
StatusCode TrackExtrapolatingAlg::GetTrackStateAtCalo(Cyber::Track * track, 
                                                Cyber::TrackState & trk_state_at_calo){
  std::vector<TrackState> input_trackstates = track->getTrackStates("Input");

  if(settings.map_intPars["Input_track"] == 0){
    for(int its=0; its<input_trackstates.size(); its++){
      if(input_trackstates[its].location==Cyber::TrackState::AtCalorimeter){
        trk_state_at_calo=input_trackstates[its];
        break;
      }
    }
  }
  else if((settings.map_intPars["Input_track"] == 1)){
    for(int its=0; its<input_trackstates.size(); its++){
      if(input_trackstates[its].location==Cyber::TrackState::AtIP){
        trk_state_at_calo=input_trackstates[its];
        break;
      }
    }
  }
  else{
    std::cout << "Error, wrong source of input tracks for TrackExtrapolatingAlg!" << std:: endl;
  }
  
  return StatusCode::SUCCESS;
  
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
StatusCode TrackExtrapolatingAlg::Extrapolate( const std::vector<float> & ECAL_layer_radius, const std::vector<float> & ECAL_layer_z, 
                          const std::vector<float> & HCAL_layer_radius, const std::vector<float> & HCAL_layer_z,
                          const Cyber::TrackState & CALO_trk_state, Cyber::Track* p_track){
  // Extrapolate points to circles with specific radius
  float rho = GetRho(CALO_trk_state);
  TVector2 center = GetCenterOfCircle(CALO_trk_state, rho);
  float alpha0 = GetRefAlpha0(CALO_trk_state, center); 

  std::vector<float> ECAL_barrel_delta_phi, ECAL_endcap_delta_phi;
  GetDeltaPhi(rho, center, alpha0, ECAL_layer_radius, ECAL_layer_z, CALO_trk_state, ECAL_barrel_delta_phi, ECAL_endcap_delta_phi);
  std::vector<float> HCAL_barrel_delta_phi, HCAL_endcap_delta_phi;
  GetDeltaPhi(rho, center, alpha0, HCAL_layer_radius, HCAL_layer_z, CALO_trk_state, HCAL_barrel_delta_phi, HCAL_endcap_delta_phi);

  std::vector<TVector3> ECAL_ext_points = GetExtrapoPoints("ECAL", rho, center, alpha0, CALO_trk_state, 
                                                          ECAL_barrel_delta_phi, ECAL_endcap_delta_phi);
  std::vector<TVector3> HCAL_ext_points = GetExtrapoPoints("HCAL", rho, center, alpha0, CALO_trk_state, 
                                                          HCAL_barrel_delta_phi, HCAL_endcap_delta_phi);  

  // Sort Extrapolated points
  std::vector<TrackState> t_ECAL_states;
  for(int ip=0; ip<ECAL_ext_points.size(); ip++){
    TrackState t_state = CALO_trk_state;
    t_state.location = Cyber::TrackState::AtOther;
    t_state.referencePoint = ECAL_ext_points[ip];
    // Note GetExtrapolatedPhi0 is not same as the definition of phi0 in TrackState
    t_state.phi0 = GetExtrapolatedPhi0(CALO_trk_state.Kappa, CALO_trk_state.phi0, center, ECAL_ext_points[ip]);
    t_ECAL_states.push_back(t_state);
  }
  // std::sort(t_ECAL_states.begin(), t_ECAL_states.end(), SortByPhi0);
  p_track->setTrackStates("Ecal", t_ECAL_states);

  std::vector<TrackState> t_HCAL_states;
  for(int ip=0; ip<HCAL_ext_points.size(); ip++){
    TrackState t_state = CALO_trk_state;
    t_state.location = Cyber::TrackState::AtOther;
    t_state.referencePoint = HCAL_ext_points[ip];
    // Note GetExtrapolatedPhi0 is not same as the definition of phi0 in TrackState
    t_state.phi0 = GetExtrapolatedPhi0(CALO_trk_state.Kappa, CALO_trk_state.phi0, center, HCAL_ext_points[ip]);
    t_HCAL_states.push_back(t_state);
  }
  // std::sort(t_HCAL_states.begin(), t_HCAL_states.end(), SortByPhi0);
  p_track->setTrackStates("Hcal", t_HCAL_states);

  return StatusCode::SUCCESS;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
float TrackExtrapolatingAlg::GetRho(const Cyber::TrackState & trk_state){
  float rho = Abs(1000. / (0.3*settings.map_floatPars["B_field"]*trk_state.Kappa));
  return rho;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
TVector2 TrackExtrapolatingAlg::GetCenterOfCircle(const Cyber::TrackState & trk_state, const float & rho){
  float phi;
  if(trk_state.Kappa>=0) phi = trk_state.phi0 - Pi()/2;
  else phi = trk_state.phi0 + Pi()/2;

  float xc = trk_state.referencePoint.X() + ((rho+trk_state.D0)*Cos(phi));
  float yc = trk_state.referencePoint.Y() + ((rho+trk_state.D0)*Sin(phi));
  TVector2 center(xc, yc);
  return center;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
float TrackExtrapolatingAlg::GetRefAlpha0(const Cyber::TrackState & trk_state, const TVector2 & center){
  float deltaX = trk_state.referencePoint.X() - center.X();
  float deltaY = trk_state.referencePoint.Y() - center.Y();
  float alpha0 = ATan2(deltaY, deltaX);
  return alpha0;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
StatusCode TrackExtrapolatingAlg::GetDeltaPhi(float rho, TVector2 center, float alpha0, 
const std::vector<float> & layer_radius, const std::vector<float> & layer_z, const Cyber::TrackState & CALO_trk_state,
vector<float> & barrel_delta_phi, vector<float> & endcap_delta_phi){
  // Barrel
  for(int il=0; il<layer_radius.size(); il++){
    float param0 = (Power(layer_radius[il], 2) - center.Mod2() - Power(rho, 2)) / (2 * rho * center.Mod());
    if(Abs(param0)>1) continue;
    float t_as = ASin(param0);

    float param1 = ATan2(center.X(), center.Y());

    float layer_delta_phi1 = t_as - alpha0 - param1;
    float layer_delta_phi2;
    if(t_as<0) layer_delta_phi2 = -Pi() - t_as - alpha0 - param1;
    else layer_delta_phi2 = Pi() - t_as - alpha0 - param1;

    while(layer_delta_phi1>Pi()) layer_delta_phi1 = layer_delta_phi1 - 2*Pi();
    while(layer_delta_phi1<-Pi()) layer_delta_phi1 = layer_delta_phi1 + 2*Pi();
    while(layer_delta_phi2>Pi()) layer_delta_phi2 = layer_delta_phi2 - 2*Pi();
    while(layer_delta_phi2<-Pi()) layer_delta_phi2 = layer_delta_phi2 + 2*Pi();

    if(CALO_trk_state.Kappa < 0){ 
      barrel_delta_phi.push_back(layer_delta_phi1); 
    }
    else if(CALO_trk_state.Kappa > 0){
      barrel_delta_phi.push_back(layer_delta_phi2); 
    }
    else{
      std::cout << "TrackExtrapolatingAlg: Error! Kappa=0!" << std::endl;
    }

  }

  // Endcap
  float z0 = CALO_trk_state.referencePoint.Z() + CALO_trk_state.Z0;
  for(int iz=0; iz<layer_z.size(); iz++){
    float layer_delta_phi;
    if (CALO_trk_state.tanLambda>0){
      layer_delta_phi = ( layer_z[iz] - z0 ) / ( rho * CALO_trk_state.tanLambda );
    }
    else{
      layer_delta_phi = ( -layer_z[iz] - z0 ) / ( rho * CALO_trk_state.tanLambda );
    }

    if(CALO_trk_state.Kappa < 0){
      endcap_delta_phi.push_back(layer_delta_phi);
    }
    else if(CALO_trk_state.Kappa > 0){
      endcap_delta_phi.push_back(-layer_delta_phi);
    }
    else{
      std::cout << "TrackExtrapolatingAlg: Error! Kappa=0!" << std::endl;
    }

  }

  return StatusCode::SUCCESS;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
std::vector<TVector3> TrackExtrapolatingAlg::GetExtrapoPoints(std::string calo_name, 
                                         float rho, TVector2 center, float alpha0, 
                                         const Cyber::TrackState & CALO_trk_state,
                                         const std::vector<float>& barrel_delta_phi, const std::vector<float>& endcap_delta_phi){
  std::vector<TVector3> barrel_ext_points;
  for(int ip=0; ip<barrel_delta_phi.size(); ip++){
    float x = center.X() + (rho*Cos(alpha0+barrel_delta_phi[ip]));
    float y = center.Y() + (rho*Sin(alpha0+barrel_delta_phi[ip]));
    float z;
    if(CALO_trk_state.Kappa > 0){
      z = CALO_trk_state.referencePoint.Z() + CALO_trk_state.Z0 - 
              (barrel_delta_phi[ip]*rho*CALO_trk_state.tanLambda);
    }else{
      z = CALO_trk_state.referencePoint.Z() + CALO_trk_state.Z0 + 
              (barrel_delta_phi[ip]*rho*CALO_trk_state.tanLambda);
    }
    TVector3 extp(x,y,z);
    
    if(calo_name=="ECAL"){
      if( Abs(z)<settings.map_floatPars["ECAL_half_length"] && 
          Sqrt(x*x+y*y) < settings.map_floatPars["ECAL_outermost_distance"] && 
          Sqrt(x*x+y*y) > settings.map_floatPars["ECAL_innermost_distance"]){
        barrel_ext_points.push_back(extp);
      }
    }
    else if(calo_name=="HCAL"){
      if( Abs(z)<settings.map_floatPars["HCAL_half_length"] && 
          Sqrt(x*x+y*y) < settings.map_floatPars["HCAL_outermost_distance"] && 
          Sqrt(x*x+y*y) > settings.map_floatPars["HCAL_innermost_distance"]){
        barrel_ext_points.push_back(extp);
      }
    }
    else continue;
  }

  std::vector<TVector3> endcap_ext_points;
  for(int ip=0; ip<endcap_delta_phi.size(); ip++){
    float x = center.X() + (rho*Cos(alpha0+endcap_delta_phi[ip]));
    float y = center.Y() + (rho*Sin(alpha0+endcap_delta_phi[ip]));
    float z;
    if(CALO_trk_state.Kappa > 0){
      z = CALO_trk_state.referencePoint.Z() + CALO_trk_state.Z0 - 
              (endcap_delta_phi[ip]*rho*CALO_trk_state.tanLambda);
    }else{
      z = CALO_trk_state.referencePoint.Z() + CALO_trk_state.Z0 + 
              (endcap_delta_phi[ip]*rho*CALO_trk_state.tanLambda);
    }
    TVector3 extp(x,y,z);

    if(calo_name=="ECAL"){
      if( Abs(z)>settings.map_floatPars["ECAL_endcap_zmin"] && 
          Abs(z)<settings.map_floatPars["ECAL_endcap_zmax"] && 
          Sqrt(x*x+y*y) < settings.map_floatPars["ECAL_outermost_distance"]){
        endcap_ext_points.push_back(extp);
      }
    }
    else if(calo_name=="HCAL"){
      if( Abs(z)>settings.map_floatPars["HCAL_endcap_zmin"] && 
          Abs(z)<settings.map_floatPars["HCAL_endcap_zmax"] && 
          Sqrt(x*x+y*y) < settings.map_floatPars["HCAL_outermost_distance"]){
        endcap_ext_points.push_back(extp);
      }
    }
    else{
      std::cout << "TrackExtrapolatingAlg: Error! Wrong calorimeter name!" << std::endl;
    }
  }

  std::vector<TVector3> ext_points;
  if(barrel_ext_points.size()>0 && endcap_ext_points.size()==0){
    ext_points = barrel_ext_points;
  }
  else if(barrel_ext_points.size()==0 && endcap_ext_points.size()>0){
    ext_points = endcap_ext_points;
  }
  else if(barrel_ext_points.size()>0 && endcap_ext_points.size()>0){
    ext_points = barrel_ext_points;

    float barrel_delta_z = TMath::Abs(barrel_ext_points[barrel_ext_points.size()-1].Z() - barrel_ext_points[barrel_ext_points.size()-2].Z());
    float endcap_barrel_delta_z = TMath::Abs(endcap_ext_points[0].Z() - barrel_ext_points[barrel_ext_points.size()-1].Z());
    if(endcap_barrel_delta_z < barrel_delta_z){
      ext_points.insert(ext_points.end(), endcap_ext_points.begin(), endcap_ext_points.end());
    }
  } 
  else{
    std::cout << "TrackExtrapolatingAlg: Error! No extrapolated points!" << std::endl;
  }

  return ext_points;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
bool TrackExtrapolatingAlg::IsReturn(float rho, TVector2 & center){
  float farest = rho + center.Mod();
  if (farest < settings.map_floatPars["ECAL_outermost_distance"]/Cos(Pi()/settings.map_intPars["Nmodule"])+100) {return true;}
  else{return false;}
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
float TrackExtrapolatingAlg::GetExtrapolatedPhi0(float Kappa, float ECAL_phi0, TVector2 center, TVector3 ext_point){
  // Note: phi0 of extrapolated points is (phi of velocity at extrapolated point) - (phi of velocity at ECAL front face)
  TVector2 ext_point_xy(ext_point.X(), ext_point.Y());
  TVector2 ext2center = center - ext_point_xy;
  float ext_phi0;
  if(Kappa>=0) ext_phi0 =  ext2center.Phi() + TMath::Pi()/2.;
  else ext_phi0 = ext2center.Phi() - TMath::Pi()/2.;

  float phi0 = ext_phi0 - ECAL_phi0;
  while(phi0 < -Pi())  phi0 = phi0 + 2*Pi();
  while(phi0 >  Pi())  phi0 = phi0 - 2*Pi();
  return phi0;
}


#endif
