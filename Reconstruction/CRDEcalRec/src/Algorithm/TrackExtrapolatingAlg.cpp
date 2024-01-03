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
    settings.map_floatPars["ECAL_innermost_distance"] = 1860;
  if(settings.map_floatPars.find("ECAL_outermost_distance")==settings.map_floatPars.end()) 
    settings.map_floatPars["ECAL_outermost_distance"] = 2140;
  if(settings.map_intPars.find("ECAL_Nlayers")==settings.map_intPars.end()) 
    settings.map_intPars["ECAL_Nlayers"] = 28;
  if(settings.map_floatPars.find("ECAL_layer_width")==settings.map_floatPars.end())
    settings.map_floatPars["ECAL_layer_width"] = 10;
  if(settings.map_floatPars.find("ECAL_half_length")==settings.map_floatPars.end())
    settings.map_floatPars["ECAL_half_length"] = 3300;
  // HCAL parameters
  if(settings.map_floatPars.find("HCAL_innermost_distance")==settings.map_floatPars.end()) 
    settings.map_floatPars["HCAL_innermost_distance"] = 2150;
  if(settings.map_floatPars.find("HCAL_outermost_distance")==settings.map_floatPars.end()) 
    settings.map_floatPars["HCAL_outermost_distance"] = 3230;
  if(settings.map_intPars.find("HCAL_Nlayers")==settings.map_intPars.end()) 
    settings.map_intPars["HCAL_Nlayers"] = 35;
  if(settings.map_floatPars.find("HCAL_layer_width")==settings.map_floatPars.end())
    settings.map_floatPars["HCAL_layer_width"] = 26.61;
  if(settings.map_floatPars.find("HCAL_sensitive_distance")==settings.map_floatPars.end())
    settings.map_floatPars["HCAL_sensitive_distance"] = 22.81;  // distance between sensitive material and front face of each layer
  if(settings.map_floatPars.find("HCAL_half_length")==settings.map_floatPars.end())
    settings.map_floatPars["HCAL_half_length"] = 4480;

  if(settings.map_intPars.find("Nmodule")==settings.map_intPars.end()) 
    settings.map_intPars["Nmodule"] = 8;
  if(settings.map_floatPars.find("B_field")==settings.map_floatPars.end())
    settings.map_floatPars["B_field"] = 3.0;

  if(settings.map_intPars.find("Input_track")==settings.map_intPars.end())
    settings.map_intPars["Input_track"] = 0;    // 0: reconstructed tracks.  1: MC particle track

  return StatusCode::SUCCESS;
};


StatusCode TrackExtrapolatingAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  std::cout<<"Initialize TrackExtrapolatingAlg"<<std::endl;

  return StatusCode::SUCCESS;
};


StatusCode TrackExtrapolatingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){
//std::cout<<"---oooOO0OOooo--- Excuting TrackExtrapolatingAlg ---oooOO0OOooo---"<<std::endl;

  std::vector<std::shared_ptr<PandoraPlus::Track>>* p_tracks = &(m_datacol.TrackCol);
//std::cout<<"  Track size: "<<p_tracks->size()<<std::endl;

//std::cout<<"  GetPlaneNormalVector() "<<std::endl;
  std::vector<TVector2> normal_vectors;
  GetPlaneNormalVector(normal_vectors);

//std::cout<<"  GetLayerPoints() "<<std::endl;
  std::vector<std::vector<TVector2>> ECAL_layer_points;    // 8 modules, 28 layer points in each modules
  std::vector<std::vector<TVector2>> HCAL_layer_points;    // 8 modules, 35 layer points in each modules
  GetLayerPoints(normal_vectors, ECAL_layer_points, HCAL_layer_points);

  for(int itrk=0; itrk<p_tracks->size(); itrk++){
    // Only tracks that reach ECAL should be processed.
    if(!IsReachECAL( p_tracks->at(itrk).get() )) continue;

//std::cout<<"  GetTrackStateAtCalo() "<<std::endl;
    // get track state at calorimeter
    PandoraPlus::TrackState CALO_trk_state;
    GetTrackStateAtCalo(p_tracks->at(itrk).get(), CALO_trk_state);
    
//std::cout<<"  ExtrapolateByLayer() "<<std::endl;
    ExtrapolateByLayer(normal_vectors, ECAL_layer_points, HCAL_layer_points, CALO_trk_state, p_tracks->at(itrk).get());
  } // end loop tracks


  p_tracks = nullptr;
  return StatusCode::SUCCESS;
}; // RunAlgorithm end


StatusCode TrackExtrapolatingAlg::ClearAlgorithm(){
  std::cout<<"End run TrackExtrapolatingAlg. Clean it."<<std::endl;

  return StatusCode::SUCCESS;
};


// StatusCode TrackExtrapolatingAlg::SelfAlg1(){
//   std::cout<<"  Processing SelfAlg1: print Par1 = "<<settings.map_floatPars["Par1"]<<std::endl;

//   return StatusCode::SUCCESS;
// };



// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
StatusCode TrackExtrapolatingAlg::GetPlaneNormalVector(std::vector<TVector2> & normal_vectors){
  for(int im=0; im<settings.map_intPars["Nmodule"]; im++){    
    TVector2 t_vec(0, 1);
    t_vec = t_vec.Rotate(1.*im/settings.map_intPars["Nmodule"]*2*Pi());
    normal_vectors.push_back(t_vec);
  }

  return StatusCode::SUCCESS;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
StatusCode TrackExtrapolatingAlg::GetLayerPoints(const std::vector<TVector2> & normal_vectors, 
                          std::vector<std::vector<TVector2>> & ECAL_layer_points,
                          std::vector<std::vector<TVector2>> & HCAL_layer_points){
  // ECAL
  for(int im=0; im< normal_vectors.size(); im++){
    std::vector<TVector2> t_points;
    for(int il=0; il<settings.map_intPars["ECAL_Nlayers"]; il++){
      TVector2 t_p;
      float dist = settings.map_floatPars["ECAL_innermost_distance"] + 
                   settings.map_floatPars["ECAL_layer_width"]*il + settings.map_floatPars["ECAL_layer_width"]/2;
      float xx = (normal_vectors[im].X()/normal_vectors[im].Mod()) * dist;
      float yy = (normal_vectors[im].Y()/normal_vectors[im].Mod()) * dist;
      t_p.SetX(xx); t_p.SetY(yy);
      t_points.push_back(t_p);
    }
    ECAL_layer_points.push_back(t_points);
  }

  // HCAL
  for(int im=0; im< normal_vectors.size(); im++){
    std::vector<TVector2> t_points;
    for(int il=0; il<settings.map_intPars["HCAL_Nlayers"]; il++){
      TVector2 t_p;
      float dist = settings.map_floatPars["HCAL_innermost_distance"] + 
                   settings.map_floatPars["HCAL_layer_width"]*il + settings.map_floatPars["HCAL_sensitive_distance"];
      float xx = (normal_vectors[im].X()/normal_vectors[im].Mod()) * dist;
      float yy = (normal_vectors[im].Y()/normal_vectors[im].Mod()) * dist;
      t_p.SetX(xx); t_p.SetY(yy);
      t_points.push_back(t_p);
    }
    HCAL_layer_points.push_back(t_points);
  }

  return StatusCode::SUCCESS;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
bool TrackExtrapolatingAlg::IsReachECAL(PandoraPlus::Track * track){
  if(settings.map_intPars["Input_track"] == 0){
    // The track is reconstructed in ECAL. If the track reach ECAL, it should have track state at calorimeter
    std::vector<TrackState> input_trackstates = track->getTrackStates("Input");
    int count=0;
    TVector3 t_vec;
    for(int i=0; i<input_trackstates.size(); i++){
      if(input_trackstates[i].location==PandoraPlus::TrackState::AtCalorimeter){
        count++;
        t_vec = input_trackstates[i].referencePoint;
        break;
      }
    }
    if(count==0){
  std::cout<<"      the track has no track state at calorimeter"<<std::endl;
      return false;
    } 
    if( Abs(Abs(t_vec.Z())-settings.map_floatPars["ECAL_half_length"]) < 1e-6 ){
  std::cout<<"      the track escape from endcap"<<std::endl;
      return false;
    } 
    

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
      if(input_trackstates[i].location==PandoraPlus::TrackState::AtIP)
        IP_trk_state = input_trackstates[i];
      break;
    }

    TVector3 ref_point = IP_trk_state.referencePoint;
    double rho = GetRho(IP_trk_state);
    double r_max = TMath::Sqrt(ref_point.X()*ref_point.X() + ref_point.Y()*ref_point.Y()) + rho*2;

    if(r_max<settings.map_floatPars["ECAL_innermost_distance"]){ return false; }

    return true;
  }
  else{
    std::cout << "Error, wrong source of input tracks for TrackExtrapolatingAlg!" << std:: endl;
    return false;
  }
  
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
StatusCode TrackExtrapolatingAlg::GetTrackStateAtCalo(PandoraPlus::Track * track, 
                                                PandoraPlus::TrackState & trk_state_at_calo){
  std::vector<TrackState> input_trackstates = track->getTrackStates("Input");

  if(settings.map_intPars["Input_track"] == 0){
    for(int its=0; its<input_trackstates.size(); its++){
      if(input_trackstates[its].location==PandoraPlus::TrackState::AtCalorimeter){
        trk_state_at_calo=input_trackstates[its];
        break;
      }
    }
  }
  else if((settings.map_intPars["Input_track"] == 1)){
    for(int its=0; its<input_trackstates.size(); its++){
      if(input_trackstates[its].location==PandoraPlus::TrackState::AtIP){
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
StatusCode TrackExtrapolatingAlg::ExtrapolateByLayer(const std::vector<TVector2> & normal_vectors,
                              const std::vector<std::vector<TVector2>> & ECAL_layer_points, 
                              const std::vector<std::vector<TVector2>> & HCAL_layer_points, 
                              const PandoraPlus::TrackState & CALO_trk_state, 
                              PandoraPlus::Track* p_track){
  float rho = GetRho(CALO_trk_state);
  TVector2 center = GetCenterOfCircle(CALO_trk_state, rho);
  float alpha0 = GetRefAlpha0(CALO_trk_state, center); 
  // bool is_rtbk = IsReturn(rho, center);
  // Evaluate delta_phi
  std::vector<std::vector<float>> ECAL_delta_phi = GetDeltaPhi(rho, center, alpha0, normal_vectors, ECAL_layer_points, CALO_trk_state);
  std::vector<std::vector<float>> HCAL_delta_phi = GetDeltaPhi(rho, center, alpha0, normal_vectors, HCAL_layer_points, CALO_trk_state);
  
  // extrpolated track points. Will be stored as reference point in TrackState
  std::vector<TVector3> ECAL_ext_points = GetExtrapoPoints("ECAL", rho, center, alpha0, CALO_trk_state, ECAL_delta_phi);
  std::vector<TVector3> HCAL_ext_points = GetExtrapoPoints("HCAL", rho, center, alpha0, CALO_trk_state, HCAL_delta_phi);
  
  // Sort Extrapolated points
  std::vector<TrackState> t_ECAL_states;
  for(int ip=0; ip<ECAL_ext_points.size(); ip++){
    TrackState t_state = CALO_trk_state;
    t_state.location = PandoraPlus::TrackState::AtOther;
    t_state.referencePoint = ECAL_ext_points[ip];
    // Note GetExtrapolatedPhi0 is not same as the definition of phi0 in TrackState
    t_state.phi0 = GetExtrapolatedPhi0(CALO_trk_state.Kappa, CALO_trk_state.phi0, center, ECAL_ext_points[ip]);
    t_ECAL_states.push_back(t_state);
  }
  std::sort(t_ECAL_states.begin(), t_ECAL_states.end(), SortByPhi0);
  p_track->setTrackStates("Ecal", t_ECAL_states);

  std::vector<TrackState> t_HCAL_states;
  for(int ip=0; ip<HCAL_ext_points.size(); ip++){
    TrackState t_state = CALO_trk_state;
    t_state.location = PandoraPlus::TrackState::AtOther;
    t_state.referencePoint = HCAL_ext_points[ip];
    // Note GetExtrapolatedPhi0 is not same as the definition of phi0 in TrackState
    t_state.phi0 = GetExtrapolatedPhi0(CALO_trk_state.Kappa, CALO_trk_state.phi0, center, HCAL_ext_points[ip]);
    t_HCAL_states.push_back(t_state);
  }
  std::sort(t_HCAL_states.begin(), t_HCAL_states.end(), SortByPhi0);
  p_track->setTrackStates("Hcal", t_HCAL_states);
  
  return StatusCode::SUCCESS;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
float TrackExtrapolatingAlg::GetRho(const PandoraPlus::TrackState & trk_state){
  float rho = Abs(1000. / (0.3*settings.map_floatPars["B_field"]*trk_state.Kappa));
  return rho;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
TVector2 TrackExtrapolatingAlg::GetCenterOfCircle(const PandoraPlus::TrackState & trk_state, const float & rho){
  float phi;
  if(trk_state.Kappa>=0) phi = trk_state.phi0 - Pi()/2;
  else phi = trk_state.phi0 + Pi()/2;

  float xc = trk_state.referencePoint.X() + ((rho+trk_state.D0)*Cos(phi));
  float yc = trk_state.referencePoint.Y() + ((rho+trk_state.D0)*Sin(phi));
  TVector2 center(xc, yc);
  return center;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
float TrackExtrapolatingAlg::GetRefAlpha0(const PandoraPlus::TrackState & trk_state, const TVector2 & center){
  float deltaX = trk_state.referencePoint.X() - center.X();
  float deltaY = trk_state.referencePoint.Y() - center.Y();
  float alpha0 = ATan2(deltaY, deltaX);
  return alpha0;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
std::vector<std::vector<float>> TrackExtrapolatingAlg::GetDeltaPhi(float rho, TVector2 center, float alpha0,
                              const std::vector<TVector2> & normal_vectors,
                              const std::vector<std::vector<TVector2>> & layer_points, 
                              const PandoraPlus::TrackState & CALO_trk_state){
  std::vector<std::vector<float>> delta_phi;
  for(int im=0; im<normal_vectors.size(); im++){  // for each module
    std::vector<float> t_delta_phi;
    float beta = ATan2(normal_vectors[im].X(), normal_vectors[im].Y());
    float denominator = rho * Sqrt( (normal_vectors[im].X()*normal_vectors[im].X()) +
                                    (normal_vectors[im].Y()*normal_vectors[im].Y()) );
    for(int il=0; il<layer_points[im].size(); il++){  // for each layer
      float numerator = (layer_points[im][il].X()*normal_vectors[im].X()) -
                        (center.X()*normal_vectors[im].X()) + 
                        (layer_points[im][il].Y()*normal_vectors[im].Y()) -
                        (center.Y()*normal_vectors[im].Y()) ;
      if(Abs(numerator/denominator)>1) continue;
      float t_as = ASin(numerator/denominator);
      float t_ab = ASin(Sin(alpha0+beta));

      float t_dphi1 = t_as - t_ab;
      // float t_dphi2 = Pi()-t_as - t_ab;
      if(CALO_trk_state.Kappa < 0){ 
        t_delta_phi.push_back(t_dphi1); 
        // if(is_rtbk){
        //   t_delta_phi.push_back(t_dphi2);
        // }
      }
      else{
        t_delta_phi.push_back(-t_dphi1); 
        // if(is_rtbk){
        //   t_delta_phi.push_back(-t_dphi2);
        // }
      }
    }
    delta_phi.push_back(t_delta_phi);
  }

  return delta_phi;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
std::vector<TVector3> TrackExtrapolatingAlg::GetExtrapoPoints(std::string calo_name, 
                                         float rho, TVector2 center, float alpha0, 
                                         const PandoraPlus::TrackState & CALO_trk_state,
                                         const std::vector<std::vector<float>>& delta_phi){
  std::vector<TVector3> ext_points;
  for(int im=0; im<delta_phi.size(); im++){
    for(int ip=0; ip<delta_phi[im].size(); ip++){
      float x = center.X() + (rho*Cos(alpha0+delta_phi[im][ip]));
      float y = center.Y() + (rho*Sin(alpha0+delta_phi[im][ip]));
      float z;
      if(CALO_trk_state.Kappa > 0){
        z = CALO_trk_state.referencePoint.Z() + CALO_trk_state.Z0 - 
                (delta_phi[im][ip]*rho*CALO_trk_state.tanLambda);
      }else{
        z = CALO_trk_state.referencePoint.Z() + CALO_trk_state.Z0 + 
                (delta_phi[im][ip]*rho*CALO_trk_state.tanLambda);
      }

      TVector3 t_vec3(x,y,z);
      float rotate_angle = (6-im)*Pi()/4;
      t_vec3.RotateZ(rotate_angle);

      if(calo_name=="ECAL"){
        if(Abs(z)>settings.map_floatPars["ECAL_half_length"]) continue;
        if(Sqrt(x*x+y*y) > settings.map_floatPars["ECAL_outermost_distance"]/Cos(Pi()/settings.map_intPars["Nmodule"])+100) continue;
        if(Sqrt(x*x+y*y) < settings.map_floatPars["ECAL_innermost_distance"]) continue;
        if(settings.map_floatPars["ECAL_outermost_distance"]*Sqrt(2)-t_vec3.X() < t_vec3.Y()) continue;
        if(t_vec3.X()-settings.map_floatPars["ECAL_innermost_distance"]*Sqrt(2) > t_vec3.Y()) continue;
      }
      else if(calo_name=="HCAL"){
        if(Abs(z)>settings.map_floatPars["HCAL_half_length"]) continue;
        if(Sqrt(x*x+y*y) > settings.map_floatPars["HCAL_outermost_distance"]) continue;
        if(Sqrt(x*x+y*y) < settings.map_floatPars["HCAL_innermost_distance"]) continue;
        if(t_vec3.X()-settings.map_floatPars["HCAL_innermost_distance"]*Sqrt(2) > t_vec3.Y()) continue;
      }
      else continue;
      
      TVector3 extp(x,y,z);
      ext_points.push_back(extp);
    }
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
