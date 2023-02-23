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
  // if(settings.map_floatPars.find("Par1")==settings.map_floatPars.end()) settings.map_floatPars["Par1"] = 0;
  if(settings.map_floatPars.find("innermost_distance")==settings.map_floatPars.end()) 
    settings.map_floatPars["innermost_distance"] = 1860;
  if(settings.map_floatPars.find("outermost_distance")==settings.map_floatPars.end()) 
    settings.map_floatPars["outermost_distance"] = 2140;
  if(settings.map_intPars.find("Nlayers")==settings.map_intPars.end()) 
    settings.map_intPars["Nlayers"] = 28;
  if(settings.map_floatPars.find("layer_width")==settings.map_floatPars.end())
    settings.map_floatPars["layer_width"] = 10;
  if(settings.map_floatPars.find("ECAL_half_length")==settings.map_floatPars.end())
    settings.map_floatPars["ECAL_half_length"] = 2450;
  if(settings.map_intPars.find("Nmodule")==settings.map_intPars.end()) 
    settings.map_intPars["Nmodule"] = 8;
  if(settings.map_floatPars.find("B_field")==settings.map_floatPars.end())
    settings.map_floatPars["B_field"] = 3.0;

  return StatusCode::SUCCESS;
};


StatusCode TrackExtrapolatingAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  std::cout<<"Initialize TrackExtrapolatingAlg"<<std::endl;

  return StatusCode::SUCCESS;
};


StatusCode TrackExtrapolatingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){
std::cout<<"---oooOO0OOooo---Excute TrackExtrapolatingAlg---oooOO0OOooo---"<<std::endl;

  std::vector<PandoraPlus::Track*>* p_tracks = &(m_datacol.TrackCol);
std::cout<<"  Track size: "<<p_tracks->size()<<std::endl;

  std::vector<TVector2> normal_vectors;
  GetPlaneNormalVector(normal_vectors);

  std::vector<std::vector<TVector2>> layer_points;    // 8 modules, 28 layer points in each modules
  GetLayerPoints(normal_vectors, layer_points);

  for(int itrk=0; itrk<p_tracks->size(); itrk++){
// std::cout<<"    --------------------Processing track "<<itrk<<"--------------------"<<std::endl;
    // Only tracks that reach ECAL should be processed.
    if(!IsReachECAL(p_tracks->at(itrk))) continue;

    // get track state at calorimeter
    PandoraPlus::TrackState ECAL_trk_state;
    GetECALTrackState(p_tracks->at(itrk), ECAL_trk_state);
    
// std::cout<<"      track location: "<<ECAL_trk_state.location<<std::endl;
// std::cout<<"      track reference point: ("<<ECAL_trk_state.referencePoint.X()<<", "
//                                          <<ECAL_trk_state.referencePoint.Y()<<", "
//                                          <<ECAL_trk_state.referencePoint.Z()<<")"<<std::endl;

    ExtrapolateByLayer(normal_vectors, layer_points, ECAL_trk_state, p_tracks->at(itrk));
    GetTrackPoints(ECAL_trk_state, p_tracks->at(itrk));
  } // end loop tracks

  // std::cout<<"  Run SelfAlg1"<<std::endl;
  // SelfAlg1();

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
                          std::vector<std::vector<TVector2>> & layer_points){
  for(int im=0; im< normal_vectors.size(); im++){
    std::vector<TVector2> t_points;
    for(int il=0; il<settings.map_intPars["Nlayers"]; il++){
      TVector2 t_p;
      float dist = settings.map_floatPars["layer_width"]*il + settings.map_floatPars["innermost_distance"] + 
                   settings.map_floatPars["layer_width"]/2;
      float xx = (normal_vectors[im].X()/normal_vectors[im].Mod()) * dist;
      float yy = (normal_vectors[im].Y()/normal_vectors[im].Mod()) * dist;
      t_p.SetX(xx); t_p.SetY(yy);
      t_points.push_back(t_p);
    }
    layer_points.push_back(t_points);
  }

  return StatusCode::SUCCESS;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
bool TrackExtrapolatingAlg::IsReachECAL(PandoraPlus::Track * track){
  int count=0;
  TVector3 t_vec;
  for(int i=0; i<track->trackStates_size(); i++){
    if(track->getTrackStates()[i].location==PandoraPlus::TrackState::AtCalorimeter){
      count++;
      t_vec = track->getTrackStates()[i].referencePoint;
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


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
StatusCode TrackExtrapolatingAlg::GetECALTrackState(PandoraPlus::Track * track, 
                             PandoraPlus::TrackState & ECAL_trk_state){
  for(int its=0; its<track->trackStates_size();its++){
      if(track->getTrackStates(its).location==PandoraPlus::TrackState::AtCalorimeter){
        ECAL_trk_state=track->getTrackStates(its);
        break;
      }
  }
  return StatusCode::SUCCESS;
}


// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
StatusCode TrackExtrapolatingAlg::ExtrapolateByLayer(const std::vector<TVector2> & normal_vectors,
                              const std::vector<std::vector<TVector2>> & layer_points, 
                              const PandoraPlus::TrackState & ECAL_trk_state, 
                              PandoraPlus::Track* p_track){
    float rho = GetRho(ECAL_trk_state);
    TVector2 center = GetCenterOfCircle(ECAL_trk_state, rho);
    float alpha0 = GetRefAlpha0(ECAL_trk_state, center); 
    bool is_rtbk = IsReturn(rho, center);
    // Evaluate delta_phi
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
        float t_dphi2 = Pi()-t_as - t_ab;
        if(ECAL_trk_state.Kappa < 0){ 
          t_delta_phi.push_back(t_dphi1); 
          if(is_rtbk){
            t_delta_phi.push_back(t_dphi2);
          }
        }
        else{
          t_delta_phi.push_back(-t_dphi1); 
          if(is_rtbk){
            t_delta_phi.push_back(-t_dphi2);
          }
        }
        
      }
      delta_phi.push_back(t_delta_phi);
    }

    // extrpolated track points. Will be stored as reference point in TrackState
    std::vector<TVector3> ext_points;
    for(int im=0; im<delta_phi.size(); im++){
      for(int ip=0; ip<delta_phi[im].size(); ip++){
        float x = center.X() + (rho*Cos(alpha0+delta_phi[im][ip]));
        float y = center.Y() + (rho*Sin(alpha0+delta_phi[im][ip]));
        float z;
        if(ECAL_trk_state.Kappa > 0){
          z = ECAL_trk_state.referencePoint.Z() + ECAL_trk_state.Z0 - 
                  (delta_phi[im][ip]*rho*ECAL_trk_state.tanLambda);
        }else{
          z = ECAL_trk_state.referencePoint.Z() + ECAL_trk_state.Z0 + 
                  (delta_phi[im][ip]*rho*ECAL_trk_state.tanLambda);
        }

        TVector3 t_vec3(x,y,z);
        float rotate_angle = (6-im)*Pi()/4;
        t_vec3.RotateZ(rotate_angle);

        if(Sqrt(x*x+y*y) > settings.map_floatPars["outermost_distance"]/Cos(Pi()/settings.map_intPars["Nmodule"])+100) continue;
        if(Sqrt(x*x+y*y) < settings.map_floatPars["innermost_distance"]) continue;
        if(settings.map_floatPars["outermost_distance"]*Sqrt(2)-t_vec3.X() < t_vec3.Y()) continue;
        if(t_vec3.X()-settings.map_floatPars["innermost_distance"]*Sqrt(2) > t_vec3.Y()) continue;
        
        TVector3 extp(x,y,z);
        ext_points.push_back(extp);
      }
    }
    
    // Add track state to the track
    for(int ip=0; ip<ext_points.size(); ip++){
      TrackState t_state = ECAL_trk_state;
      t_state.location = PandoraPlus::TrackState::AtOther;
      t_state.referencePoint = ext_points[ip];
      p_track->AddTrackState(t_state);
    }
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
bool TrackExtrapolatingAlg::IsReturn(float rho, TVector2 & center){
  float farest = rho + center.Mod();
  if (farest < settings.map_floatPars["outermost_distance"]/Cos(Pi()/settings.map_intPars["Nmodule"])+100) {return true;}
  else{return false;}
}

// ...oooOO0OOooo......oooOO0OOooo......oooOO0OOooo...
StatusCode TrackExtrapolatingAlg::GetTrackPoints(const PandoraPlus::TrackState & ECAL_trk_state, 
                                                 PandoraPlus::Track* p_track){
  float rho = GetRho(ECAL_trk_state);
  TVector2 center = GetCenterOfCircle(ECAL_trk_state, rho);
  float alpha0 = GetRefAlpha0(ECAL_trk_state, center); 

  // the angle from center point(x_C, y_c) of the helix circle to IP(0, 0)
  float IP_angle = ATan2( -center.Y(), -center.X() );

  vector<float> delta_phi;

  // negative charged particle
  // IP_angle should be < alpha0
  if(ECAL_trk_state.Kappa < 0){
    while ( IP_angle > alpha0 ){
      IP_angle = IP_angle - 2*Pi();
    }
  }
  // positive charged particle
  // IP_angle should be > alpha0
  else{
    while (IP_angle < alpha0){
      IP_angle = IP_angle + 2*Pi();
    }
  }

  for(int i = 0; i <= 100; i++){
    float t_phi = (IP_angle - alpha0) / 100. * i;
    delta_phi.push_back(t_phi);
  }

  std::vector<TVector3> ext_points;
  for(int ip=0; ip<delta_phi.size(); ip++){
    float x = center.X() + (rho*Cos(alpha0+delta_phi[ip]));
    float y = center.Y() + (rho*Sin(alpha0+delta_phi[ip]));
    float z;
    if(ECAL_trk_state.Kappa > 0){
      z = ECAL_trk_state.referencePoint.Z() + ECAL_trk_state.Z0 - 
              (delta_phi[ip]*rho*ECAL_trk_state.tanLambda);
    }else{
      z = ECAL_trk_state.referencePoint.Z() + ECAL_trk_state.Z0 + 
              (delta_phi[ip]*rho*ECAL_trk_state.tanLambda);
    }
    TVector3 t_vec3(x,y,z);
    ext_points.push_back(t_vec3);
  }

  for(int ip=0; ip<ext_points.size(); ip++){
    TrackState t_state = ECAL_trk_state;
    t_state.location = PandoraPlus::TrackState::AtOther;
    t_state.referencePoint = ext_points[ip];
    p_track->AddTrackState(t_state);
  }

  return StatusCode::SUCCESS;
}

#endif
