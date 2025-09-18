#ifndef TRACKMATCHING_C
#define TRACKMATCHING_C

#include "Algorithm/TrackMatchingAlg.h"


StatusCode TrackMatchingAlg::ReadSettings(Cyber::Settings& m_settings){
  settings = m_settings;
  // ECAL geometry settings
  // Note: Bar half length is also geometry parameter, but obtained from the function GetBarHalfLength()
  if(settings.map_floatPars.find("localmax_area")==settings.map_floatPars.end())
    settings.map_floatPars["localmax_area"] = 10; // unit: mm
  if(settings.map_floatPars.find("ConeNearByDistance")==settings.map_floatPars.end())
    settings.map_floatPars["ConeNearByDistance"] = 100;
  if(settings.map_floatPars.find("ConeMatchingCut_pT")==settings.map_floatPars.end())
    settings.map_floatPars["ConeMatchingCut_pT"] = 2.0; // GeV. If pT of a track < ConeMatchingCut_pT, use Cone matching
  if(settings.map_intPars.find("Max_Seed_Point")==settings.map_intPars.end())
    settings.map_intPars["Max_Seed_Point"] = 4;
  if(settings.map_floatPars.find("ConeSeedDistance")==settings.map_floatPars.end())
    settings.map_floatPars["ConeSeedDistance"] = 20;
  if(settings.map_floatPars.find("th_ConeTheta")==settings.map_floatPars.end())    
    settings.map_floatPars["th_ConeTheta"] = TMath::Pi()/4.;
  if(settings.map_floatPars.find("th_ConeR")==settings.map_floatPars.end())        
    settings.map_floatPars["th_ConeR"] = 50;
  if(settings.map_stringPars.find("ReadinLocalMaxName")==settings.map_stringPars.end())
    settings.map_stringPars["ReadinLocalMaxName"] = "AllLocalMax";
  if(settings.map_stringPars.find("OutputLongiClusName")==settings.map_stringPars.end()) 
    settings.map_stringPars["OutputLongiClusName"] = "TrackAxis"; 

  return StatusCode::SUCCESS;
}


StatusCode TrackMatchingAlg::Initialize( CyberDataCol& m_datacol ){
  m_TrackCol.clear();
  p_HalfClusterV = nullptr;
  p_HalfClusterU = nullptr;
  // m_trackAxisVCol.clear();
  // m_trackAxisUCol.clear();

  for(int itrk=0; itrk<m_datacol.TrackCol.size(); itrk++ ) m_TrackCol.push_back(m_datacol.TrackCol[itrk].get());
  p_HalfClusterU = &(m_datacol.map_HalfCluster["HalfClusterColU"]);
  p_HalfClusterV = &(m_datacol.map_HalfCluster["HalfClusterColV"]);

  return StatusCode::SUCCESS;
}


StatusCode TrackMatchingAlg::RunAlgorithm( CyberDataCol& m_datacol ){
  //std::cout << "---oooOO0OOooo---Excuting TrackMatchingAlg---oooOO0OOooo---"<<std::endl;
  // Associate tracks to HalfClusters.
  // This association is a many-to-many relationship: 
  //    One HalCluster may have multiple tracks; 
  //    One track may pass through multiple HalfClusters.
  //cout<<"track size: "<<m_TrackCol.size()<<", HFClusterU size "<<p_HalfClusterU->size()<<", HFClusterV size "<<p_HalfClusterV->size()<<endl;


  for(int itrk=0; itrk<m_TrackCol.size(); itrk++){  // loop tracks
    if(m_TrackCol[itrk]->getTrackStates("Ecal").size()==0) continue;

    // Get extrapolated points of the track. These points are sorted by the track
    std::vector<TVector3> extrapo_points;
    GetExtrpoECALPoints(m_TrackCol[itrk], extrapo_points);
    if(extrapo_points.size()==0) continue;

    double pT = TMath::Abs(1. / m_TrackCol[itrk]->getTrackStates("Ecal")[0].Kappa);
    if (pT >= settings.map_floatPars["ConeMatchingCut_pT"]){
      for(int ihc=0; ihc<p_HalfClusterV->size(); ihc++){  // loop HalfClusterV
        // Get local max of the HalfCluster
        std::vector<const Cyber::Calo1DCluster*> localMaxColV = p_HalfClusterV->at(ihc).get()->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

        // Track axis candidate.
        std::shared_ptr<Cyber::CaloHalfCluster> t_track_axis = std::make_shared<Cyber::CaloHalfCluster>();
        CreateTrackAxis(extrapo_points, localMaxColV, t_track_axis.get());

        // If the track does not match the Halfcluster, the track axis candidate will have no 1DCluster
        if(t_track_axis->getCluster().size()==0)
          continue;

        t_track_axis->addAssociatedTrack(m_TrackCol[itrk]);
        t_track_axis->setType(10000); //Track-type axis. 
        m_TrackCol[itrk]->addAssociatedHalfClusterV( p_HalfClusterV->at(ihc).get() );
        m_datacol.map_HalfCluster["bkHalfCluster"].push_back(t_track_axis);
        p_HalfClusterV->at(ihc).get()->addHalfCluster(settings.map_stringPars["OutputLongiClusName"], t_track_axis.get());
      }  // end loop HalfClusterV

      for(int ihc=0; ihc<p_HalfClusterU->size(); ihc++){  // loop HalfClusterU
        // Get local max of the HalfCluster
        std::vector<const Cyber::Calo1DCluster*> localMaxColU = p_HalfClusterU->at(ihc).get()->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

        // Track axis candidate.
        std::shared_ptr<Cyber::CaloHalfCluster> t_track_axis = std::make_shared<Cyber::CaloHalfCluster>();
        CreateTrackAxis(extrapo_points, localMaxColU, t_track_axis.get());

        // If the track do not match the Halfcluster, the track axis candidate will have no 1DCluster
        if(t_track_axis->getCluster().size()==0)
          continue;
        
        t_track_axis->addAssociatedTrack(m_TrackCol[itrk]);
        t_track_axis->setType(10000); //Track-type axis. 
        m_TrackCol[itrk]->addAssociatedHalfClusterU( p_HalfClusterU->at(ihc).get() );
        m_datacol.map_HalfCluster["bkHalfCluster"].push_back(t_track_axis);
        p_HalfClusterU->at(ihc).get()->addHalfCluster(settings.map_stringPars["OutputLongiClusName"], t_track_axis.get());
      }  // end loop HalfClusterU
    }
    else{  // pT < settings.map_floatPars["ConeMatchingCut_pT"]
      // Get local max and HalfCluster near the extrapolated points
      std::vector<Cyber::CaloHalfCluster*> t_nearbyHalfClustersV;  t_nearbyHalfClustersV.clear();
      std::vector<Cyber::CaloHalfCluster*> t_nearbyHalfClustersU;  t_nearbyHalfClustersU.clear();
      std::vector<const Cyber::Calo1DCluster*> t_nearbyLocalMaxV;     t_nearbyLocalMaxV.clear();
      std::vector<const Cyber::Calo1DCluster*> t_nearbyLocalMaxU;     t_nearbyLocalMaxU.clear();
      GetNearby(p_HalfClusterV, extrapo_points, t_nearbyHalfClustersV, t_nearbyLocalMaxV);
      GetNearby(p_HalfClusterU, extrapo_points, t_nearbyHalfClustersU, t_nearbyLocalMaxU);

      // V plane
      std::vector<const Cyber::Calo1DCluster*> t_cone_axisV; t_cone_axisV.clear();
      LongiConeLinking(extrapo_points, t_nearbyLocalMaxV, t_cone_axisV);
      CreatConeAxis(m_datacol, m_TrackCol[itrk], t_nearbyHalfClustersV, t_cone_axisV);

      // U plane
      std::vector<const Cyber::Calo1DCluster*> t_cone_axisU; t_cone_axisU.clear();
      LongiConeLinking(extrapo_points, t_nearbyLocalMaxU, t_cone_axisU);
      CreatConeAxis(m_datacol, m_TrackCol[itrk], t_nearbyHalfClustersU, t_cone_axisU);

    }

  }  

  //Loop track to check the associated cluster: merge clusters if they are associated to the same track.
  std::vector<Cyber::CaloHalfCluster*> tmp_deleteClus; tmp_deleteClus.clear();
  for(auto &itrk : m_TrackCol){
    std::vector<Cyber::CaloHalfCluster*> m_matchedUCol = itrk->getAssociatedHalfClustersU();
    std::vector<Cyber::CaloHalfCluster*> m_matchedVCol = itrk->getAssociatedHalfClustersV();

    for(int imc=0; imc<m_matchedVCol.size(); imc++){
      int N_trk_axis = m_matchedVCol[imc]->getHalfClusterMap()[settings.map_stringPars["OutputLongiClusName"]].size() ;
    }

    if( m_matchedUCol.size()>1 ){
      for(int i=1; i<m_matchedUCol.size(); i++){ 
        m_matchedUCol[0]->mergeHalfCluster( m_matchedUCol[i] );
        tmp_deleteClus.push_back(m_matchedUCol[i]);
      }
    }
    if( m_matchedVCol.size()>1 ){
      for(int i=1; i<m_matchedVCol.size(); i++){
        m_matchedVCol[0]->mergeHalfCluster( m_matchedVCol[i] );
        tmp_deleteClus.push_back(m_matchedVCol[i]);
      }
    }

  }

  //Check vector: clean the merged clusters
  for(int ihc=0; ihc<p_HalfClusterU->size(); ihc++){
    if( find(tmp_deleteClus.begin(), tmp_deleteClus.end(), p_HalfClusterU->at(ihc).get())!=tmp_deleteClus.end() ){
      p_HalfClusterU->erase(p_HalfClusterU->begin()+ihc);
      ihc--;
    }
  }

  for(int ihc=0; ihc<p_HalfClusterV->size(); ihc++){
    if( find(tmp_deleteClus.begin(), tmp_deleteClus.end(), p_HalfClusterV->at(ihc).get())!=tmp_deleteClus.end() ){
      p_HalfClusterV->erase(p_HalfClusterV->begin()+ihc);
      ihc--;
    }
  }

//for(int ic=0; ic<p_HalfClusterU->size(); ic++){
//printf("HalfClusterU #%d: energy %.4f, track axis size %d \n", ic, p_HalfClusterU->at(ic)->getEnergy(), p_HalfClusterU->at(ic)->getHalfClusterCol(settings.map_stringPars["OutputLongiClusName"]).size() );
//}
//for(int ic=0; ic<p_HalfClusterV->size(); ic++){
//printf("HalfClusterV #%d: energy %.4f, track axis size %d \n", ic, p_HalfClusterV->at(ic)->getEnergy(), p_HalfClusterV->at(ic)->getHalfClusterCol(settings.map_stringPars["OutputLongiClusName"]).size() );
//}

  return StatusCode::SUCCESS;
}


StatusCode TrackMatchingAlg::ClearAlgorithm(){
  m_TrackCol.clear();
  p_HalfClusterV = nullptr;
  p_HalfClusterU = nullptr;
  // m_trackAxisVCol.clear();
  // m_trackAxisUCol.clear();

  return StatusCode::SUCCESS;
}


StatusCode TrackMatchingAlg::GetExtrpoECALPoints(const Cyber::Track* track, std::vector<TVector3>& extrapo_points){
  std::vector<Cyber::TrackState> ecal_track_states = track->getTrackStates("Ecal");
  for(int its=0; its<ecal_track_states.size(); its++){
    extrapo_points.push_back(ecal_track_states[its].referencePoint);
  }

  return StatusCode::SUCCESS;
}


StatusCode TrackMatchingAlg::CreateTrackAxis(vector<TVector3>& extrapo_points, std::vector<const Cyber::Calo1DCluster*>& localMaxCol,
                             Cyber::CaloHalfCluster* t_track_axis){                           
  if(localMaxCol.size()==0 || extrapo_points.size()==0)
    return StatusCode::SUCCESS;
  int t_slayer = localMaxCol[0]->getSlayer();
  int t_system = localMaxCol[0]->getSystem();
  
  if(t_system==Cyber::CaloUnit::System_Barrel){  // Barrel
    if(t_slayer==1){  // V plane (xy plane)
      for(int ipt=0; ipt<extrapo_points.size(); ipt++){
      for(int ilm=0; ilm<localMaxCol.size(); ilm++){
        // distance from the extrpolated point to the center of the local max bar
        TVector3 distance = extrapo_points[ipt] - localMaxCol[ilm]->getPos();
        if( TMath::Abs(distance.Z()) < (localMaxCol[ilm]->getBars()[0]->getBarLength())/2. &&
            distance.Perp() < Cyber::CaloUnit::barsize ) { 
          t_track_axis->addUnit(localMaxCol[ilm]);
        }
        else { continue; }
      }}
    }
    else{  // U plane (r-phi plane)
      for(int ipt=0; ipt<extrapo_points.size(); ipt++){
      for(int ilm=0; ilm<localMaxCol.size(); ilm++){
        TVector3 lm_pos = localMaxCol[ilm]->getPos();
        float barLength = localMaxCol[ilm]->getBars()[0]->getBarLength();
        if( fabs(extrapo_points[ipt].z()-lm_pos.z()) < Cyber::CaloUnit::barsize &&  
            fabs(extrapo_points[ipt].Phi()-lm_pos.Phi()) < barLength/2./Cyber::CaloUnit::ecal_innerR &&  
            fabs(extrapo_points[ipt].Perp()-lm_pos.Perp()) < Cyber::CaloUnit::barsize ){
          t_track_axis->addUnit(localMaxCol[ilm]);
        }
        else { continue; }
      }}
    }
  }
  else if(t_system==Cyber::CaloUnit::System_Endcap){ // Endcap
    if(t_slayer==0){ // U plane
      for(int ipt=0; ipt<extrapo_points.size(); ipt++){
      for(int ilm=0; ilm<localMaxCol.size(); ilm++){
        TVector3 lm_pos = localMaxCol[ilm]->getPos();
        float barLength = localMaxCol[ilm]->getBars()[0]->getBarLength();
        TVector3 distance = extrapo_points[ipt] - lm_pos;
        if( fabs(distance.z()) < Cyber::CaloUnit::barsize &&  
            fabs(distance.x()) < Cyber::CaloUnit::barsize &&
            fabs(distance.y()) < barLength/2.){
          t_track_axis->addUnit(localMaxCol[ilm]);
        }
        else { continue; }
      }}
    }
    else{ // V plane
      for(int ipt=0; ipt<extrapo_points.size(); ipt++){
      for(int ilm=0; ilm<localMaxCol.size(); ilm++){
        TVector3 lm_pos = localMaxCol[ilm]->getPos();
        float barLength = localMaxCol[ilm]->getBars()[0]->getBarLength();
        TVector3 distance = extrapo_points[ipt] - lm_pos;
        if( fabs(distance.z()) < Cyber::CaloUnit::barsize &&  
            fabs(distance.y()) < Cyber::CaloUnit::barsize &&
            fabs(distance.x()) < barLength/2.){  
          t_track_axis->addUnit(localMaxCol[ilm]);
        }
        else { continue; }
      }}
    }
  }
  

  return StatusCode::SUCCESS;
}


StatusCode TrackMatchingAlg::GetNearby(const std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>* p_HalfCluster, 
                                       const std::vector<TVector3>& extrapo_points, 
                                       std::vector<Cyber::CaloHalfCluster*>& t_nearbyHalfClusters, 
                                       std::vector<const Cyber::Calo1DCluster*>& t_nearbyLocalMax){


  if(p_HalfCluster->size()==0 || extrapo_points.size()==0)  return StatusCode::SUCCESS; 

  std::set<Cyber::CaloHalfCluster*> set_nearbyHalfClusters;
  int slayer = p_HalfCluster->at(0).get()->getSlayer();

  if(slayer==1){  // V plane
    for(int ihc=0; ihc<p_HalfCluster->size(); ihc++){
      int system = p_HalfCluster->at(ihc).get()->getBars()[0]->getSystem();
      std::vector<const Cyber::Calo1DCluster*> localMaxCol = p_HalfCluster->at(ihc).get()->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);
//cout<<"    HFclus #"<<ihc<<": system "<<system<<", localMax size "<<localMaxCol.size()<<endl;
      for(int ilm=0; ilm<localMaxCol.size(); ilm++){
      for(int ipt=0; ipt<extrapo_points.size(); ipt++){
        TVector3 distance(extrapo_points[ipt] - localMaxCol[ilm]->getPos());
        float barLength = localMaxCol[ilm]->getBars()[0]->getBarLength();
        if(system==Cyber::CaloUnit::System_Barrel){  // Barrel
          if(TMath::Abs(distance.Z()) < barLength/2. && 
              distance.Perp() < settings.map_floatPars["ConeNearByDistance"] ){  
            t_nearbyLocalMax.push_back(localMaxCol[ilm]);
            set_nearbyHalfClusters.insert(p_HalfCluster->at(ihc).get());
            break;
          }
        }
        else if(system==Cyber::CaloUnit::System_Endcap){
          if( TMath::Sqrt(distance.z()*distance.z() + distance.y()*distance.y()) < settings.map_floatPars["ConeNearByDistance"] && 
              fabs(distance.x()) < barLength/2.){
            t_nearbyLocalMax.push_back(localMaxCol[ilm]);
            set_nearbyHalfClusters.insert(p_HalfCluster->at(ihc).get());
            break;
          }
        }
      }}
    }
  }
  else if(slayer==0){
    for(int ihc=0; ihc<p_HalfCluster->size(); ihc++){
      int system = p_HalfCluster->at(ihc).get()->getBars()[0]->getSystem();
      std::vector<const Cyber::Calo1DCluster*> localMaxCol = p_HalfCluster->at(ihc).get()->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);
//cout<<"    HFclus #"<<ihc<<": system "<<system<<", localMax size "<<localMaxCol.size()<<endl;
      for(int ilm=0; ilm<localMaxCol.size(); ilm++){
      for(int ipt=0; ipt<extrapo_points.size(); ipt++){
        TVector3 lm_pos = localMaxCol[ilm]->getPos();
        float barLength = localMaxCol[ilm]->getBars()[0]->getBarLength();
        TVector3 distance = extrapo_points[ipt] - lm_pos;
        if(system==Cyber::CaloUnit::System_Barrel){  // Barrel
          if( fabs(extrapo_points[ipt].z()-lm_pos.z()) < settings.map_floatPars["ConeNearByDistance"] &&
              fabs(extrapo_points[ipt].Phi()-lm_pos.Phi()) < barLength/2./Cyber::CaloUnit::ecal_innerR &&
              fabs(extrapo_points[ipt].Perp()-lm_pos.Perp()) < settings.map_floatPars["ConeNearByDistance"] ){
            t_nearbyLocalMax.push_back(localMaxCol[ilm]);
            set_nearbyHalfClusters.insert(p_HalfCluster->at(ihc).get());
            break;
          }
        }
        else if(system==Cyber::CaloUnit::System_Endcap){
          if( TMath::Sqrt(distance.z()*distance.z() + distance.x()*distance.x()) < settings.map_floatPars["ConeNearByDistance"] &&  
              fabs(distance.y()) < barLength/2.){
            t_nearbyLocalMax.push_back(localMaxCol[ilm]);
            set_nearbyHalfClusters.insert(p_HalfCluster->at(ihc).get());
            break;
          }
        }
      }}
    }
  }
  
  t_nearbyHalfClusters.assign(set_nearbyHalfClusters.begin(), set_nearbyHalfClusters.end());
  return StatusCode::SUCCESS; 
}


StatusCode TrackMatchingAlg::LongiConeLinking(const std::vector<TVector3>& extrapo_points,  
                                              std::vector<const Cyber::Calo1DCluster*>& nearbyLocalMax, 
                                              std::vector<const Cyber::Calo1DCluster*>& cone_axis){
  if(nearbyLocalMax.size()==0 || extrapo_points.size()==0) return StatusCode::SUCCESS;
  int slayer = nearbyLocalMax[0]->getSlayer();
  // int min_point = settings.map_intPars["Max_Seed_Point"];
  // if (extrapo_points.size()<min_point) min_point = extrapo_points.size();
  
  std::vector<const Cyber::Calo1DCluster*> barrel_localMax, endcap_localMax;
  for(int ilm=0; ilm<nearbyLocalMax.size(); ilm++){
    if(nearbyLocalMax[ilm]->getSystem()==Cyber::CaloUnit::System_Barrel) barrel_localMax.push_back(nearbyLocalMax[ilm]);
    else if(nearbyLocalMax[ilm]->getSystem()==Cyber::CaloUnit::System_Endcap) endcap_localMax.push_back(nearbyLocalMax[ilm]);
  }
  // Seed finding for barrel
  std::vector<const Cyber::Calo1DCluster*> barrel_cone_axis;
  if(slayer==1){  // If V plane
    // for(int ip=0; ip<min_point; ip++){
    for(int ip=0; ip<extrapo_points.size(); ip++){
      double min_distance = 99999;
      int seed_candidate_index = -1;

      for(int il=0;il<barrel_localMax.size(); il++){
        TVector3 distance = extrapo_points[ip] - barrel_localMax[il]->getPos();
        double distance_2d = distance.Perp();
        if(TMath::Abs(distance.Z()) < (barrel_localMax[il]->getBars()[0]->getBarLength())/2.
           && distance_2d < settings.map_floatPars["ConeSeedDistance"]
           && distance_2d < min_distance)
        {
          seed_candidate_index = il;
          min_distance = distance_2d;
        }
      }
      if (seed_candidate_index<0) continue;

      barrel_cone_axis.push_back(barrel_localMax[seed_candidate_index]);
      barrel_localMax.erase(barrel_localMax.begin() + seed_candidate_index);
      break;
    }
  }
  else{  // If U plane
    // for(int ip=0; ip<min_point; ip++){
    for(int ip=0; ip<extrapo_points.size(); ip++){
      double min_distance = 99999;
      int seed_candidate_index = -1;

      for(int il=0;il<barrel_localMax.size(); il++){
        TVector3 lm_pos = barrel_localMax[il]->getPos();
        float barLength = barrel_localMax[il]->getBars()[0]->getBarLength();
        float distance_2d = sqrt( pow(extrapo_points[ip].z()-lm_pos.z(), 2) + pow(extrapo_points[ip].Perp()-lm_pos.Perp(), 2) );
        if( fabs(extrapo_points[ip].Phi()-lm_pos.Phi()) < barLength/2./Cyber::CaloUnit::ecal_innerR && 
            distance_2d < settings.map_floatPars["ConeSeedDistance"] &&
            distance_2d < min_distance)
        {
          seed_candidate_index = il;
          min_distance = distance_2d;
        }
      }
      if (seed_candidate_index<0) continue;

      barrel_cone_axis.push_back(barrel_localMax[seed_candidate_index]);
      barrel_localMax.erase(barrel_localMax.begin() + seed_candidate_index);
      break;
    }
  }

  // Seed finding for endcap
  std::vector<const Cyber::Calo1DCluster*> endcap_cone_axis;
  if(slayer==0){  // U plane (bars parralel to y-axis)
    // for(int ip=0; ip<min_point; ip++){
    for(int ip=0; ip<extrapo_points.size(); ip++){
      double min_distance = 99999;
      int seed_candidate_index = -1;

      for(int il=0;il<endcap_localMax.size(); il++){
        TVector3 lm_pos = endcap_localMax[il]->getPos();
        float barLength = endcap_localMax[il]->getBars()[0]->getBarLength();
        TVector3 distance = extrapo_points[ip] - lm_pos;
        double distance_2d = sqrt( distance.z()*distance.z() + distance.x()*distance.x() );
        if( fabs(distance.y()) < barLength/2. && 
            distance_2d < settings.map_floatPars["ConeSeedDistance"] &&
            distance_2d < min_distance)
        {
          seed_candidate_index = il;
          min_distance = distance_2d;
        }
      }
      if (seed_candidate_index<0) continue;

      endcap_cone_axis.push_back(endcap_localMax[seed_candidate_index]);
      endcap_localMax.erase(endcap_localMax.begin() + seed_candidate_index);
      break;
    }
  }
  else{  // V plane (bars parralel to x-axis)
    // for(int ip=0; ip<min_point; ip++){
    for(int ip=0; ip<extrapo_points.size(); ip++){
      double min_distance = 99999;
      int seed_candidate_index = -1;

      for(int il=0;il<endcap_localMax.size(); il++){
        TVector3 lm_pos = endcap_localMax[il]->getPos();
        float barLength = endcap_localMax[il]->getBars()[0]->getBarLength();
        TVector3 distance = extrapo_points[ip] - lm_pos;
        double distance_2d = sqrt( distance.z()*distance.z() + distance.y()*distance.y() );
        if( fabs(distance.x()) < barLength/2. && 
            distance_2d < settings.map_floatPars["ConeSeedDistance"] &&
            distance_2d < min_distance)
        {
          seed_candidate_index = il;
          min_distance = distance_2d;
        }
      }
      if (seed_candidate_index<0) continue;

      endcap_cone_axis.push_back(endcap_localMax[seed_candidate_index]);
      endcap_localMax.erase(endcap_localMax.begin() + seed_candidate_index);
      break;
    }

  }

  if (barrel_cone_axis.size()==0 && endcap_cone_axis.size()==0) return StatusCode::SUCCESS;

  // Linking for barrel
  while(barrel_localMax.size()>0){
    if(barrel_cone_axis.size()==0) break;
    const Cyber::Calo1DCluster* shower_in_axis = barrel_cone_axis.back();
    if(!shower_in_axis) break; 
    if(isStopLinking(extrapo_points, shower_in_axis)) break;

    double min_distance = 9999;
    int shower_candidate_index = -1;

    for(int il=0; il<barrel_localMax.size(); il++){
      TVector2 relR = GetProjectedRelR(shower_in_axis, barrel_localMax[il]);  //Return vec: 1->2.
      TVector2 clusaxis = GetProjectedAxis(extrapo_points, shower_in_axis);

      double delta_phi = relR.DeltaPhi(clusaxis);
      double delta_distance = (relR - (clusaxis*2)).Mod();

      if( delta_phi<settings.map_floatPars["th_ConeTheta"] 
          && relR.Mod()<settings.map_floatPars["th_ConeR"] 
          && delta_distance<min_distance){
        shower_candidate_index = il;
        min_distance = delta_distance;
      }
    }
    if (shower_candidate_index<0) break;

    barrel_cone_axis.push_back(barrel_localMax[shower_candidate_index]);
    barrel_localMax.erase(barrel_localMax.begin() + shower_candidate_index);
  }


  // Linking for endcap
  while(endcap_localMax.size()>0){
    if(endcap_cone_axis.size()==0) break;
    const Cyber::Calo1DCluster* shower_in_axis = endcap_cone_axis.back();
    if(!shower_in_axis) break; 
    if(isStopLinking(extrapo_points, shower_in_axis)) break;

    double min_distance = 9999;
    int shower_candidate_index = -1;

    for(int il=0; il<endcap_localMax.size(); il++){
      TVector2 relR = GetProjectedRelR(shower_in_axis, endcap_localMax[il]);  //Return vec: 1->2.
      TVector2 clusaxis = GetProjectedAxis(extrapo_points, shower_in_axis);

      double delta_phi = relR.DeltaPhi(clusaxis);
      double delta_distance = (relR - (clusaxis*2)).Mod();

      if( delta_phi<settings.map_floatPars["th_ConeTheta"] 
          && relR.Mod()<settings.map_floatPars["th_ConeR"] 
          && delta_distance<min_distance){
        shower_candidate_index = il;
        min_distance = delta_distance;
      }
    }
    if (shower_candidate_index<0) break;

    endcap_cone_axis.push_back(endcap_localMax[shower_candidate_index]);
    endcap_localMax.erase(endcap_localMax.begin() + shower_candidate_index);
  }

  cone_axis.insert(cone_axis.end(), barrel_cone_axis.begin(), barrel_cone_axis.end());
  cone_axis.insert(cone_axis.end(), endcap_cone_axis.begin(), endcap_cone_axis.end());

  // nearbyLocalMax.clear();
  // nearbyLocalMax.insert(nearbyLocalMax.end(), barrel_localMax.begin(), barrel_localMax.end());
  // nearbyLocalMax.insert(nearbyLocalMax.end(), endcap_localMax.begin(), endcap_localMax.end());

  return StatusCode::SUCCESS;
}


bool TrackMatchingAlg::isStopLinking( const std::vector<TVector3>& extrapo_points, 
                                      const Cyber::Calo1DCluster* final_cone_hit){

  int slayer = final_cone_hit->getSlayer();
  int system = final_cone_hit->getSystem();
  if(system==Cyber::CaloUnit::System_Barrel){
    if(slayer==1){
      TVector3 f_distance = extrapo_points.back() - final_cone_hit->getPos();
      double f_distance_2d = f_distance.Perp(); 
      for(int i=0; i<extrapo_points.size(); i++){
        TVector3 distance = extrapo_points[i] - final_cone_hit->getPos();
        double distance_2d = distance.Perp(); 
        if (distance_2d < f_distance_2d) return false;
      }
    }
    else{
      TVector3 lm_pos = final_cone_hit->getPos();
      float f_distance_2d = sqrt( pow(extrapo_points.back().z()-lm_pos.z(), 2) + pow(extrapo_points.back().Perp()-lm_pos.Perp(), 2) );
      for(int i=0; i<extrapo_points.size(); i++){
        float distance_2d = sqrt( pow(extrapo_points[i].z()-lm_pos.z(), 2) + pow(extrapo_points[i].Perp()-lm_pos.Perp(), 2) );
        if (distance_2d < f_distance_2d) return false;
      }
    }
  }
  else if(system==Cyber::CaloUnit::System_Endcap){
    if(slayer==0){
      TVector3 lm_pos = final_cone_hit->getPos();
      TVector3 f_distance = extrapo_points.back() - lm_pos;
      double f_distance_2d = sqrt( f_distance.z()*f_distance.z() + f_distance.x()*f_distance.x() );
      for(int i=0; i<extrapo_points.size(); i++){
        TVector3 distance = extrapo_points[i] - lm_pos;
        double distance_2d = sqrt( distance.z()*distance.z() + distance.x()*distance.x() );
        if (distance_2d < f_distance_2d) return false;
      }
    }
    else{
      TVector3 lm_pos = final_cone_hit->getPos();
      TVector3 f_distance = extrapo_points.back() - lm_pos;
      double f_distance_2d = sqrt( f_distance.z()*f_distance.z() + f_distance.y()*f_distance.y() );
      for(int i=0; i<extrapo_points.size(); i++){
        TVector3 distance = extrapo_points[i] - lm_pos;
        double distance_2d = sqrt( distance.z()*distance.z() + distance.y()*distance.y() );
        if (distance_2d < f_distance_2d) return false;
      }
    }
  }
  
  return true;   
}


TVector2 TrackMatchingAlg::GetProjectedRelR( const Cyber::Calo1DCluster* m_shower1, const Cyber::Calo1DCluster* m_shower2 ){
  if(m_shower1->getSystem()==Cyber::CaloUnit::System_Barrel){  // For Barrel
    if(m_shower1->getSlayer()==1){ //For V-bars
      TVector3 vec = m_shower2->getPos() - m_shower1->getPos();
      TVector2 vec2d(vec.x(), vec.y());
      return vec2d;
    }
    else{  //For U-bars
      TVector3 vec = m_shower2->getPos() - m_shower1->getPos();
      TVector2 vec2d(vec.Perp(), vec.z());
      return vec2d;
    }
  }
  else if(m_shower1->getSystem()==Cyber::CaloUnit::System_Endcap){  // For Endcap
    if(m_shower1->getSlayer()==0){ //For U-bars
      TVector3 vec = m_shower2->getPos() - m_shower1->getPos();
      TVector2 vec2d(vec.z(), vec.x());
      return vec2d;
    }
    else{  //For V-bars
      TVector3 vec = m_shower2->getPos() - m_shower1->getPos();
      TVector2 vec2d(vec.z(), vec.y());
      return vec2d;
    }
  }  
}


TVector2 TrackMatchingAlg::GetProjectedAxis(const std::vector<TVector3>& extrapo_points, const Cyber::Calo1DCluster* m_shower){
  int min_index=0;
  TVector2 distance(999., 999.);
  if(m_shower->getSystem()==Cyber::CaloUnit::System_Barrel){
    if( m_shower->getSlayer()==1 ){  // V plane
      for(int i=0; i<extrapo_points.size(); i++){
        TVector2 t_distance(m_shower->getPos().x()-extrapo_points[i].x(), m_shower->getPos().y()-extrapo_points[i].y());
        if(t_distance.Mod()<distance.Mod()){
          distance = t_distance;
          min_index = i;
        }
      }

      if(min_index < extrapo_points.size()-1){
        TVector2 axis(extrapo_points[min_index+1].x()-extrapo_points[min_index].x(), extrapo_points[min_index+1].y()-extrapo_points[min_index].y());
        return axis;
      }else{
        TVector2 axis(extrapo_points[min_index].x()-extrapo_points[min_index-1].x(), extrapo_points[min_index].y()-extrapo_points[min_index-1].y());
        return axis;
      }

    }else{  // U plane
      for(int i=0; i<extrapo_points.size(); i++){
        TVector3 dist3d = m_shower->getPos() - extrapo_points[i];
        //dist3d.RotateZ( TMath::Pi()/4.*(6-m_shower->getTowerID()[0][0]) );
        TVector2 t_distance(dist3d.Perp(), dist3d.z());
        if(t_distance.Mod()<distance.Mod()){
          distance = t_distance;
          min_index = i;
        }
      }

      if(min_index < extrapo_points.size()-1){
        TVector3 vec = extrapo_points[min_index+1] - extrapo_points[min_index];
        TVector2 axis(vec.Perp(), vec.z());
        return axis;
      }else{
        TVector3 vec = extrapo_points[min_index] - extrapo_points[min_index-1];
        TVector2 axis(vec.Perp(), vec.z());
        return axis;
      }
    }
  }
  else if (m_shower->getSystem()==Cyber::CaloUnit::System_Endcap){
    if( m_shower->getSlayer()==0 ){  // U plane
      for(int i=0; i<extrapo_points.size(); i++){
        TVector3 dist3d = m_shower->getPos() - extrapo_points[i];
        TVector2 t_distance(dist3d.z(), dist3d.x());
        if(t_distance.Mod()<distance.Mod()){
          distance = t_distance;
          min_index = i;
        }
      }

      if(min_index < extrapo_points.size()-1){
        TVector3 vec = extrapo_points[min_index+1] - extrapo_points[min_index];
        TVector2 axis(vec.z(), vec.x());
        return axis;
      }else{
        TVector3 vec = extrapo_points[min_index] - extrapo_points[min_index-1];
        TVector2 axis(vec.z(), vec.x());
        return axis;
      }

    }else{  // V plane
      for(int i=0; i<extrapo_points.size(); i++){
        TVector3 dist3d = m_shower->getPos() - extrapo_points[i];
        TVector2 t_distance(dist3d.z(), dist3d.y());
        if(t_distance.Mod()<distance.Mod()){
          distance = t_distance;
          min_index = i;
        }
      }

      if(min_index < extrapo_points.size()-1){
        TVector3 vec = extrapo_points[min_index+1] - extrapo_points[min_index];
        TVector2 axis(vec.z(), vec.y());
        return axis;
      }else{
        TVector3 vec = extrapo_points[min_index] - extrapo_points[min_index-1];
        TVector2 axis(vec.z(), vec.y());
        return axis;
      }
    }
  }
  

  
}


StatusCode TrackMatchingAlg::CreatConeAxis(CyberDataCol& m_datacol, Cyber::Track* track, std::vector<Cyber::CaloHalfCluster*>& nearbyHalfClusters, 
                                           std::vector<const Cyber::Calo1DCluster*>& cone_axis){
  if(nearbyHalfClusters.size()==0 || cone_axis.size()==0) return StatusCode::SUCCESS; 

  for(int ihc=0; ihc<nearbyHalfClusters.size(); ihc++){
    std::vector<const Cyber::Calo1DCluster*> localMaxCol = nearbyHalfClusters[ihc]->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);
    // Track axis candidate.
    // Cyber::CaloHalfCluster* t_track_axis = new Cyber::CaloHalfCluster();
    std::shared_ptr<Cyber::CaloHalfCluster> t_track_axis = std::make_shared<Cyber::CaloHalfCluster>();
    for(int ica=0; ica<cone_axis.size(); ica++){
      if( find(localMaxCol.begin(), localMaxCol.end(), cone_axis[ica]) != localMaxCol.end()){
        t_track_axis->addUnit(cone_axis[ica]);
      }
    }

    // If the track does not match the Halfcluster, the track axis candidate will have no 1DCluster
    if(t_track_axis->getCluster().size()==0)
      continue;
    
    t_track_axis->addAssociatedTrack(track);
    t_track_axis->setType(10000); //Track-type axis. 
    
    if(nearbyHalfClusters[ihc]->getSlayer()==1){
      track->addAssociatedHalfClusterV( nearbyHalfClusters[ihc] );
    }
    else{
      track->addAssociatedHalfClusterU( nearbyHalfClusters[ihc] );
    }
    m_datacol.map_HalfCluster["bkHalfCluster"].push_back(t_track_axis);
    nearbyHalfClusters[ihc]->addHalfCluster(settings.map_stringPars["OutputLongiClusName"], t_track_axis.get());
    
    
    
  }

  return StatusCode::SUCCESS; 
}

#endif
