#ifndef TRACKMATCHING_C
#define TRACKMATCHING_C

#include "Algorithm/TrackMatchingAlg.h"


StatusCode TrackMatchingAlg::ReadSettings(PandoraPlus::Settings& m_settings){
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
  //if(settings.map_intPars.find("Nmodule")==settings.map_intPars.end()) 
  //  settings.map_intPars["Nmodule"] = 8;
  if(settings.map_stringPars.find("ReadinLocalMaxName")==settings.map_stringPars.end())
    settings.map_stringPars["ReadinLocalMaxName"] = "AllLocalMax";
  if(settings.map_stringPars.find("OutputLongiClusName")==settings.map_stringPars.end()) 
    settings.map_stringPars["OutputLongiClusName"] = "TrackAxis"; 

  return StatusCode::SUCCESS;
}


StatusCode TrackMatchingAlg::Initialize( PandoraPlusDataCol& m_datacol ){
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


StatusCode TrackMatchingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){
  std::cout << "---oooOO0OOooo---Excuting TrackMatchingAlg---oooOO0OOooo---"<<std::endl;
  // Associate tracks to HalfClusters.
  // This association is a many-to-many relationship: 
  //    One HalCluster may have multiple tracks; 
  //    One track may pass through multiple HalfClusters.
  for(int itrk=0; itrk<m_TrackCol.size(); itrk++){  // loop tracks
    if(m_TrackCol[itrk]->getTrackStates("Ecal").size()==0) continue;

    // Get extrapolated points of the track. These points are sorted by the track
    std::vector<TVector3> extrapo_points;
    GetExtrpoECALPoints(m_TrackCol[itrk], extrapo_points);
    if(extrapo_points.size()==0) continue;

    double pT = TMath::Abs(1. / m_TrackCol[itrk]->getTrackStates("Ecal")[0].Kappa);
    if (pT >= settings.map_floatPars["ConeMatchingCut_pT"]){
      std::cout << "For track " << itrk << ", pT = " << pT <<", match directly" << std::endl;

      for(int ihc=0; ihc<p_HalfClusterV->size(); ihc++){  // loop HalfClusterV
        // Get local max of the HalfCluster
        std::vector<const PandoraPlus::Calo1DCluster*> localMaxColV = p_HalfClusterV->at(ihc).get()->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

        // Track axis candidate.
        std::shared_ptr<PandoraPlus::CaloHalfCluster> t_track_axis = std::make_shared<PandoraPlus::CaloHalfCluster>();
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
        std::vector<const PandoraPlus::Calo1DCluster*> localMaxColU = p_HalfClusterU->at(ihc).get()->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

        // Track axis candidate.
        std::shared_ptr<PandoraPlus::CaloHalfCluster> t_track_axis = std::make_shared<PandoraPlus::CaloHalfCluster>();
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
      std::cout << "For track " << itrk << ", pT = " << pT <<", using Cone method" << std::endl;

      // Get local max and HalfCluster near the extrapolated points
      std::vector<PandoraPlus::CaloHalfCluster*> t_nearbyHalfClustersV;  t_nearbyHalfClustersV.clear();
      std::vector<PandoraPlus::CaloHalfCluster*> t_nearbyHalfClustersU;  t_nearbyHalfClustersU.clear();
      std::vector<const PandoraPlus::Calo1DCluster*> t_nearbyLocalMaxV;     t_nearbyLocalMaxV.clear();
      std::vector<const PandoraPlus::Calo1DCluster*> t_nearbyLocalMaxU;     t_nearbyLocalMaxU.clear();
      GetNearby(p_HalfClusterV, extrapo_points, t_nearbyHalfClustersV, t_nearbyLocalMaxV);
      GetNearby(p_HalfClusterU, extrapo_points, t_nearbyHalfClustersU, t_nearbyLocalMaxU);

      // V plane
      std::vector<const PandoraPlus::Calo1DCluster*> t_cone_axisV; t_cone_axisV.clear();
      LongiConeLinking(extrapo_points, t_nearbyLocalMaxV, t_cone_axisV);
      CreatConeAxis(m_datacol, m_TrackCol[itrk], t_nearbyHalfClustersV, t_cone_axisV);

      // U plane
      // Sort local max by their modules
      std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > m_orderedLocalMaxU;  // key: module of the bar
      m_orderedLocalMaxU.clear();
      for(int is=0; is<t_nearbyLocalMaxU.size(); is++)
        m_orderedLocalMaxU[t_nearbyLocalMaxU[is]->getTowerID()[0][0]].push_back(t_nearbyLocalMaxU[is]);
      // linking 
      std::vector<const PandoraPlus::Calo1DCluster*> merged_cone_axisU; merged_cone_axisU.clear();
      for (auto it = m_orderedLocalMaxU.begin(); it != m_orderedLocalMaxU.end(); ++it){
        std::vector<const PandoraPlus::Calo1DCluster*> moduled_localMaxU = it->second;
        std::vector<const PandoraPlus::Calo1DCluster*> t_cone_axisU; t_cone_axisU.clear();
        LongiConeLinking(extrapo_points, moduled_localMaxU, t_cone_axisU);
        merged_cone_axisU.insert(merged_cone_axisU.end(), t_cone_axisU.begin(), t_cone_axisU.end());
      }
      CreatConeAxis(m_datacol, m_TrackCol[itrk], t_nearbyHalfClustersU, merged_cone_axisU);
    }

  

  }  // end loop tracks

  //Loop track to check the associated cluster: merge clusters if they are associated to the same track.
  std::vector<PandoraPlus::CaloHalfCluster*> tmp_deleteClus; tmp_deleteClus.clear();
  for(auto &itrk : m_TrackCol){
    std::vector<PandoraPlus::CaloHalfCluster*> m_matchedUCol = itrk->getAssociatedHalfClustersU();
    std::vector<PandoraPlus::CaloHalfCluster*> m_matchedVCol = itrk->getAssociatedHalfClustersV();

    // std::cout << "yyy: Before merge, m_matchedVCol.size() = " << m_matchedVCol.size() << std::endl;
    for(int imc=0; imc<m_matchedVCol.size(); imc++){
      // std::cout<<"  yyy: m_matchedVCol["<<imc<<"]->getCluster().size()="<<m_matchedVCol[imc]->getCluster().size()<<std::endl;
      // std::cout<<"  yyy: m_matchedVCol["<<imc<<"]->getAssociatedTracks().size()=" << m_matchedVCol[imc]->getAssociatedTracks().size() << std::endl;
      int N_trk_axis = m_matchedVCol[imc]->getHalfClusterMap()["TrackAxis"].size() ;
      // std::cout<<"  yyy: m_matchedVCol["<<imc<<"]->getHalfClusterMap()[TrackAxis].size()="<< N_trk_axis << std::endl;
      for(int itk=0; itk<N_trk_axis; itk++){
        // std::cout<<"    yyy: for TrackAxis " << itk 
        // << ", Nlm = " << m_matchedVCol[imc]->getHalfClusterMap()["TrackAxis"][itk]->getCluster().size()
        // << ", N trk = " << m_matchedVCol[imc]->getHalfClusterMap()["TrackAxis"][itk]->getAssociatedTracks().size() << std::endl;
      }
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

    /////////////////////////////////////////////
    // if(m_matchedVCol.size()>0){
    //   // std::cout<<"yyy: After merge, m_matchedVCol[0]->getCluster().size()=" << m_matchedVCol[0]->getCluster().size() << std::endl;
    //   // std::cout<<"yyy: After merge, m_matchedVCol[0]->getAssociatedTracks().size()=" << m_matchedVCol[0]->getAssociatedTracks().size() << std::endl;
    //   int N_trk_axis = m_matchedVCol[0]->getHalfClusterMap()["TrackAxis"].size() ;
    //   // std::cout<<"yyy: After merge, m_matchedVCol[0]->getHalfClusterMap()[TrackAxis].size()="<< N_trk_axis << std::endl;
    //   for(int nn=0; nn<N_trk_axis; nn++){
    //     std::cout<<"yyy: for TrackAxis " << nn 
    //     << ", Nlm = " << m_matchedVCol[0]->getHalfClusterMap()["TrackAxis"][nn]->getCluster().size()
    //     << ", N trk = " << m_matchedVCol[0]->getHalfClusterMap()["TrackAxis"][nn]->getAssociatedTracks().size() << std::endl;

    //   }
    // }
    /////////////////////////////////////////////
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


StatusCode TrackMatchingAlg::GetExtrpoECALPoints(const PandoraPlus::Track* track, std::vector<TVector3>& extrapo_points){
  std::vector<PandoraPlus::TrackState> ecal_track_states = track->getTrackStates("Ecal");
  for(int its=0; its<ecal_track_states.size(); its++){
    extrapo_points.push_back(ecal_track_states[its].referencePoint);
  }

  return StatusCode::SUCCESS;
}


StatusCode TrackMatchingAlg::CreateTrackAxis(vector<TVector3>& extrapo_points, std::vector<const PandoraPlus::Calo1DCluster*>& localMaxCol,
                             PandoraPlus::CaloHalfCluster* t_track_axis){                           
  if(localMaxCol.size()==0 || extrapo_points.size()==0)
    return StatusCode::SUCCESS;
  int t_slayer = localMaxCol[0]->getSlayer();
  
  if(t_slayer==1){  // V plane (xy plane)
    for(int ipt=0; ipt<extrapo_points.size(); ipt++){
    for(int ilm=0; ilm<localMaxCol.size(); ilm++){
      // distance from the extrpolated point to the center of the local max bar
      TVector3 distance = extrapo_points[ipt] - localMaxCol[ilm]->getPos();
      if(TMath::Abs(distance.Z()) < (localMaxCol[ilm]->getBars()[0]->getBarLength())/2.
         //&& TMath::Sqrt(distance.X()*distance.X() + distance.Y()*distance.Y()) < settings.map_floatPars["localmax_area"] ) { 
         && TMath::Sqrt(distance.X()*distance.X() + distance.Y()*distance.Y()) < PandoraPlus::CaloUnit::barsize ) { 
        t_track_axis->addUnit(localMaxCol[ilm]);
      }
      else { continue; }
    }}
  }
  else{  // U plane (r-phi plane)
    for(int ipt=0; ipt<extrapo_points.size(); ipt++){
    for(int ilm=0; ilm<localMaxCol.size(); ilm++){
      // distance from the extrpolated point to the center of the local max bar
      TVector3 distance = extrapo_points[ipt] - localMaxCol[ilm]->getPos();
      // rotate the distance to module 6
      distance.RotateZ( TMath::Pi()/4.*(6-localMaxCol[ilm]->getTowerID()[0][0]) );
      if(TMath::Abs(distance.Y()) < (localMaxCol[ilm]->getBars()[0]->getBarLength())/2.
         //&& TMath::Sqrt(distance.X()*distance.X() + distance.Z()*distance.Z()) < settings.map_floatPars["localmax_area"] ) { 
         && TMath::Sqrt(distance.X()*distance.X() + distance.Z()*distance.Z()) < PandoraPlus::CaloUnit::barsize ) { 
        t_track_axis->addUnit(localMaxCol[ilm]);
      }
      else { continue; }
    }}
  }

//std::cout << "end calling CreateTrackAxis()" << std::endl;
  return StatusCode::SUCCESS;
}


StatusCode TrackMatchingAlg::GetNearby(const std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>* p_HalfCluster, 
                                       const std::vector<TVector3>& extrapo_points, 
                                       std::vector<PandoraPlus::CaloHalfCluster*>& t_nearbyHalfClusters, 
                                       std::vector<const PandoraPlus::Calo1DCluster*>& t_nearbyLocalMax){
  // std::cout << "calling TrackMatchingAlg::GetNearby()" << std::endl;

  if(p_HalfCluster->size()==0 || extrapo_points.size()==0)  return StatusCode::SUCCESS; 

  std::set<PandoraPlus::CaloHalfCluster*> set_nearbyHalfClusters;
  int slayer = p_HalfCluster->at(0).get()->getSlayer();
  if(slayer==1){  // V plane
    for(int ihc=0; ihc<p_HalfCluster->size(); ihc++){
      std::vector<const PandoraPlus::Calo1DCluster*> localMaxCol = p_HalfCluster->at(ihc).get()->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);
      for(int ilm=0; ilm<localMaxCol.size(); ilm++){
      for(int ipt=0; ipt<extrapo_points.size(); ipt++){
        TVector3 distance(extrapo_points[ipt] - localMaxCol[ilm]->getPos());
        if(TMath::Abs(distance.Z()) < (localMaxCol[ilm]->getBars()[0]->getBarLength())/2. && 
            TMath::Sqrt(distance.X()*distance.X() + distance.Y()*distance.Y()) < settings.map_floatPars["ConeNearByDistance"] ){  
          t_nearbyLocalMax.push_back(localMaxCol[ilm]);
          set_nearbyHalfClusters.insert(p_HalfCluster->at(ihc).get());
          break;
        }
      }}
      
    }
  }else{  // U plane
    for(int ihc=0; ihc<p_HalfCluster->size(); ihc++){
      std::vector<const PandoraPlus::Calo1DCluster*> localMaxCol = p_HalfCluster->at(ihc).get()->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);
      for(int ilm=0; ilm<localMaxCol.size(); ilm++){
      for(int ipt=0; ipt<extrapo_points.size(); ipt++){
        TVector3 distance(extrapo_points[ipt] - localMaxCol[ilm]->getPos());
        distance.RotateZ( TMath::Pi()/4.*(6-localMaxCol[ilm]->getTowerID()[0][0]) );
        if(TMath::Abs(distance.Y()) < (localMaxCol[ilm]->getBars()[0]->getBarLength())/2. && 
          TMath::Sqrt(distance.X()*distance.X() + distance.Z()*distance.Z()) < settings.map_floatPars["ConeNearByDistance"] ){
          t_nearbyLocalMax.push_back(localMaxCol[ilm]);
          set_nearbyHalfClusters.insert(p_HalfCluster->at(ihc).get());
          break;
        }
      }}
    }
  }
  t_nearbyHalfClusters.assign(set_nearbyHalfClusters.begin(), set_nearbyHalfClusters.end());

  return StatusCode::SUCCESS; 
}


StatusCode TrackMatchingAlg::LongiConeLinking(const std::vector<TVector3>& extrapo_points,  
                                              std::vector<const PandoraPlus::Calo1DCluster*>& nearbyLocalMax, 
                                              std::vector<const PandoraPlus::Calo1DCluster*>& cone_axis){
  // std::cout<<"yyy: calling longiConeLinking()"<<std::endl;
  if(nearbyLocalMax.size()==0 || extrapo_points.size()==0) return StatusCode::SUCCESS;

  // Seed finding
  int slayer = nearbyLocalMax[0]->getSlayer();
  
  // int min_point = settings.map_intPars["Max_Seed_Point"];
  // if (extrapo_points.size()<min_point) min_point = extrapo_points.size();
  
  if(slayer==1){  // If V plane
    // for(int ip=0; ip<min_point; ip++){
    for(int ip=0; ip<extrapo_points.size(); ip++){
      double min_distance = 99999;
      int seed_candidate_index = -1;

      for(int il=0;il<nearbyLocalMax.size(); il++){
        TVector3 distance = extrapo_points[ip] - nearbyLocalMax[il]->getPos();
        double distance_2d = TMath::Sqrt(distance.X()*distance.X() + distance.Y()*distance.Y());
        if(TMath::Abs(distance.Z()) < (nearbyLocalMax[il]->getBars()[0]->getBarLength())/2.
           && distance_2d < settings.map_floatPars["ConeSeedDistance"]
           && distance_2d < min_distance)
        {
          seed_candidate_index = il;
          min_distance = distance_2d;
        }
      }
      if (seed_candidate_index<0) continue;

      cone_axis.push_back(nearbyLocalMax[seed_candidate_index]);
      nearbyLocalMax.erase(nearbyLocalMax.begin() + seed_candidate_index);
      break;
    }
  }
  else{  // If U plane
    // for(int ip=0; ip<min_point; ip++){
    for(int ip=0; ip<extrapo_points.size(); ip++){
      double min_distance = 99999;
      int seed_candidate_index = -1;

      for(int il=0;il<nearbyLocalMax.size(); il++){
        TVector3 distance = extrapo_points[ip] - nearbyLocalMax[il]->getPos();
        distance.RotateZ( TMath::Pi()/4.*(6-nearbyLocalMax[il]->getTowerID()[0][0]) );
        double distance_2d = TMath::Sqrt(distance.X()*distance.X() + distance.Z()*distance.Z());
        if(TMath::Abs(distance.Y()) < (nearbyLocalMax[il]->getBars()[0]->getBarLength())/2.
           && distance_2d < settings.map_floatPars["ConeSeedDistance"]
           && distance_2d < min_distance)
        {
          seed_candidate_index = il;
          min_distance = distance_2d;
        }
      }
      if (seed_candidate_index<0) continue;

      cone_axis.push_back(nearbyLocalMax[seed_candidate_index]);
      nearbyLocalMax.erase(nearbyLocalMax.begin() + seed_candidate_index);
      break;
    }
  }

  if (cone_axis.size() == 0) return StatusCode::SUCCESS;

  // std::cout<<"  yyy: after seed-finding, cone_axis.size()="<<cone_axis.size()<<endl;

  // Linking
  while(nearbyLocalMax.size()>0){
    // std::cout<<"  yyy: nearbyLocalMax.size()="<<nearbyLocalMax.size()<<", linking it!"<<std::endl;
    const PandoraPlus::Calo1DCluster* shower_in_axis = cone_axis.back();
    if(!shower_in_axis) break; 
    if(isStopLinking(extrapo_points, shower_in_axis)) break;

    // std::cout<<"  yyy: looking for a lm to link"<<std::endl;

    double min_delta = 9999;
    int shower_candidate_index = -1;

    for(int il=0; il<nearbyLocalMax.size(); il++){
      TVector2 relR = GetProjectedRelR(shower_in_axis, nearbyLocalMax[il]);  //Return vec: 1->2.
      TVector2 clusaxis = GetProjectedAxis(extrapo_points, shower_in_axis);

      double delta_phi = relR.DeltaPhi(clusaxis);
      double delta_distance = (relR - (clusaxis*2)).Mod();

      // std::cout<<"    yyy: for nearbyLocalMax["<<il<<"], "<<std::endl
      //          <<"         shower_in_axis=("<<shower_in_axis->getPos().x()<<", "
      //          <<shower_in_axis->getPos().y()<<", "
      //          <<shower_in_axis->getPos().z()<<")"<<std::endl
      //          <<"         nearbyLocalMax=("<<nearbyLocalMax[il]->getPos().x()<<", "
      //          <<nearbyLocalMax[il]->getPos().y()<<", "
      //          <<nearbyLocalMax[il]->getPos().z()<<")"<<std::endl
      //          <<"         relR = ("<<relR.X()<<", "<<relR.Y()<<")"<<std::endl
      //          <<"         clusaxis = ("<<clusaxis.X()<<", "<<clusaxis.Y()<<")"<<std::endl
      //          <<"         delta_phi = "<<delta_phi<<", delta_distance="<<delta_distance<<std::endl;
      
      if( delta_phi<settings.map_floatPars["th_ConeTheta"] 
          && relR.Mod()<settings.map_floatPars["th_ConeR"] 
          && delta_distance<min_delta)
      {
        // std::cout<<"    yyy: nearbyLocalMax["<<il<<"] renewed"<<std::endl;
        shower_candidate_index = il;
        min_delta = delta_distance;
      }
    }
    if (shower_candidate_index<0) break;

    cone_axis.push_back(nearbyLocalMax[shower_candidate_index]);
    nearbyLocalMax.erase(nearbyLocalMax.begin() + shower_candidate_index);
  }
  
  return StatusCode::SUCCESS;
}


bool TrackMatchingAlg::isStopLinking(const std::vector<TVector3>& extrapo_points, 
                   const PandoraPlus::Calo1DCluster* final_cone_hit){
  double slayer = final_cone_hit->getSlayer();
  if(slayer==1){
    TVector3 f_distance = extrapo_points.back() - final_cone_hit->getPos();
    double f_distance_2d = TMath::Sqrt(f_distance.X()*f_distance.X() + f_distance.Y()*f_distance.Y());
    for(int i=0; i<extrapo_points.size(); i++){
      TVector3 distance = extrapo_points[i] - final_cone_hit->getPos();
      double distance_2d = TMath::Sqrt(distance.X()*distance.X() + distance.Y()*distance.Y());
      if (distance_2d < f_distance_2d) return false;
    }
  }
  else{
    TVector3 f_distance = extrapo_points.back() - final_cone_hit->getPos();
    f_distance.RotateZ( TMath::Pi()/4.*(6-final_cone_hit->getTowerID()[0][0]) );
    double f_distance_2d = TMath::Sqrt(f_distance.X()*f_distance.X() + f_distance.Z()*f_distance.Z());
    for(int i=0; i<extrapo_points.size(); i++){
      TVector3 distance = extrapo_points[i] - final_cone_hit->getPos();
      distance.RotateZ( TMath::Pi()/4.*(6-final_cone_hit->getTowerID()[0][0]) );
      double distance_2d = TMath::Sqrt(distance.X()*distance.X() + distance.Z()*distance.Z());
      if (distance_2d < f_distance_2d) return false;
    }
  }
  // std::cout<<"  yyy: calling isStopLinking(): stop!"<<std::endl;
  return true;   
}


TVector2 TrackMatchingAlg::GetProjectedRelR( const PandoraPlus::Calo1DCluster* m_shower1, const PandoraPlus::Calo1DCluster* m_shower2 ){
  TVector2 paxis1, paxis2;
  if(m_shower1->getSlayer()==1){ //For V-bars
    paxis1.Set(m_shower1->getPos().x(), m_shower1->getPos().y());
    paxis2.Set(m_shower2->getPos().x(), m_shower2->getPos().y());
    return paxis2 - paxis1;
  }
  else{  //For U-bars
    if (m_shower1->getTowerID()[0][0] != m_shower2->getTowerID()[0][0])
      std::cout << "warning: In GetProjectedRelR(), modules are different!" << std::endl;
    TVector3 vec = m_shower2->getPos() - m_shower1->getPos();
    vec.RotateZ( TMath::Pi()/4.*(6-m_shower1->getTowerID()[0][0]) );

    TVector2 vec2d(vec.x(), vec.z());

    return vec2d;
  }

  
}


TVector2 TrackMatchingAlg::GetProjectedAxis(const std::vector<TVector3>& extrapo_points, const PandoraPlus::Calo1DCluster* m_shower){
  int min_index=0;
  TVector2 distance(999., 999.);
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
      dist3d.RotateZ( TMath::Pi()/4.*(6-m_shower->getTowerID()[0][0]) );
      TVector2 t_distance(dist3d.x(), dist3d.z());
      if(t_distance.Mod()<distance.Mod()){
        distance = t_distance;
        min_index = i;
      }
    }

    if(min_index < extrapo_points.size()-1){
      double dx = extrapo_points[min_index+1].x() - extrapo_points[min_index].x();
      double dy = extrapo_points[min_index+1].y() - extrapo_points[min_index].y();
      double dz = extrapo_points[min_index+1].z() - extrapo_points[min_index].z();
      TVector3 vec(dx, dy, dz);
      vec.RotateZ( TMath::Pi()/4.*(6-m_shower->getTowerID()[0][0]) );
      TVector2 axis(vec.x(), vec.z());
      return axis;
    }else{
      double dx = extrapo_points[min_index].x() - extrapo_points[min_index-1].x();
      double dy = extrapo_points[min_index].y() - extrapo_points[min_index-1].y();
      double dz = extrapo_points[min_index].z() - extrapo_points[min_index-1].z();
      TVector3 vec(dx, dy, dz);
      vec.RotateZ( TMath::Pi()/4.*(6-m_shower->getTowerID()[0][0]) );
      TVector2 axis(vec.x(), vec.z());
      return axis;
    }
  }

  
}


StatusCode TrackMatchingAlg::CreatConeAxis(PandoraPlusDataCol& m_datacol, PandoraPlus::Track* track, std::vector<PandoraPlus::CaloHalfCluster*>& nearbyHalfClusters, 
                                           std::vector<const PandoraPlus::Calo1DCluster*>& cone_axis){
  // std::cout<<"yyy: Calling CreateConeAxis()"<<std::endl;
  if(nearbyHalfClusters.size()==0 || cone_axis.size()==0) return StatusCode::SUCCESS; 

  for(int ihc=0; ihc<nearbyHalfClusters.size(); ihc++){
    std::vector<const PandoraPlus::Calo1DCluster*> localMaxCol = nearbyHalfClusters[ihc]->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);
    // std::cout<<"    yyy: for nearbyHalfClusters["<<ihc<<"], localMax are:"
    // Track axis candidate.
    // PandoraPlus::CaloHalfCluster* t_track_axis = new PandoraPlus::CaloHalfCluster();
    std::shared_ptr<PandoraPlus::CaloHalfCluster> t_track_axis = std::make_shared<PandoraPlus::CaloHalfCluster>();
    for(int ica=0; ica<cone_axis.size(); ica++){
      if( find(localMaxCol.begin(), localMaxCol.end(), cone_axis[ica]) != localMaxCol.end()){
        t_track_axis->addUnit(cone_axis[ica]);
        // std::cout<<"    add Unit from cone_axis to track_axis: " << cone_axis[ica]->getPos().x() << ", "
        //          << cone_axis[ica]->getPos().y() << ", " << cone_axis[ica]->getPos().z() << std::endl;
      }
    }

    // // If the track does not match the Halfcluster, the track axis candidate will have no 1DCluster
    if(t_track_axis->getCluster().size()==0)
      continue;
    
    t_track_axis->addAssociatedTrack(track);
    t_track_axis->setType(10000); //Track-type axis. 
    
    if(nearbyHalfClusters[ihc]->getSlayer()==1){
      track->addAssociatedHalfClusterV( nearbyHalfClusters[ihc] );
      // std::cout<<"    yyy: track->addAssociatedHalfClusterV"<<std::endl;
    }
    else{
      track->addAssociatedHalfClusterU( nearbyHalfClusters[ihc] );
      // std::cout<<"    yyy: track->addAssociatedHalfClusterU"<<std::endl;
    }
    m_datacol.map_HalfCluster["bkHalfCluster"].push_back(t_track_axis);
    nearbyHalfClusters[ihc]->addHalfCluster(settings.map_stringPars["OutputLongiClusName"], t_track_axis.get());
    // std::cout<<"    yyy: nearbyHalfClusters["<<ihc<<"]->addHalfCluster, bars in the nearbyHalfClusters:"<<std::endl;
    // for(int ii=0; ii<nearbyHalfClusters[ihc]->getCluster().size(); ii++){
    //   std::cout<<"         "<<nearbyHalfClusters[ihc]->getCluster()[ii]->getPos().x()<<", "
    //                         <<nearbyHalfClusters[ihc]->getCluster()[ii]->getPos().y()<<", "
    //                         <<nearbyHalfClusters[ihc]->getCluster()[ii]->getPos().z()<<std::endl;
    // }
    // std::cout<<"          bars in the t_track_axis:"<<std::endl;
    // for(int ii=0; ii<t_track_axis->getCluster().size(); ii++){
    //   std::cout<<"         "<<t_track_axis->getCluster()[ii]->getPos().x()<<", "
    //                         <<t_track_axis->getCluster()[ii]->getPos().y()<<", "
    //                         <<t_track_axis->getCluster()[ii]->getPos().z()<<std::endl;
    // }
    
    
    
  }

  return StatusCode::SUCCESS; 
}

#endif
