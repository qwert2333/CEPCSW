#ifndef TRACKMATCHING_C
#define TRACKMATCHING_C

#include "Algorithm/TrackMatchingAlg.h"


StatusCode TrackMatchingAlg::ReadSettings(PandoraPlus::Settings& m_settings){
  settings = m_settings;
  // ECAL geometry settings
  // Note: Bar half length is also geometry parameter, but obtained from the function GetBarHalfLength()
  if(settings.map_floatPars.find("localmax_area")==settings.map_floatPars.end())
    settings.map_floatPars["localmax_area"] = 10; // unit: mm
  if(settings.map_intPars.find("Nmodule")==settings.map_intPars.end()) 
    settings.map_intPars["Nmodule"] = 8;

  if(settings.map_stringPars.find("OutputLongiClusName")==settings.map_stringPars.end()) 
    settings.map_stringPars["OutputLongiClusName"] = "TrackAxis"; 

  return StatusCode::SUCCESS;
}


StatusCode TrackMatchingAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  m_TrackCol.clear();
  m_HalfClusterV.clear();
  m_HalfClusterU.clear();
  m_trackAxisVCol.clear();
  m_trackAxisUCol.clear();

  m_TrackCol = m_datacol.TrackCol;
  m_HalfClusterU = m_datacol.map_HalfCluster["HalfClusterColU"];
  m_HalfClusterV = m_datacol.map_HalfCluster["HalfClusterColV"];

  return StatusCode::SUCCESS;
}


StatusCode TrackMatchingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){
  std::cout << "---oooOO0OOooo---Excuting TrackMatchingAlg---oooOO0OOooo---"<<std::endl;
  // Associate tracks to HalfClusters.
  // This association is a many-to-many relationship: 
  //    One HalCluster may have multiple tracks; 
  //    One track may pass through multiple HalfClusters.
  for(int ihc=0; ihc<m_HalfClusterV.size(); ihc++){  // loop HalfClusterV
    m_trackAxisVCol.clear();

    // Get local max of the HalfCluster
    std::vector<const PandoraPlus::Calo1DCluster*> localMaxColV = m_HalfClusterV[ihc]->getLocalMaxCol("AllLocalMax");

    for(int itrk=0; itrk<m_TrackCol.size(); itrk++){  // loop tracks
      // Get extrapolated points of the track. These points are sorted by the track
      std::vector<TVector3> extrapo_points;
      GetExtrpoPoints(m_TrackCol[itrk], extrapo_points);

      // Track axis candidate.
      PandoraPlus::CaloHalfCluster* t_track_axis = new PandoraPlus::CaloHalfCluster();
      CreateTrackAxis(extrapo_points, localMaxColV, t_track_axis);

      // If the track do not match the Halfcluster, the track axis candidate will have no 1DCluster
      if(t_track_axis->getCluster().size()==0){
        delete t_track_axis;
        continue;
      }
      
      t_track_axis->addAssociatedTrack(m_TrackCol[itrk]);
      m_trackAxisVCol.push_back(t_track_axis);

    }  // end loop tracks

    m_HalfClusterV[ihc]->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], m_trackAxisVCol);

  }  // end loop HalfClusterV

  for(int ihc=0; ihc<m_HalfClusterU.size(); ihc++){  // loop HalfClusterU
    m_trackAxisUCol.clear();

    // Get local max of the HalfCluster
    std::vector<const PandoraPlus::Calo1DCluster*> localMaxColU = m_HalfClusterU[ihc]->getLocalMaxCol("AllLocalMax");

    for(int itrk=0; itrk<m_TrackCol.size(); itrk++){  // loop tracks
      // Get extrapolated points of the track. These points are sorted by the track
      std::vector<TVector3> extrapo_points;
      GetExtrpoPoints(m_TrackCol[itrk], extrapo_points);

      // Track axis candidate.
      PandoraPlus::CaloHalfCluster* t_track_axis = new PandoraPlus::CaloHalfCluster();
      CreateTrackAxis(extrapo_points, localMaxColU, t_track_axis);

      // If the track do not match the Halfcluster, the track axis candidate will have no 1DCluster
      if(t_track_axis->getCluster().size()==0){
        delete t_track_axis;
        continue;
      }
      
      t_track_axis->addAssociatedTrack(m_TrackCol[itrk]);
      m_trackAxisUCol.push_back(t_track_axis);

    }  // end loop tracks

    m_HalfClusterU[ihc]->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], m_trackAxisUCol);

  }



  // // Program check
  // std::cout << "yyy: check TrackMatchingAlg." << std::endl;
  // for(int ihc=0; ihc<m_HalfClusterV.size(); ihc++){
  //   std::vector<const CaloHalfCluster*> check_track_axis;
  //   check_track_axis = m_HalfClusterV[ihc]->getHalfClusterCol("TrackAxis");
  //   std::cout << setprecision(6);
  //   std::cout << "  HalfclusterV[" << ihc << "], E=" << m_HalfClusterV[ihc]->getEnergy() 
  //             << ", Position=(" << setw(10) << m_HalfClusterV[ihc]->getPos().X() << ", " 
  //                               << setw(10) << m_HalfClusterV[ihc]->getPos().Y() << ", " 
  //                               << setw(10) << m_HalfClusterV[ihc]->getPos().Z() << ") " 
  //             << std::endl;
  //   if(check_track_axis.size() == 0) std::cout << "    No track" << std::endl;
  //   else{
  //     std::cout << "    " << check_track_axis.size() << " track(s) matched to this Halfcluster. The Track axes are:"  << std::endl;
  //     for(int ita=0; ita<check_track_axis.size(); ita++){
  //       std::cout << setw(10) << "x/mm" << setw(10) << "y/mm" << setw(10) << "z/mm" << setw(10) << "E/GeV" << std::endl;
  //       for(int ilm=0; ilm<check_track_axis[ita]->getBars().size(); ilm++){
  //         std::cout << setprecision(6);
  //         std::cout << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getPos().X()
  //                   << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getPos().Y()
  //                   << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getPos().Z()
  //                   << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getEnergy()
  //                   << std::endl;
  //       }
  //       std::cout << "    The track extrapolated points are:" << std::endl;
  //       std::vector<const PandoraPlus::Track*> check_track = check_track_axis[ita]->getAssotiatedTracks();
  //       if(check_track.size()!=1) std::cout << "check_track.size()!=1-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
  //       std::vector<TVector3> check_extrapo_points;
  //       GetExtrpoPoints(check_track[0], check_extrapo_points);
  //       std::cout << setw(10) << "x/mm" << setw(10) << "y/mm" << setw(10) << "z/mm" << std::endl;
  //       for(int ipt=0; ipt<check_extrapo_points.size(); ipt++){
  //         std::cout << setprecision(6);
  //         std::cout << setw(10) << check_extrapo_points[ipt].X()
  //                   << setw(10) << check_extrapo_points[ipt].Y()
  //                   << setw(10) << check_extrapo_points[ipt].Z()
  //                   << std::endl;
  //       }
  //     }
  //   }
  // }
  // for(int ihc=0; ihc<m_HalfClusterU.size(); ihc++){
  //   std::vector<const CaloHalfCluster*> check_track_axis;
  //   check_track_axis = m_HalfClusterU[ihc]->getHalfClusterCol("TrackAxis");
  //   std::cout << setprecision(6);
  //   std::cout << "  HalfclusterU[" << ihc << "], E=" << m_HalfClusterU[ihc]->getEnergy() 
  //             << ", Position=(" << setw(10) << m_HalfClusterU[ihc]->getPos().X() << ", " 
  //                               << setw(10) << m_HalfClusterU[ihc]->getPos().Y() << ", " 
  //                               << setw(10) << m_HalfClusterU[ihc]->getPos().Z() << ") " 
  //             << std::endl;
  //   if(check_track_axis.size() == 0) std::cout << "    No track" << std::endl;
  //   else{
  //     std::cout << "    " << check_track_axis.size() << " track(s) matched to this Halfcluster. The Track axes are:"  << std::endl;
  //     for(int ita=0; ita<check_track_axis.size(); ita++){
  //       std::cout << setw(10) << "x/mm" << setw(10) << "y/mm" << setw(10) << "z/mm" << setw(10) << "E/GeV" << std::endl;
  //       for(int ilm=0; ilm<check_track_axis[ita]->getBars().size(); ilm++){
  //         std::cout << setprecision(6);
  //         std::cout << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getPos().X()
  //                   << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getPos().Y()
  //                   << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getPos().Z()
  //                   << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getEnergy()
  //                   << std::endl;
  //       }
  //       std::cout << "    The track extrapolated points are:" << std::endl;
  //       std::vector<const PandoraPlus::Track*> check_track = check_track_axis[ita]->getAssotiatedTracks();
  //       if(check_track.size()!=1) std::cout << "check_track.size()!=1-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
  //       std::vector<TVector3> check_extrapo_points;
  //       GetExtrpoPoints(check_track[0], check_extrapo_points);
  //       std::cout << setw(10) << "x/mm" << setw(10) << "y/mm" << setw(10) << "z/mm" << std::endl;
  //       for(int ipt=0; ipt<check_extrapo_points.size(); ipt++){
  //         std::cout << setprecision(6);
  //         std::cout << setw(10) << check_extrapo_points[ipt].X()
  //                   << setw(10) << check_extrapo_points[ipt].Y()
  //                   << setw(10) << check_extrapo_points[ipt].Z()
  //                   << std::endl;
  //       }
  //     }
  //   }
  // }

  return StatusCode::SUCCESS;
}


StatusCode TrackMatchingAlg::ClearAlgorithm(){
  m_TrackCol.clear();
  m_HalfClusterV.clear();
  m_HalfClusterU.clear();
  m_trackAxisVCol.clear();
  m_trackAxisUCol.clear();

  return StatusCode::SUCCESS;
}


StatusCode TrackMatchingAlg::GetExtrpoPoints(const PandoraPlus::Track* track, std::vector<TVector3>& extrapo_points){
  for(int its=0; its<track->trackStates_size(); its++){
    if(track->getTrackStates(its).location==PandoraPlus::TrackState::AtOther)
      extrapo_points.push_back(track->getTrackStates(its).referencePoint);
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
         && TMath::Sqrt(distance.X()*distance.X() + distance.Y()*distance.Y()) < settings.map_floatPars["localmax_area"] ) { 
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
         && TMath::Sqrt(distance.X()*distance.X() + distance.Z()*distance.Z()) < settings.map_floatPars["localmax_area"] ) { 
        t_track_axis->addUnit(localMaxCol[ilm]);
      }
      else { continue; }
    }}
  }

std::cout << "end calling CreateTrackAxis()" << std::endl;
  return StatusCode::SUCCESS;
}



#endif