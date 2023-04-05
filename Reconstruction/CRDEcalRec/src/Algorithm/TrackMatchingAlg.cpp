#ifndef TRACKMATCHING_C
#define TRACKMATCHING_C

#include "Algorithm/TrackMatchingAlg.h"


StatusCode TrackMatchingAlg::ReadSettings(PandoraPlus::Settings& m_settings){
  settings = m_settings;
  // ECAL geometry settings
  // Note: Bar half length is also geometry parameter, but obtained from the function GetBarHalfLength()
  if(settings.map_floatPars.find("localmax_area")==settings.map_floatPars.end())
    settings.map_floatPars["localmax_area"] = 10; // unit: mm
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
  m_trackAxisVCol.clear();
  m_trackAxisUCol.clear();

  m_TrackCol = m_datacol.TrackCol;
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
  for(int ihc=0; ihc<p_HalfClusterV->size(); ihc++){  // loop HalfClusterV
    m_trackAxisVCol.clear();

    // Get local max of the HalfCluster
    std::vector<const PandoraPlus::Calo1DCluster*> localMaxColV = p_HalfClusterV->at(ihc)->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

    for(int itrk=0; itrk<m_TrackCol.size(); itrk++){  // loop tracks
      // Get extrapolated points of the track. These points are sorted by the track
      std::vector<TVector3> extrapo_points;
      GetExtrpoPoints(m_TrackCol[itrk], extrapo_points);

      // Track axis candidate.
      PandoraPlus::CaloHalfCluster* t_track_axis = new PandoraPlus::CaloHalfCluster();
      CreateTrackAxis(extrapo_points, localMaxColV, t_track_axis);

      // If the track does not match the Halfcluster, the track axis candidate will have no 1DCluster
      if(t_track_axis->getCluster().size()==0){
        delete t_track_axis;
        continue;
      }
      
      t_track_axis->addAssociatedTrack(m_TrackCol[itrk]);
      t_track_axis->setType(0); //Track-type axis. 
      m_TrackCol[itrk]->addAssociatedHalfClusterV( p_HalfClusterV->at(ihc) );
      m_trackAxisVCol.push_back(t_track_axis);

    }  // end loop tracks

    p_HalfClusterV->at(ihc)->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], m_trackAxisVCol);

  }  // end loop HalfClusterV

  for(int ihc=0; ihc<p_HalfClusterU->size(); ihc++){  // loop HalfClusterU
    m_trackAxisUCol.clear();

    // Get local max of the HalfCluster
    std::vector<const PandoraPlus::Calo1DCluster*> localMaxColU = p_HalfClusterU->at(ihc)->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

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
      t_track_axis->setType(0); //Track-type axis. 
      m_TrackCol[itrk]->addAssociatedHalfClusterU( p_HalfClusterU->at(ihc) );
      m_trackAxisUCol.push_back(t_track_axis);

    }  // end loop tracks

    p_HalfClusterU->at(ihc)->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], m_trackAxisUCol);

  }


printf("TrackMatch: Readin cluster size [%d, %d] \n", p_HalfClusterU->size(), p_HalfClusterV->size());
  //Loop track to check the associated cluster: merge clusters if they are associated to the same track.
  std::vector<PandoraPlus::CaloHalfCluster*> tmp_deleteClus; tmp_deleteClus.clear();
  for(auto &itrk : m_TrackCol){
    std::vector<PandoraPlus::CaloHalfCluster*> m_matchedUCol = itrk->getAssociatedHalfClustersU();
    std::vector<PandoraPlus::CaloHalfCluster*> m_matchedVCol = itrk->getAssociatedHalfClustersV();
//cout<<"  Check in track: matched cluster size ("<<m_matchedUCol.size()<<", "<<m_matchedVCol.size()<<")"<<endl;
//for(int i=0; i<m_matchedUCol.size(); i++) 
//  printf("    %p \n", m_matchedUCol[i]);
//for(int i=0; i<m_matchedVCol.size(); i++) 
//  printf("    %p \n", m_matchedVCol[i]);


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
cout<<"Saved tmp cluster: size = "<<tmp_deleteClus.size()<<endl;
//for(int i=0; i<tmp_deleteClus.size(); i++) printf("  address %p \n", tmp_deleteClus[i] );

  //Check vector: clean the merged clusters
  for(int ihc=0; ihc<p_HalfClusterU->size(); ihc++){
//printf("  Readin HalfClusU address: %p \n", p_HalfClusterU->at(ihc));
    if( find(tmp_deleteClus.begin(), tmp_deleteClus.end(), p_HalfClusterU->at(ihc))!=tmp_deleteClus.end() ){
      p_HalfClusterU->erase(p_HalfClusterU->begin()+ihc);
      ihc--;
    }
  }

  for(int ihc=0; ihc<p_HalfClusterV->size(); ihc++){
//printf("  Readin HalfClusV address: %p \n", p_HalfClusterV->at(ihc));
    if( find(tmp_deleteClus.begin(), tmp_deleteClus.end(), p_HalfClusterV->at(ihc))!=tmp_deleteClus.end() ){
      p_HalfClusterV->erase(p_HalfClusterV->begin()+ihc);
      ihc--;
    }
  }

/*
printf("TrackMatch: After match cluster size [%d, %d] \n", p_HalfClusterU->size(), p_HalfClusterV->size());
for(int i=0; i<p_HalfClusterV->size(); i++){
  cout<<"  In HalfClusterV #"<<i<<": axis size = "<<p_HalfClusterV->at(i)->getHalfClusterCol("TrackAxis").size()<<endl;
  for(int iax=0; iax<p_HalfClusterV->at(i)->getHalfClusterCol("TrackAxis").size(); iax++)
    cout<<"  Axis #"<<iax<<" associated track size: "<<p_HalfClusterV->at(i)->getHalfClusterCol("TrackAxis")[iax]->getAssociatedTracks().size()<<endl;
cout<<endl;
}


   // Program check
   std::cout << "yyy: check TrackMatchingAlg." << std::endl;
   for(int ihc=0; ihc<p_HalfClusterV->size(); ihc++){
     std::vector<const CaloHalfCluster*> check_track_axis;
     check_track_axis = p_HalfClusterV->at(ihc)->getHalfClusterCol("TrackAxis");
     std::cout << setprecision(6);
     std::cout << "  HalfclusterV[" << ihc << "], E=" << p_HalfClusterV->at(ihc)->getEnergy() 
               << ", Position=(" << setw(10) << p_HalfClusterV->at(ihc)->getPos().X() << ", " 
                                 << setw(10) << p_HalfClusterV->at(ihc)->getPos().Y() << ", " 
                                 << setw(10) << p_HalfClusterV->at(ihc)->getPos().Z() << ") " 
               << std::endl;
     if(check_track_axis.size() == 0) std::cout << "    No track" << std::endl;
     else{
       std::cout << "    " << check_track_axis.size() << " track(s) matched to this Halfcluster. The Track axes are:"  << std::endl;
       for(int ita=0; ita<check_track_axis.size(); ita++){
         std::cout << setw(10) << "x/mm" << setw(10) << "y/mm" << setw(10) << "z/mm" << setw(10) << "E/GeV" << std::endl;
         for(int ilm=0; ilm<check_track_axis[ita]->getBars().size(); ilm++){
           std::cout << setprecision(6);
           std::cout << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getPos().X()
                     << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getPos().Y()
                     << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getPos().Z()
                     << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getEnergy()
                     << std::endl;
         }
         std::cout << "    The track extrapolated points are:" << std::endl;
         std::vector<const PandoraPlus::Track*> check_track = check_track_axis[ita]->getAssociatedTracks();
         if(check_track.size()!=1) std::cout << "check_track.size()!=1-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
         std::vector<TVector3> check_extrapo_points;
         GetExtrpoPoints(check_track[0], check_extrapo_points);
         std::cout << setw(10) << "x/mm" << setw(10) << "y/mm" << setw(10) << "z/mm" << std::endl;
         for(int ipt=0; ipt<check_extrapo_points.size(); ipt++){
           std::cout << setprecision(6);
           std::cout << setw(10) << check_extrapo_points[ipt].X()
                     << setw(10) << check_extrapo_points[ipt].Y()
                     << setw(10) << check_extrapo_points[ipt].Z()
                     << std::endl;
         }
       }
     }
   }
   for(int ihc=0; ihc<p_HalfClusterU->size(); ihc++){
     std::vector<const CaloHalfCluster*> check_track_axis;
     check_track_axis = p_HalfClusterU->at(ihc)->getHalfClusterCol("TrackAxis");
     std::cout << setprecision(6);
     std::cout << "  HalfclusterU[" << ihc << "], E=" << p_HalfClusterU->at(ihc)->getEnergy() 
               << ", Position=(" << setw(10) << p_HalfClusterU->at(ihc)->getPos().X() << ", " 
                                 << setw(10) << p_HalfClusterU->at(ihc)->getPos().Y() << ", " 
                                 << setw(10) << p_HalfClusterU->at(ihc)->getPos().Z() << ") " 
               << std::endl;
     if(check_track_axis.size() == 0) std::cout << "    No track" << std::endl;
     else{
       std::cout << "    " << check_track_axis.size() << " track(s) matched to this Halfcluster. The Track axes are:"  << std::endl;
       for(int ita=0; ita<check_track_axis.size(); ita++){
         std::cout << setw(10) << "x/mm" << setw(10) << "y/mm" << setw(10) << "z/mm" << setw(10) << "E/GeV" << std::endl;
         for(int ilm=0; ilm<check_track_axis[ita]->getBars().size(); ilm++){
           std::cout << setprecision(6);
           std::cout << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getPos().X()
                     << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getPos().Y()
                     << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getPos().Z()
                     << setw(10) << check_track_axis[ita]->getCluster()[ilm]->getEnergy()
                     << std::endl;
         }
         std::cout << "    The track extrapolated points are:" << std::endl;
         std::vector<const PandoraPlus::Track*> check_track = check_track_axis[ita]->getAssociatedTracks();
         if(check_track.size()!=1) std::cout << "check_track.size()!=1-=-=-=-=-=-=-=-=-=-=-=-=-" << std::endl;
         std::vector<TVector3> check_extrapo_points;
         GetExtrpoPoints(check_track[0], check_extrapo_points);
         std::cout << setw(10) << "x/mm" << setw(10) << "y/mm" << setw(10) << "z/mm" << std::endl;
         for(int ipt=0; ipt<check_extrapo_points.size(); ipt++){
           std::cout << setprecision(6);
           std::cout << setw(10) << check_extrapo_points[ipt].X()
                     << setw(10) << check_extrapo_points[ipt].Y()
                     << setw(10) << check_extrapo_points[ipt].Z()
                     << std::endl;
         }
       }
     }
   }
*/

  return StatusCode::SUCCESS;
}


StatusCode TrackMatchingAlg::ClearAlgorithm(){
  m_TrackCol.clear();
  p_HalfClusterV = nullptr;
  p_HalfClusterU = nullptr;
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



#endif
