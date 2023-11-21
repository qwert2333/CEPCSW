#ifndef _PFOCREATING_ALG_C
#define _PFOCREATING_ALG_C

#include "Algorithm/PFOCreatingAlg.h"

StatusCode PFOCreatingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;
  if(settings.map_floatPars.find("delta_phi_cut")==settings.map_floatPars.end())
    settings.map_floatPars["delta_phi_cut"] = 30./180.*TMath::Pi(); 
  if(settings.map_floatPars.find("delta_cosTheta_cut")==settings.map_floatPars.end())
    settings.map_floatPars["delta_cosTheta_cut"] = 0.25; 
  if(settings.map_floatPars.find("extrPoint_HCALHit_distCut")==settings.map_floatPars.end())
    settings.map_floatPars["extrPoint_HCALHit_distCut"] = 100; 
  if(settings.map_floatPars.find("nearby_clus_angleCut")==settings.map_floatPars.end())
    settings.map_floatPars["nearby_clus_angleCut"] = 0.15; 
  return StatusCode::SUCCESS;
};

StatusCode PFOCreatingAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  m_tracks.clear();
  m_ecal_clusters.clear();
  m_hcal_clusters.clear();
  m_pfobjects.clear();

  for(int it=0; it<m_datacol.TrackCol.size(); it++){
    m_tracks.push_back( m_datacol.TrackCol[it].get() );
  }
  for(int ie=0; ie<m_datacol.map_CaloCluster["EcalCluster"].size(); ie++){
    m_ecal_clusters.push_back( m_datacol.map_CaloCluster["EcalCluster"][ie].get() );
  }
  for(int ih=0; ih<m_datacol.map_CaloCluster["HCALCluster"].size(); ih++){
    m_hcal_clusters.push_back( m_datacol.map_CaloCluster["HCALCluster"][ih].get() );
  }

  return StatusCode::SUCCESS;
};

StatusCode PFOCreatingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){
  std::cout << "yyy: Running PFOCreatingAlg" << std::endl;

  if(m_tracks.size()==0 && m_ecal_clusters.size()==0 && m_hcal_clusters.size()==0){
    std::cout << "  yyy: No tracks, no clusters in ECAL and HCAL. End PFOCreatingAlg" << std::endl;
    return StatusCode::SUCCESS;
  }

  // Create PFO with ECAl clusters. If a ECAL cluster is a charged cluster, connect HCAL clusters using extrapolated points
  for(int ie=0; ie<m_ecal_clusters.size(); ie++){
    std::vector<const PandoraPlus::Track*> ecal_cls_track = m_ecal_clusters[ie]->getAssociatedTracks();
    if(ecal_cls_track.size()>1){
      std::cout << "Error! " << ecal_cls_track.size() << " tracks associated to one ECAL cluster!" << std::endl;
      continue;
    }

    // Charged cluster in ECAL (A cluster with a track)
    if(ecal_cls_track.size()==1){  
      // cout << "  yyy: ECAL cluster " << ie << " is a track cluster." << endl;
      // cout << "       (x, y, z, E) = ( " 
      //       << std::fixed << std::setprecision(2) << m_ecal_clusters[ie]->getShowerCenter().x() << ", "
      //       << std::fixed << std::setprecision(2) << m_ecal_clusters[ie]->getShowerCenter().y() << ", "
      //       << std::fixed << std::setprecision(2) << m_ecal_clusters[ie]->getShowerCenter().z() << ", "
      //       << std::fixed << std::setprecision(3) << m_ecal_clusters[ie]->getLongiE() << ")" << endl;

      std::vector<PandoraPlus::Calo3DCluster*> hcal_clus_candidate;
      hcal_clus_candidate.clear();
      GetChargedHCALCandidates(ecal_cls_track[0], m_hcal_clusters, hcal_clus_candidate);

      // cout << "    yyy: Number of connected HCAL clusters: " << hcal_clus_candidate.size() << endl;
      // cout << "         Their (x, y, z, E):" << endl;
      // for(int ihc=0; ihc<hcal_clus_candidate.size(); ihc++){
      //   cout << "           (" 
      //         << std::fixed << std::setprecision(2) << hcal_clus_candidate[ihc]->getHitCenter().x() << ", "
      //         << std::fixed << std::setprecision(2) << hcal_clus_candidate[ihc]->getHitCenter().y() << ", "
      //         << std::fixed << std::setprecision(2) << hcal_clus_candidate[ihc]->getHitCenter().z() << ", "
      //         << std::fixed << std::setprecision(6) << hcal_clus_candidate[ihc]->getHitsE() << ") " << endl; 
      // }

      std::shared_ptr<PandoraPlus::PFObject> tmp_pfo = std::make_shared<PandoraPlus::PFObject>();
      tmp_pfo->addTrack(ecal_cls_track[0]);
      tmp_pfo->addECALCluster(m_ecal_clusters[ie]);
      for(int ic=0; ic<hcal_clus_candidate.size(); ic++){
        tmp_pfo->addHCALCluster(hcal_clus_candidate[ic]);
      }
      m_pfobjects.push_back(tmp_pfo);

      CleanUsedElements(hcal_clus_candidate, m_hcal_clusters);
      CleanUsedElements(ecal_cls_track, m_tracks);
    }

    // Neutral cluster in ECAL (A cluster without track)
    else if(ecal_cls_track.size()==0){
      // cout << "  yyy: ECAL cluster " << ie << " is a neutral cluster." << endl;
      // cout << "       (x, y, z, E) = ( " 
      //       << std::fixed << std::setprecision(2) << m_ecal_clusters[ie]->getShowerCenter().x() << ", "
      //       << std::fixed << std::setprecision(2) << m_ecal_clusters[ie]->getShowerCenter().y() << ", "
      //       << std::fixed << std::setprecision(2) << m_ecal_clusters[ie]->getShowerCenter().z() << ", "
      //       << std::fixed << std::setprecision(3) << m_ecal_clusters[ie]->getLongiE() << ")" << endl;

      // Create PFO with only a ECAL cluster
      std::shared_ptr<PandoraPlus::PFObject> tmp_pfo = std::make_shared<PandoraPlus::PFObject>();
      tmp_pfo->addECALCluster(m_ecal_clusters[ie]);
      m_pfobjects.push_back(tmp_pfo);
      
      
    }
    else{
      cout << "  yyy: Error: Wrong number of tracks" << endl;
    }
  }

  
  // Connect left HCAL clusters to PFO. These PFO are made of ECAL clusters( and tracks) now
  for(int ih=0; ih<m_hcal_clusters.size(); ih++){
    TVector3 hcal_clus_pos = m_hcal_clusters[ih]->getHitCenter();
    double min_angle = 999.0;
    int pfo_index = -1;
    
    for(int ip=0; ip<m_pfobjects.size(); ip++){
      std::vector<const PandoraPlus::Calo3DCluster*> pf_ecal_clus = m_pfobjects[ip].get()->getECALClusters();
      for(int ie=0; ie<pf_ecal_clus.size(); ie++){
        TVector3 ecal_clus_pos = pf_ecal_clus[ie]->getShowerCenter();

        double angle = hcal_clus_pos.Angle(ecal_clus_pos);
        if(angle<settings.map_floatPars["nearby_clus_angleCut"] && angle<min_angle){
          min_angle=angle;
          pfo_index = ip;
        }
      } 
    }

    if(pfo_index>=0){
      m_pfobjects[pfo_index].get()->addHCALCluster(m_hcal_clusters[ih]);
    }
    else{
      std::shared_ptr<PandoraPlus::PFObject> tmp_pfo = std::make_shared<PandoraPlus::PFObject>();
      tmp_pfo->addHCALCluster(m_hcal_clusters[ih]);
      m_pfobjects.push_back(tmp_pfo);
    }

  }

  // Create PFO with only tracks
  for(int it=0; it<m_tracks.size(); it++){
    std::shared_ptr<PandoraPlus::PFObject> tmp_pfo = std::make_shared<PandoraPlus::PFObject>();
    tmp_pfo->addTrack(m_tracks[it]);
    m_pfobjects.push_back(tmp_pfo);
  }
  

  m_datacol.map_PFObjects["outputPFO"] = m_pfobjects;

  

  return StatusCode::SUCCESS;
};

StatusCode PFOCreatingAlg::ClearAlgorithm(){
  m_tracks.clear();
  m_ecal_clusters.clear();
  m_hcal_clusters.clear();
  m_pfobjects.clear();

  return StatusCode::SUCCESS;
}

StatusCode PFOCreatingAlg::GetChargedHCALCandidates(const PandoraPlus::Track* _track,
                                      std::vector<PandoraPlus::Calo3DCluster*>& _hcal_clusters,
                                      std::vector<PandoraPlus::Calo3DCluster*>& _hcal_clus_candidate)
{
  std::vector<TrackState> hcal_trk_states = _track->getTrackStates("Hcal");
  if(hcal_trk_states.size()==0)
    return StatusCode::SUCCESS;


  std::vector<TVector3> extrpolated_points;
  for(int it=0; it<hcal_trk_states.size(); it++){
    extrpolated_points.push_back(hcal_trk_states[it].referencePoint);
  }
  
  for(int ih=0; ih<_hcal_clusters.size(); ih++){
    TVector3 clus_center = _hcal_clusters[ih]->getHitCenter();
    TVector3 distance = clus_center - extrpolated_points[0];
    if(distance.Mag()>1000) continue;  // yyy: harcode for 1000. If the cluster is too far away from the extrpolated points in HCAL, it is obviously not a candidate

    bool is_candidate = false;
    for(int ihit=0; ihit<_hcal_clusters[ih]->getCaloHits().size(); ihit++){
      if(is_candidate) break;
      TVector3 hit_pos = _hcal_clusters[ih]->getCaloHits()[ihit]->getPosition();
      for(int ie=0; ie<extrpolated_points.size(); ie++){
        TVector3 dist = hit_pos - extrpolated_points[ie];
        if(dist.Mag()<settings.map_floatPars["extrPoint_HCALHit_distCut"]){
          is_candidate = true;
          break;
        }
      }
    }

    if(is_candidate){
      _hcal_clus_candidate.push_back(_hcal_clusters[ih]);
    }
  }
} 

StatusCode PFOCreatingAlg::GetNearbyHCALCandidates( PandoraPlus::Calo3DCluster* _ecal_cluster,
                                  std::vector<PandoraPlus::Calo3DCluster*>& _hcal_clusters,
                                  std::vector<PandoraPlus::Calo3DCluster*>& _hcal_clus_candidate)
{
  TVector3 ecal_pos = _ecal_cluster->getShowerCenter();
  for(int ih=0; ih<_hcal_clusters.size(); ih++){
    TVector3 hcal_pos = _hcal_clusters[ih]->getHitCenter();

    double angle = ecal_pos.Angle(hcal_pos);
    if(angle<settings.map_floatPars["nearby_clus_angleCut"]){
      _hcal_clus_candidate.push_back(_hcal_clusters[ih]);
    }
  }

  return StatusCode::SUCCESS;
}

bool PFOCreatingAlg::isReachOuterMostECAL(PandoraPlus::Calo3DCluster* _ecal_cluster)
{
  // A neutral cluster may deposits energy into HCAL if its cluster in ECAL reach the outermost boundary.
  //   In this case, it hits the last layer of ECAL or hit the boundary of different modules.
  const CaloHalfCluster* p_HFClusterU = _ecal_cluster->getHalfClusterUCol("LinkedLongiCluster")[0];
  const CaloHalfCluster* p_HFClusterV = _ecal_cluster->getHalfClusterVCol("LinkedLongiCluster")[0];
  int endLayer = max(p_HFClusterU->getEndDlayer(), p_HFClusterV->getEndDlayer());
  if (endLayer>13){  // yyy: hardcode! number of layer in ECAL may be different from 14 for other design
    cout << "    yyy: The ECAL cluster reach the outermost ECAL. endLayer=" <<  endLayer << endl;
    return true;
  }

  std::vector< std::vector<int> > cluster_towers = _ecal_cluster->getTowerID();
  set<int> towerID;
  for(int it=0; it<cluster_towers.size(); it++){
    towerID.insert(cluster_towers[it][0]);
  }
  if(towerID.size()>1){
    return true;

    cout << "    yyy: The ECAL cluster deposits in over one module. towerID = (";
    for (int i : towerID) {
      cout << i << ", ";
    }cout << ")" << endl;
  }

  cout << "    yyy: The ECAL cluster does not reach the outermost ECAL. EndLayer=" <<  endLayer 
        << ", towerID = (";
  for (int i : towerID) {
    cout << i << ", ";
  }cout << ")" << endl;

  return false;
}

template<typename T1, typename T2>
StatusCode PFOCreatingAlg::CleanUsedElements(std::vector<T1>& _used_elements,
                                                    std::vector<T2>& _left_elements)
{
  for(int i=0; i<_used_elements.size(); i++){
    auto it = std::find(_left_elements.begin(), _left_elements.end(), _used_elements[i]);
    if (it != _left_elements.end()) {
      _left_elements.erase(it);
    }
  }

  return StatusCode::SUCCESS;
}
template<typename T1, typename T2> StatusCode CleanUsedElement(T1 _used_element,
                                                    std::vector<T2>& _left_elements)
{
  auto it = std::find(_left_elements.begin(), _left_elements.end(), _used_element);

  if (it != _left_elements.end()) {
      _left_elements.erase(it);
  }

  return StatusCode::SUCCESS;
}

StatusCode PFOCreatingAlg::CreateLeftPFO(std::vector<PandoraPlus::Track*>& _tracks,
                            std::vector<PandoraPlus::Calo3DCluster*>& _hcal_clusters,
                            std::vector<std::shared_ptr<PandoraPlus::PFObject>>& _pfobjects)
{
  for(int it=0; it<_tracks.size(); it++){
    std::shared_ptr<PandoraPlus::PFObject> tmp_pfo = std::make_shared<PandoraPlus::PFObject>();
    tmp_pfo->addTrack(_tracks[it]);
    _pfobjects.push_back(tmp_pfo);
  }

  for(int ih=0; ih<_hcal_clusters.size(); ih++){
    std::shared_ptr<PandoraPlus::PFObject> tmp_pfo = std::make_shared<PandoraPlus::PFObject>();
    tmp_pfo->addHCALCluster(_hcal_clusters[ih]);
    _pfobjects.push_back(tmp_pfo);
  }
  return StatusCode::SUCCESS;
} 

#endif