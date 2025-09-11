#ifndef _TRACKCLUSTERCONNECTING_ALG_C
#define _TRACKCLUSTERCONNECTING_ALG_C

#include "Algorithm/TrackClusterConnectingAlg.h"

StatusCode TrackClusterConnectingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_stringPars.find("ReadinECALClusterName")==settings.map_stringPars.end()) settings.map_stringPars["ReadinECALClusterName"] = "EcalCluster";
  if(settings.map_stringPars.find("ReadinHCALClusterName")==settings.map_stringPars.end()) settings.map_stringPars["ReadinHCALClusterName"] = "HcalCluster";
  if(settings.map_floatPars.find("ECALChargedCalib")==settings.map_floatPars.end()) settings.map_floatPars["ECALChargedCalib"] = 1.26;
  if(settings.map_floatPars.find("HCALChargedCalib")==settings.map_floatPars.end()) settings.map_floatPars["HCALChargedCalib"] = 4.0;
  if(settings.map_floatPars.find("ECALNeutralCalib")==settings.map_floatPars.end()) settings.map_floatPars["ECALNeutralCalib"] = 1.0;
  if(settings.map_floatPars.find("HCALNeutralCalib")==settings.map_floatPars.end()) settings.map_floatPars["HCALNeutralCalib"] = 4.0;

  if(settings.map_floatPars.find("th_ChFragEn")==settings.map_floatPars.end()) settings.map_floatPars["th_ChFragEn"] = 2.;
  if(settings.map_floatPars.find("th_ChFragDepth")==settings.map_floatPars.end()) settings.map_floatPars["th_ChFragDepth"] = 100.;
  if(settings.map_floatPars.find("th_ChFragMinR")==settings.map_floatPars.end()) settings.map_floatPars["th_ChFragMinR"] = 200.;
  if(settings.map_floatPars.find("th_HcalMatchingR")==settings.map_floatPars.end()) settings.map_floatPars["th_HcalMatchingR"] = 100.;
  if(settings.map_floatPars.find("th_NeutralHcalMergeTheta")==settings.map_floatPars.end()) settings.map_floatPars["th_NeutralHcalMergeTheta"] = 0.3;
  if(settings.map_floatPars.find("th_NeutralHcalMergeR")==settings.map_floatPars.end()) settings.map_floatPars["th_NeutralHcalMergeR"] = 300.;

  if(settings.map_floatPars.find("th_MIPEnergy")==settings.map_floatPars.end()) settings.map_floatPars["th_MIPEnergy"] = 0.5;
  if(settings.map_floatPars.find("th_AbsorbCone")==settings.map_floatPars.end()) settings.map_floatPars["th_AbsorbCone"] = 0.8;
  if(settings.map_intPars.find("minNhit")==settings.map_intPars.end())                settings.map_intPars["minNhit"] = 5;


  if(settings.map_stringPars.find("OutputMergedECALCluster")==settings.map_stringPars.end()) settings.map_stringPars["OutputMergedECALCluster"] = "TrkMergedECAL";
  if(settings.map_stringPars.find("OutputCombPFO")==settings.map_stringPars.end()) settings.map_stringPars["OutputCombPFO"] = "outputPFO";


  return StatusCode::SUCCESS;
};

StatusCode TrackClusterConnectingAlg::Initialize( CyberDataCol& m_datacol ){
  m_EcalClusters.clear();
  m_HcalClusters.clear();
  m_tracks.clear();
  m_absorbedEcal.clear();
  m_PFObjects.clear();
  m_bkCol.Clear();

  for(int ic=0; ic<m_datacol.map_CaloCluster[settings.map_stringPars["ReadinECALClusterName"]].size(); ic++){
    m_EcalClusters.push_back( m_datacol.map_CaloCluster[settings.map_stringPars["ReadinECALClusterName"]][ic].get() );
  }
  for(int ic=0; ic<m_datacol.map_CaloCluster[settings.map_stringPars["ReadinHCALClusterName"]].size(); ic++){
    m_HcalClusters.push_back( m_datacol.map_CaloCluster[settings.map_stringPars["ReadinHCALClusterName"]][ic].get() );
  }
  for(int itrk=0; itrk<m_datacol.TrackCol.size(); itrk++){
    m_tracks.push_back( m_datacol.TrackCol[itrk].get() );
  }

  m_bkCol.EnergyCorrSvc = m_datacol.EnergyCorrSvc;

cout<<"Readin Track size: "<<m_tracks.size()<<", ECAL cluster size: "<<m_EcalClusters.size()<<", HCAL cluster size "<<m_HcalClusters.size()<<endl;
cout<<"Print track"<<endl;
for(int i=0; i<m_tracks.size(); i++)
  cout<<"Track #"<<i<<": P = "<<m_tracks[i]->getMomentum()<<", Pt = "<<m_tracks[i]->getPt()<<endl;
cout<<"Print all ECAL cluster "<<endl;
for(int ic=0; ic<m_EcalClusters.size(); ic++){
  cout<<"    ECAL Cluster #"<<ic<<": En = "<<m_EcalClusters[ic]->getLongiE()<<", track size "<<m_EcalClusters[ic]->getAssociatedTracks().size();
  if(m_EcalClusters[ic]->getAssociatedTracks().size()>0) cout<<", Leading track P = "<<m_EcalClusters[ic]->getAssociatedTracks()[0]->getMomentum()<<endl;
  else cout<<endl;
}


  return StatusCode::SUCCESS;
};

StatusCode TrackClusterConnectingAlg::RunAlgorithm( CyberDataCol& m_datacol ){
  //Readin: tracks, ECAL clusters and HCAL clusters. 
  //Output: PFObject


  //1. Merge ECAL clusters. 
  //Possible to add some cluster ID and merge functions. 
  m_absorbedEcal.clear();
  EcalChFragAbsorption(m_EcalClusters, m_tracks, m_absorbedEcal);
//cout<<"  TrackClusterConnectingAlg: After ECAL charged fragment absorption: cluster size "<<m_absorbedEcal.size()<<endl;
//cout<<"Print merged ECAL cluster "<<endl;
//for(int ic=0; ic<m_absorbedEcal.size(); ic++){
//  cout<<"    ECAL Cluster #"<<ic<<": En = "<<m_absorbedEcal[ic]->getLongiE()<<", track size "<<m_absorbedEcal[ic]->getAssociatedTracks().size();
//  if(m_absorbedEcal[ic]->getAssociatedTracks().size()>0) cout<<", Leading track P = "<<m_absorbedEcal[ic]->getAssociatedTracks()[0]->getMomentum()<<endl;
//  else cout<<endl;
//}

  //2. Create PFObject with ECAL cluster and track
  std::vector<const Cyber::Calo3DCluster*> tmp_constClus; 
  for(int ic=0; ic<m_absorbedEcal.size(); ic++) tmp_constClus.push_back(m_absorbedEcal[ic].get());
  PFOCreating(tmp_constClus, m_tracks, m_PFObjects);

//cout<<"  TrackClusterConnectingAlg: created PFO: "<<m_PFObjects.size()<<endl;
//for(int i=0; i<m_PFObjects.size(); i++){
//  cout<<"    PFO #"<<i<<": track size "<<m_PFObjects[i]->getTracks().size()<<", leading P "<<m_PFObjects[i]->getTrackMomentum();
//  cout<<", ECAL cluster size "<<m_PFObjects[i]->getECALClusters().size()<<", totE "<<m_PFObjects[i]->getECALClusterEnergy();
//  cout<<", HCAL cluster size "<<m_PFObjects[i]->getHCALClusters().size()<<", totE "<<m_PFObjects[i]->getHCALClusterEnergy()<<endl;
//}

//cout<<"Print all HCAL cluster"<<endl;
//double totE_Hcal = 0;
//for(int i=0; i<m_HcalClusters.size(); i++){
//  printf("  HCAL Cluster #%d: En = %.6f, position (%.3f, %.3f, %.3f), address %p \n", i, 
//            m_HcalClusters[i]->getHitsE(), 
//            m_HcalClusters[i]->getHitCenter().x(), 
//            m_HcalClusters[i]->getHitCenter().y(), 
//            m_HcalClusters[i]->getHitCenter().z(), 
//            m_HcalClusters[i] );
//  totE_Hcal += m_HcalClusters[i]->getHitsE();
//}
//cout<<"Hcal cluster total E "<<totE_Hcal<<endl;

  //3. Add HCAL clusters into the PFObject. 
  std::sort(m_PFObjects.begin(), m_PFObjects.end(), compTrkP);
  HcalExtrapolatingMatch(m_HcalClusters, m_PFObjects);

  //cout<<"  TrackClusterConnectingAlg: PFO size after HCAL matching: "<<m_PFObjects.size()<<endl;
  //for(int i=0; i<m_PFObjects.size(); i++){
  //  cout<<"    PFO #"<<i<<": track size "<<m_PFObjects[i]->getTracks().size()<<", leading P "<<m_PFObjects[i]->getTrackMomentum();
  //  cout<<", ECAL cluster size "<<m_PFObjects[i]->getECALClusters().size()<<", totE "<<m_PFObjects[i]->getECALClusterEnergy();
  //  cout<<", HCAL cluster size "<<m_PFObjects[i]->getHCALClusters().size()<<", totE "<<m_PFObjects[i]->getHCALClusterEnergy()<<endl;
  //  for(int j=0; j<m_PFObjects[i]->getHCALClusters().size(); j++){
  //    printf("      HCAL #%d: energy %.6f, Nhit %d, address %p \n", j, m_PFObjects[i]->getHCALClusters()[j]->getHitsE(), m_PFObjects[i]->getHCALClusters()[j]->getCaloHits().size(), m_PFObjects[i]->getHCALClusters()[j]);
  //  }
  //} 

  m_datacol.map_CaloCluster[ settings.map_stringPars["OutputMergedECALCluster"] ] = m_absorbedEcal;
  m_datacol.map_PFObjects[settings.map_stringPars["OutputCombPFO"]] = m_PFObjects;

  m_datacol.map_CaloCluster["bk3DCluster"].insert( m_datacol.map_CaloCluster["bk3DCluster"].end(), m_bkCol.map_CaloCluster["bk3DCluster"].begin(), m_bkCol.map_CaloCluster["bk3DCluster"].end() );
  m_datacol.map_PFObjects["bkPFO"].insert( m_datacol.map_PFObjects["bkPFO"].end(), m_bkCol.map_PFObjects["bkPFO"].begin(), m_bkCol.map_PFObjects["bkPFO"].end() );
  m_datacol.map_PFObjects["bkPFO"].insert( m_datacol.map_PFObjects["bkPFO"].end(), m_PFObjects.begin(), m_PFObjects.end() );

  return StatusCode::SUCCESS;
};

StatusCode TrackClusterConnectingAlg::ClearAlgorithm(){
  m_EcalClusters.clear();
  m_HcalClusters.clear();
  m_tracks.clear();
  m_absorbedEcal.clear();
  m_PFObjects.clear();
  m_bkCol.Clear();

  return StatusCode::SUCCESS;
};



StatusCode TrackClusterConnectingAlg::PFOCreating( std::vector<const Cyber::Calo3DCluster*>& m_clusters, 
                                                   std::vector<const Cyber::Track*>& m_trks,
                                                   std::vector<std::shared_ptr<Cyber::PFObject>>& m_PFOs ){

//cout<<"  PFOCreating: Track size "<<m_trks.size()<<", Cluster size "<<m_clusters.size()<<endl;
//for(int ic=0; ic<m_clusters.size(); ic++){
//  cout<<"    ECAL Cluster #"<<ic<<": En = "<<m_clusters[ic]->getLongiE()<<", track size "<<m_clusters[ic]->getAssociatedTracks().size();
//  if(m_clusters[ic]->getAssociatedTracks().size()>0) cout<<", Leading track P = "<<m_clusters[ic]->getAssociatedTracks()[0]->getMomentum()<<endl;
//  else cout<<endl;
//}

  std::vector<const Cyber::Track*> m_leftTrks = m_trks; 
  for(int ic=0; ic<m_clusters.size(); ic++){
    std::shared_ptr<Cyber::PFObject> m_newPFO = std::make_shared<Cyber::PFObject>();

    m_newPFO->addECALCluster( m_clusters[ic] );

    std::vector<const Cyber::Track*> m_trkInClus = m_clusters[ic]->getAssociatedTracks();
    if(m_trkInClus.size()!=0){
      m_newPFO->addTrack( m_trkInClus[0] );
      m_newPFO->setPID( m_trkInClus[0]->getPID() );
      auto iter = find(m_leftTrks.begin(), m_leftTrks.end(), m_trkInClus[0]);
      if( iter!=m_leftTrks.end() )
        m_leftTrks.erase(iter);
    }

    m_PFOs.push_back(m_newPFO);
    m_bkCol.map_PFObjects["bkPFO"].push_back(m_newPFO);
  }

  for(int itrk=0; itrk<m_leftTrks.size(); itrk++){
    std::shared_ptr<Cyber::PFObject> m_newPFO = std::make_shared<Cyber::PFObject>();
    m_newPFO->addTrack( m_leftTrks[itrk] );
    m_newPFO->setPID( m_leftTrks[itrk]->getPID() );
    m_PFOs.push_back(m_newPFO);
    m_bkCol.map_PFObjects["bkPFO"].push_back(m_newPFO);
  }

  return StatusCode::SUCCESS;
}


StatusCode TrackClusterConnectingAlg::EcalChFragAbsorption( std::vector<const Cyber::Calo3DCluster*>& m_clusters, 
                                                            std::vector<const Cyber::Track*>& m_trks, 
                                                            std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_newclusCol){
//cout<<"  In EcalChFragAbsorption: Input track size "<<m_trks.size()<<", cluster size "<<m_clusters.size()<<endl;
//for(int ic=0; ic<m_clusters.size(); ic++){
//  cout<<"    ECAL Cluster #"<<ic<<": En = "<<m_clusters[ic]->getLongiE()<<", track size "<<m_clusters[ic]->getAssociatedTracks().size()<<endl;
//}


  //1. Absorb neutral clusters to the nearby tracks
  std::map<const Cyber::Track*, int> m_matchedTrkMap; 
  for(int ic=0; ic<m_clusters.size(); ic++){
    for(int itrk=0; itrk<m_clusters[ic]->getAssociatedTracks().size(); itrk++){
      if( find(m_trks.begin(), m_trks.end(), m_clusters[ic]->getAssociatedTracks()[itrk])!=m_trks.end() )
        m_matchedTrkMap[m_clusters[ic]->getAssociatedTracks()[itrk]] = ic;
    }
  }
//cout<<"  Matched track size: "<<m_matchedTrkMap.size()<<endl;

  for(int ic=0; ic<m_clusters.size(); ic++){
    if(m_clusters[ic]->getAssociatedTracks().size()!=0) continue; 

    double clusDepth = m_clusters[ic]->getDepthToECALSurface();
    double clusEn = m_clusters[ic]->getLongiE();

    double minR2trk = 9999;
    int index = -1;
    for(int itrk=0; itrk<m_trks.size(); itrk++){
      double tmp_minR = GetMinR2Trk( m_clusters[ic], m_trks[itrk]);
      if(tmp_minR<minR2trk){
        minR2trk = tmp_minR;
        index = itrk;
      }
    }
//cout<<"    Clus #"<<ic<<": depth "<<clusDepth<<", En "<<clusEn<<", minR "<<minR2trk<<", index "<<index<<endl;
    //if(index<0){
    //  std::cout<<"ERROR: can not find closest track "<<endl;
    //  continue;
    //}

    if( clusEn<settings.map_floatPars["th_ChFragEn"] && clusDepth>settings.map_floatPars["th_ChFragDepth"] && minR2trk<settings.map_floatPars["th_ChFragMinR"]){
      const Cyber::Track* p_selTrk = m_trks[index]; //Closest track to this cluster. 

      if( m_matchedTrkMap.find(p_selTrk)==m_matchedTrkMap.end() ){ //This track does not match to any existing charged cluster
        std::shared_ptr<Cyber::Calo3DCluster> m_newclus = m_clusters[ic]->Clone();
        m_newclus->addAssociatedTrack(p_selTrk);
        m_newclusCol.push_back( m_newclus );
        m_bkCol.map_CaloCluster["bk3DCluster"].push_back(m_newclus);
      }
      else{ 
        int tmp_index = m_matchedTrkMap[p_selTrk];
        std::shared_ptr<Cyber::Calo3DCluster> m_newclus = m_clusters[tmp_index]->Clone();
        m_newclus->mergeCluster(m_clusters[ic]);
        m_newclusCol.push_back( m_newclus );
        m_matchedTrkMap.erase(p_selTrk);
        m_bkCol.map_CaloCluster["bk3DCluster"].push_back(m_newclus);
      }
    }
    else{
      std::shared_ptr<Cyber::Calo3DCluster> m_newclus = m_clusters[ic]->Clone();
      m_newclusCol.push_back( m_newclus );
      m_bkCol.map_CaloCluster["bk3DCluster"].push_back(m_newclus);
    }

  }

  for(auto iter: m_matchedTrkMap){
    std::shared_ptr<Cyber::Calo3DCluster> m_newclus = m_clusters[iter.second]->Clone();
    m_newclusCol.push_back( m_newclus );
    m_bkCol.map_CaloCluster["bk3DCluster"].push_back(m_newclus);    
  }


  //Merge clusters if linked to the same track
  for(int ic=0; ic<m_newclusCol.size() && m_newclusCol.size()>1; ic++){
    if(m_newclusCol[ic].get()->getAssociatedTracks().size()==0) continue;
    std::vector<const Cyber::Track*> m_trkCol = m_newclusCol[ic].get()->getAssociatedTracks();

    for(int jc=ic+1; jc<m_newclusCol.size(); jc++){
      if(m_newclusCol[jc].get()->getAssociatedTracks().size()==0) continue;

      for(int itrk=0; itrk<m_newclusCol[jc].get()->getAssociatedTracks().size(); itrk++){
        if( find(m_trkCol.begin(), m_trkCol.end(), m_newclusCol[jc].get()->getAssociatedTracks()[itrk])!= m_trkCol.end() ){
          m_newclusCol[ic].get()->mergeCluster( m_newclusCol[jc].get() );
          m_newclusCol.erase(m_newclusCol.begin()+jc);
          jc--;
          if(jc<ic) jc=ic;
        }
        break;
      }
    }
  }
  for(int ic=0; ic<m_newclusCol.size(); ic++) m_newclusCol[ic].get()->getLinkedMCPfromHFCluster("LinkedLongiCluster");



//cout<<"After nearby absorption: Print ECAL cluster "<<endl;
//for(int ic=0; ic<m_newclusCol.size(); ic++){
//  cout<<"    ECAL Cluster #"<<ic<<": En = "<<m_newclusCol[ic]->getLongiE()<<", track size "<<m_newclusCol[ic]->getAssociatedTracks().size();
//  if(m_newclusCol[ic]->getAssociatedTracks().size()>0) cout<<", Leading track P = "<<m_newclusCol[ic]->getAssociatedTracks()[0]->getMomentum()<<endl;
//  else cout<<endl;
//}

  //2. Find the shower vertex, absorb nearby neutral clusters (in a cone) into it. 
  std::vector<std::shared_ptr<Cyber::Calo3DCluster>> m_newChCluster;
  std::vector<std::shared_ptr<Cyber::Calo3DCluster>> m_newNeuCluster;
  for(int icl=0; icl<m_newclusCol.size(); icl++){
    if(m_newclusCol[icl]->getAssociatedTracks().size()==0 ) m_newNeuCluster.push_back(m_newclusCol[icl]);
    else m_newChCluster.push_back(m_newclusCol[icl]);
  }

  for(int icl=0; icl<m_newChCluster.size(); icl++){
    TVector3 cent = m_newChCluster[icl]->getShowerCenter();
    double tmp_Ecl = m_newChCluster[icl]->getLongiE();

    //Veto mip clusters and Eclus>Ptrk
    if(tmp_Ecl<settings.map_floatPars["th_MIPEnergy"]) continue;
    if(tmp_Ecl>m_newChCluster[icl]->getAssociatedTracks()[0]->getMomentum()) continue;

    //Absorb neutral clusters in a cone angle
    for(int jcl=0; jcl<m_newNeuCluster.size(); jcl++){
      //Do not absorb: En, Nhit, start layer, 


      TVector3 vec_pNeu = m_newNeuCluster[jcl]->getShowerCenter();
      if( cent.Angle(vec_pNeu-cent)<settings.map_floatPars["th_AbsorbCone"] ){
        m_newChCluster[icl]->mergeCluster(m_newNeuCluster[jcl].get());
        auto iter = find(m_newclusCol.begin(), m_newclusCol.end(), m_newNeuCluster[jcl]);
        m_newclusCol.erase(iter);
        m_newNeuCluster.erase(m_newNeuCluster.begin()+jcl);
        jcl--;
      }
    }
  }

  //Do cluster energy correction
  for(int ic=0; ic<m_newclusCol.size(); ic++){
    double tmp_clusE = m_newclusCol[ic]->getEnergy()*settings.map_floatPars.at("ECALNeutralCalib");
    TVector3 clus_pos = m_newclusCol[ic]->getShowerCenter();
    double tmp_phi = std::atan2(clus_pos.y(), clus_pos.x())* 180.0 / TMath::Pi(); //TODO: use TVector3 to calculate
    if (tmp_phi < 0) tmp_phi += 360.0;
    double tmp_theta = std::atan2(clus_pos.z(), clus_pos.Perp())* 180.0 / TMath::Pi() + 90;  

    double tmp_scale = m_bkCol.EnergyCorrSvc->energyCorrection(tmp_clusE, tmp_phi, tmp_theta)/tmp_clusE;
    m_newclusCol[ic]->setEnergyScale( tmp_scale );
  }
//for(int ic=0; ic<m_newclusCol.size(); ic++){
//  cout<<"    ECAL Cluster #"<<ic<<": En = "<<m_newclusCol[ic]->getLongiE()<<", track size "<<m_newclusCol[ic]->getAssociatedTracks().size()<<endl;
//}

  return StatusCode::SUCCESS;
};


StatusCode TrackClusterConnectingAlg::HcalExtrapolatingMatch(std::vector<const Cyber::Calo3DCluster*>& m_clusters, std::vector<std::shared_ptr<Cyber::PFObject>>& m_PFOs){

  if(m_clusters.size()<=0){
    std::cout<<"Error in TrackClusterConnectingAlg: No HCAL cluster. Skip this match "<<std::endl;
    return StatusCode::SUCCESS;
  }
//cout<<"  HcalExtrapolatingMatch: Readin cluster size "<<m_clusters.size()<<endl;
//cout<<"  Readin PFO: "<<m_PFOs.size()<<endl;
//for(int i=0; i<m_PFOs.size(); i++){
//  cout<<"    In PFO #"<<i<<": Ntrk "<<m_PFOs[i]->getTracks().size()<<", N Ecal "<<m_PFOs[i]->getECALClusters().size()<<", N Hcal "<<m_PFOs[i]->getHCALClusters().size()<<endl; 
//}


/*  std::vector<const Cyber::CaloHit*> m_isoHits;
  for(int ic=0; ic<m_clusters.size(); ic++){
    if(m_clusters[ic]->getCaloHits().size()<=3){
      for(int ihit=0; ihit<m_clusters[ic]->getCaloHits().size(); ihit++)
        m_isoHits.push_back(m_clusters[ic]->getCaloHits()[ihit]);

      m_clusters.erase(m_clusters.begin()+ic);
    }
  }
cout<<", real cluster size "<<m_clusters.size()<<", isolated hit size "<<m_isoHits.size()<<endl;
*/
  for(int ic=0; ic<m_clusters.size(); ic++){
    std::vector<const Cyber::CaloHit*> hcal_hits = m_clusters[ic]->getCaloHits();
//cout<<"HCAL Cluster #"<<ic<<": Nhit "<<hcal_hits.size()<<", En "<<m_clusters[ic]->getHitsE()<<endl;

    bool isInPfo = false; 
    for(int ipfo=0; ipfo<m_PFOs.size(); ipfo++){
      //Link HCAL cluster to charged PFO
      if(m_PFOs[ipfo]->getTracks().size()!=0){ 
        std::vector<TrackState> trk_points = m_PFOs[ipfo]->getTracks()[0]->getAllTrackStates();

        bool is_candidate = false;
        double minDistance = 99999;
        for(int ihit=0; ihit<hcal_hits.size(); ihit++){
          if(is_candidate) break;
          TVector3 hit_pos = hcal_hits[ihit]->getPosition();
   
          for(int ipts=0; ipts<trk_points.size(); ipts++){
            TVector3 hit_distance = hit_pos - trk_points[ipts].referencePoint;;
            if(minDistance>hit_distance.Mag())  minDistance = hit_distance.Mag();
            if(hit_distance.Mag()<settings.map_floatPars["th_HcalMatchingR"]){
              is_candidate = true;
              break;
            }
          }
        }
//cout<<"  Min distance "<<minDistance<<", is candidate "<<is_candidate<<endl;

        if(is_candidate){
//cout<<"  Pfo #"<<ipfo<<": Ntrk "<<m_PFOs[ipfo]->getTracks().size()<<", leading trk P "<<m_PFOs[ipfo]->getTrackMomentum()<<", trk state size "<<trk_points.size()<<endl;
//cout<<"  Link cluster #"<<ic<<" to pfo #"<<ipfo<<endl;
          m_PFOs[ipfo]->addHCALCluster( m_clusters[ic] );
          isInPfo = true;
          break;
        }        
      }
  
      //Link HCAL cluster to neutral PFO
      else{
        TVector3 pfo_pos(0., 0., 0.);      
        for(int jcl=0; jcl<m_PFOs[ipfo]->getECALClusters().size(); jcl++) pfo_pos += m_PFOs[ipfo]->getECALClusters()[jcl]->getShowerCenter() * m_PFOs[ipfo]->getECALClusters()[jcl]->getLongiE();
        for(int jcl=0; jcl<m_PFOs[ipfo]->getHCALClusters().size(); jcl++) pfo_pos += m_PFOs[ipfo]->getHCALClusters()[jcl]->getHitCenter() * m_PFOs[ipfo]->getHCALClusters()[jcl]->getHitsE();
        pfo_pos = pfo_pos*(1./(m_PFOs[ipfo]->getECALClusterEnergy() + m_PFOs[ipfo]->getHCALClusterEnergy()) );

        TVector3 clus_pos = m_clusters[ic]->getHitCenter();
        if( pfo_pos.Angle(clus_pos)<settings.map_floatPars["th_NeutralHcalMergeTheta"] && (clus_pos-pfo_pos).Mag()<settings.map_floatPars["th_NeutralHcalMergeR"] ){
          m_PFOs[ipfo]->addHCALCluster( m_clusters[ic] );
          isInPfo = true;
          break;
        }
      }

    }//end loop pfos
//if(isInPfo) cout<<"  Merged into PFO: Ptrk = "<<m_PFOs[index_selPfo]->getTrackMomentum()<<endl;

    //If HCAL cluster is not linked to any existing PFO: create a new one. 
    if(!isInPfo && m_clusters[ic]->getCaloHits().size()>=settings.map_intPars["minNhit"]){
//cout<<"  Create a new neutral PFO "<<endl;
      std::shared_ptr<Cyber::PFObject> m_newPFO = std::make_shared<Cyber::PFObject>();
      m_newPFO->addHCALCluster( m_clusters[ic] );     
      m_PFOs.push_back(m_newPFO);
      m_bkCol.map_PFObjects["bkPFO"].push_back(m_newPFO);
    }
    //If only have few hits: merge to closest PFO. 
    else if(!isInPfo){ 

      int index_closest = -1;
      double minDistance = 1e6;
      for(int ipfo=0; ipfo<m_PFOs.size(); ipfo++){
        TVector3 pfo_pos(0., 0., 0.);
        for(int jcl=0; jcl<m_PFOs[ipfo]->getECALClusters().size(); jcl++) pfo_pos += m_PFOs[ipfo]->getECALClusters()[jcl]->getShowerCenter() * m_PFOs[ipfo]->getECALClusters()[jcl]->getLongiE();
        for(int jcl=0; jcl<m_PFOs[ipfo]->getHCALClusters().size(); jcl++) pfo_pos += m_PFOs[ipfo]->getHCALClusters()[jcl]->getHitCenter() * m_PFOs[ipfo]->getHCALClusters()[jcl]->getHitsE();
        pfo_pos = pfo_pos*(1./(m_PFOs[ipfo]->getECALClusterEnergy() + m_PFOs[ipfo]->getHCALClusterEnergy()) );

        double tmp_dis = (pfo_pos-m_clusters[ic]->getHitCenter()).Mag();
        if(tmp_dis<minDistance){
          minDistance = tmp_dis;
          index_closest = ipfo;
        }
      }

      if(index_closest>=0) m_PFOs[index_closest]->addHCALCluster( m_clusters[ic] );
    }
  }

  return StatusCode::SUCCESS;
}


double TrackClusterConnectingAlg::GetMinR2Trk( const Cyber::Calo3DCluster* p_clus, const Cyber::Track* m_trk){
  if(!p_clus || !m_trk) return 99999;

  double minR = 99999;
  int index = -1;
  TVector3 clus_position = p_clus->getShowerCenter();
  std::vector<TrackState> trk_points = m_trk->getAllTrackStates();
  for(int i=0; i<trk_points.size(); i++){
    TVector3 hit_distance = clus_position - trk_points[i].referencePoint;
    if(hit_distance.Mag()<minR) minR = hit_distance.Mag();
  }  

  return minR;
}
#endif

