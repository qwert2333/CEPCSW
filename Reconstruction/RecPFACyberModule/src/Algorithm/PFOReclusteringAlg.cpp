#ifndef _PFORECLUSTERING_ALG_C
#define _PFORECLUSTERING_ALG_C

#include "Algorithm/PFOReclusteringAlg.h"

StatusCode PFOReclusteringAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_stringPars.find("ReadinPFOName")==settings.map_stringPars.end()) settings.map_stringPars["ReadinPFOName"] = "outputPFO";
  if(settings.map_floatPars.find("ECALChargedCalib")==settings.map_floatPars.end()) settings.map_floatPars["ECALChargedCalib"] = 1.26;
  if(settings.map_floatPars.find("ECALChargedEMCalib")==settings.map_floatPars.end()) settings.map_floatPars["ECALChargedEMCalib"] = 1.0;
  if(settings.map_floatPars.find("HCALChargedCalib")==settings.map_floatPars.end()) settings.map_floatPars["HCALChargedCalib"] = 4.0;
  if(settings.map_floatPars.find("ECALNeutralCalib")==settings.map_floatPars.end()) settings.map_floatPars["ECALNeutralCalib"] = 1.0;
  if(settings.map_floatPars.find("HCALNeutralCalib")==settings.map_floatPars.end()) settings.map_floatPars["HCALNeutralCalib"] = 4.0;
  
  if(settings.map_floatPars.find("ElecEoverP")==settings.map_floatPars.end()) settings.map_floatPars["ElecEoverP"] = 0.8;
  if(settings.map_floatPars.find("EnergyRes")==settings.map_floatPars.end()) settings.map_floatPars["EnergyRes"] = 0.4;
  if(settings.map_floatPars.find("SplitSigma")==settings.map_floatPars.end()) settings.map_floatPars["SplitSigma"] = 0.2;
  if(settings.map_floatPars.find("NeutralMergeSigma")==settings.map_floatPars.end()) settings.map_floatPars["NeutralMergeSigma"] = 0.;
  if(settings.map_floatPars.find("VirtualMergeSigma")==settings.map_floatPars.end()) settings.map_floatPars["VirtualMergeSigma"] = 0.4;
  if(settings.map_floatPars.find("MinAngleForNeuMerge")==settings.map_floatPars.end()) settings.map_floatPars["MinAngleForNeuMerge"] = 0.20;
  if(settings.map_floatPars.find("MinAngleForVirMerge")==settings.map_floatPars.end()) settings.map_floatPars["MinAngleForVirMerge"] = 0.24;


  return StatusCode::SUCCESS;
};

StatusCode PFOReclusteringAlg::Initialize( CyberDataCol& m_datacol ){
  p_PFObjects = nullptr;

  p_PFObjects = &(m_datacol.map_PFObjects[settings.map_stringPars["ReadinPFOName"]]);
  return StatusCode::SUCCESS;
};

StatusCode PFOReclusteringAlg::RunAlgorithm( CyberDataCol& m_datacol ){

  std::vector< std::shared_ptr<Cyber::PFObject> > m_chargedPFOs; 
  std::vector< std::shared_ptr<Cyber::PFObject> > m_neutralPFOs; 

  for(int ipfo=0; ipfo<p_PFObjects->size(); ipfo++){
    if(p_PFObjects->at(ipfo)->getTracks().size()==0) m_neutralPFOs.push_back( p_PFObjects->at(ipfo) );
    else m_chargedPFOs.push_back( p_PFObjects->at(ipfo) );
  }

  std::sort(m_chargedPFOs.begin(), m_chargedPFOs.end(), compTrkP);

  //If P_trk < E_cluster, create a virtual neutral PFO. 
  ReCluster_SplitFromChg(m_chargedPFOs, m_neutralPFOs);

  // double totE_Ecal = 0;
  // double totE_Hcal = 0;
  // cout<<"PFO after splitting: "<<p_PFObjects->size()<<", charged "<<m_chargedPFOs.size()<<", neutral "<<m_neutralPFOs.size()<<endl;
  // for(int i=0; i<m_neutralPFOs.size(); i++){
  //   cout<<"    PFO #"<<i<<": track size "<<m_neutralPFOs[i]->getTracks().size()<<", leading P "<<m_neutralPFOs[i]->getTrackMomentum();
  //   cout<<", ECAL cluster size "<<m_neutralPFOs[i]->getECALClusters().size()<<", totE "<<m_neutralPFOs[i]->getECALClusterEnergy()*settings.map_floatPars["ECALNeutralCalib"];
  //   cout<<", HCAL cluster size "<<m_neutralPFOs[i]->getHCALClusters().size()<<", totE "<<m_neutralPFOs[i]->getHCALClusterEnergy()*settings.map_floatPars["HCALNeutralCalib"]<<endl;
  //   totE_Ecal += m_neutralPFOs[i]->getECALClusterEnergy()*settings.map_floatPars["ECALNeutralCalib"];
  //   totE_Hcal += m_neutralPFOs[i]->getHCALClusterEnergy()*settings.map_floatPars["HCALNeutralCalib"];
  // }
  // cout<<"-----Neutral cluster Ecal total energy: "<<totE_Ecal<<", Hcal total energy: "<<totE_Hcal<<endl;
  // totE_Ecal = 0;
  // totE_Hcal = 0;
  // for(int i=0; i<m_chargedPFOs.size(); i++){
  //   cout<<"    PFO #"<<i<<": track size "<<m_chargedPFOs[i]->getTracks().size()<<", leading P "<<m_chargedPFOs[i]->getTrackMomentum();
  //   cout<<", ECAL cluster size "<<m_chargedPFOs[i]->getECALClusters().size()<<", totE "<<m_chargedPFOs[i]->getECALClusterEnergy()*settings.map_floatPars["ECALChargedCalib"];
  //   cout<<", HCAL cluster size "<<m_chargedPFOs[i]->getHCALClusters().size()<<", totE "<<m_chargedPFOs[i]->getHCALClusterEnergy()*settings.map_floatPars["HCALChargedCalib"]<<endl;
  //   totE_Ecal += m_chargedPFOs[i]->getECALClusterEnergy()*settings.map_floatPars["ECALChargedCalib"];
  //   totE_Hcal += m_chargedPFOs[i]->getHCALClusterEnergy()*settings.map_floatPars["HCALChargedCalib"];
  // }
  // cout<<"-----Charged cluster Ecal total energy: "<<totE_Ecal<<", Hcal total energy: "<<totE_Hcal<<endl;

  //If P_trk > E_cluster, merge nearby neutral PFO into the charged. 
  ReCluster_MergeToChg(m_chargedPFOs, m_neutralPFOs);

  m_datacol.map_CaloHit["bkHit"].insert( m_datacol.map_CaloHit["bkHit"].end(), m_bkCol.map_CaloHit["bkHit"].begin(), m_bkCol.map_CaloHit["bkHit"].end() );
  m_datacol.map_BarCol["bkBar"].insert( m_datacol.map_BarCol["bkBar"].end(), m_bkCol.map_BarCol["bkBar"].begin(), m_bkCol.map_BarCol["bkBar"].end() );
  m_datacol.map_1DCluster["bk1DCluster"].insert( m_datacol.map_1DCluster["bk1DCluster"].end(), m_bkCol.map_1DCluster["bk1DCluster"].begin(), m_bkCol.map_1DCluster["bk1DCluster"].end() );
  m_datacol.map_2DCluster["bk2DCluster"].insert( m_datacol.map_2DCluster["bk2DCluster"].end(), m_bkCol.map_2DCluster["bk2DCluster"].begin(), m_bkCol.map_2DCluster["bk2DCluster"].end() );
  m_datacol.map_HalfCluster["bkHalfCluster"].insert( m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_bkCol.map_HalfCluster["bkHalfCluster"].begin(), m_bkCol.map_HalfCluster["bkHalfCluster"].end() );
  m_datacol.map_CaloCluster["bk3DCluster"].insert( m_datacol.map_CaloCluster["bk3DCluster"].end(), m_bkCol.map_CaloCluster["bk3DCluster"].begin(), m_bkCol.map_CaloCluster["bk3DCluster"].end() );
  m_datacol.map_PFObjects["bkPFO"].insert( m_datacol.map_PFObjects["bkPFO"].end(), m_bkCol.map_PFObjects["bkPFO"].begin(), m_bkCol.map_PFObjects["bkPFO"].end() );


  //double totE_Ecal = 0;
  //double totE_Hcal = 0;
  //cout<<"After merge all virtual to Ch: charged "<<m_chargedPFOs.size()<<", neutral "<<m_neutralPFOs.size()<<", total "<<p_PFObjects->size()<<endl;
  //for(int i=0; i<m_neutralPFOs.size(); i++){
  //  cout<<"    PFO #"<<i<<": track size "<<m_neutralPFOs[i]->getTracks().size()<<", leading P "<<m_neutralPFOs[i]->getTrackMomentum();
  //  cout<<", ECAL cluster size "<<m_neutralPFOs[i]->getECALClusters().size()<<", totE "<<m_neutralPFOs[i]->getECALClusterEnergy()*settings.map_floatPars["ECALNeutralCalib"];
  //  cout<<", HCAL cluster size "<<m_neutralPFOs[i]->getHCALClusters().size()<<", totE "<<m_neutralPFOs[i]->getHCALClusterEnergy()*settings.map_floatPars["HCALNeutralCalib"]<<endl;
  //  totE_Ecal += m_neutralPFOs[i]->getECALClusterEnergy()*settings.map_floatPars["ECALNeutralCalib"];
  //  totE_Hcal += m_neutralPFOs[i]->getHCALClusterEnergy()*settings.map_floatPars["HCALNeutralCalib"];
  //}
  //cout<<"-----Neutral cluster Ecal total energy: "<<totE_Ecal<<", Hcal total energy: "<<totE_Hcal<<endl;
  //totE_Ecal = 0;
  //totE_Hcal = 0;
  //for(int i=0; i<m_chargedPFOs.size(); i++){
  //  cout<<"    PFO #"<<i<<": track size "<<m_chargedPFOs[i]->getTracks().size()<<", leading P "<<m_chargedPFOs[i]->getTrackMomentum();
  //  cout<<", ECAL cluster size "<<m_chargedPFOs[i]->getECALClusters().size()<<", totE "<<m_chargedPFOs[i]->getECALClusterEnergy()*settings.map_floatPars["ECALChargedCalib"];
  //  cout<<", HCAL cluster size "<<m_chargedPFOs[i]->getHCALClusters().size()<<", totE "<<m_chargedPFOs[i]->getHCALClusterEnergy()*settings.map_floatPars["HCALChargedCalib"]<<endl;
  //  totE_Ecal += m_chargedPFOs[i]->getECALClusterEnergy()*settings.map_floatPars["ECALChargedCalib"];
  //  totE_Hcal += m_chargedPFOs[i]->getHCALClusterEnergy()*settings.map_floatPars["HCALChargedCalib"];
  //}
  //cout<<"-----Charged cluster Ecal total energy: "<<totE_Ecal<<", Hcal total energy: "<<totE_Hcal<<endl;


  m_chargedPFOs.clear();
  m_neutralPFOs.clear();
  return StatusCode::SUCCESS;
};

StatusCode PFOReclusteringAlg::ClearAlgorithm(){
  p_PFObjects = nullptr;
  m_bkCol.Clear();

  return StatusCode::SUCCESS;
};


StatusCode PFOReclusteringAlg::ReCluster_MergeToChg(std::vector< std::shared_ptr<Cyber::PFObject> >& m_chargedPFOs, 
                                                    std::vector< std::shared_ptr<Cyber::PFObject> >& m_neutralPFOs ){

  //Merge real neutral PFOs
  for(int ic=0; ic<m_chargedPFOs.size(); ic++){
    if(m_chargedPFOs[ic]->getECALClusters().size()==0 && m_chargedPFOs[ic]->getHCALClusters().size()==0) continue;

    double track_energy = m_chargedPFOs[ic]->getTrackMomentum();
    double ECAL_energy = settings.map_floatPars["ECALChargedCalib"]*m_chargedPFOs[ic]->getECALClusterEnergy();
    double HCAL_energy = settings.map_floatPars["HCALChargedCalib"]*m_chargedPFOs[ic]->getHCALClusterEnergy();
    if(track_energy<0 || ECAL_energy<0 || HCAL_energy<0){
      std::cout<<"ERROR: Charged PFO info break. Ptrk "<<track_energy<<", E_ecal "<<ECAL_energy<<", E_hcal "<<HCAL_energy<<endl;
      continue;
    }

    double delta_energy = ECAL_energy + HCAL_energy - track_energy;
    double sigmaE = settings.map_floatPars["EnergyRes"] * sqrt(ECAL_energy + HCAL_energy);

//cout<<"  ReCluster_MergeToChg: In ChPFO #"<<ic;
//cout<<": ECAL cluster size "<<m_chargedPFOs[ic]->getECALClusters().size()<<", HCAL cluster size "<<m_chargedPFOs[ic]->getHCALClusters().size();
//cout<<", Ptrk = "<<track_energy<<", Eecal = "<<ECAL_energy<<", Ehcal = "<<HCAL_energy<<", deltaE = "<<delta_energy<<", sigmaE = "<<sigmaE<<endl;

    if(delta_energy - settings.map_floatPars["NeutralMergeSigma"]*sigmaE >= -1e-9) continue;
//cout<<"    Do charged PFO merge. "<<endl;

    // All HCAL clusters in the neutral PFO
    std::vector<const Cyber::Calo3DCluster*> all_neutral_HCAL_clus;
    for(int ip=0; ip<m_neutralPFOs.size(); ip++){
      std::vector<const Cyber::Calo3DCluster*> tmp_HCAL_clus = m_neutralPFOs[ip]->getHCALClusters();
      all_neutral_HCAL_clus.insert(all_neutral_HCAL_clus.end(), tmp_HCAL_clus.begin(), tmp_HCAL_clus.end());
    }
//cout<<"      Neutral HCAL cluster size "<<all_neutral_HCAL_clus.size()<<endl;
//cout<<"      Print all HCAL clusters "<<endl;
//for(int acl=0; acl<all_neutral_HCAL_clus.size(); acl++){
//printf("        Cluster #%d: energy %.8f, position (%.3f, %.3f, %.3f) \n", acl, settings.map_floatPars["HCALCalib"]*all_neutral_HCAL_clus[acl]->getHitsE(), all_neutral_HCAL_clus[acl]->getHitCenter().x(), all_neutral_HCAL_clus[acl]->getHitCenter().y(), all_neutral_HCAL_clus[acl]->getHitCenter().z() );
//}

    TVector3 trackclus_pos(0, 0, 0);
    if(m_chargedPFOs[ic]->getECALClusters().size()>0){
      trackclus_pos = m_chargedPFOs[ic]->getECALClusters()[0]->getShowerCenter();
    }
    else{
      for(int jc=0; jc<m_chargedPFOs[ic]->getHCALClusters().size(); jc++)
        trackclus_pos += m_chargedPFOs[ic]->getHCALClusters()[jc]->getHitCenter() * settings.map_floatPars["HCALChargedCalib"]*m_chargedPFOs[ic]->getHCALClusters()[jc]->getHitsE();

      trackclus_pos = trackclus_pos*(1./HCAL_energy);
    }
//printf("    charged cluster trk pos: [%.3f, %.3f, %.3f] \n", trackclus_pos.x(), trackclus_pos.y(), trackclus_pos.z());

    int loop_count = 0;
    std::vector<const Cyber::Calo3DCluster*> skip_clus;
    while(delta_energy < 0){
      loop_count++;
      if(loop_count>all_neutral_HCAL_clus.size()+10){
//cout<<"    --- Out of loop number: "<<loop_count<<endl;
        break;
      }
      double min_angle = 999.0;
      int clus_index = -1;

      for(int in=0; in<all_neutral_HCAL_clus.size(); in++){
        if(find(skip_clus.begin(), skip_clus.end(), all_neutral_HCAL_clus[in]) != skip_clus.end()) continue;

        TVector3 neutral_clus_pos = all_neutral_HCAL_clus[in]->getHitCenter();
        double pfo_angle = trackclus_pos.Angle(neutral_clus_pos);
        if (pfo_angle<settings.map_floatPars["MinAngleForNeuMerge"] && pfo_angle<min_angle){
          min_angle=pfo_angle;
          clus_index = in;
        }
//printf("      Check neutral HCAL cluster #%d: energy %.8f, position (%.3f, %.3f, %.3f), angle with ch: %.4f, current closest cluster %d, minAngle %.4f \n", in, settings.map_floatPars["HCALCalib"]*all_neutral_HCAL_clus[in]->getHitsE(), all_neutral_HCAL_clus[in]->getHitCenter().x(), all_neutral_HCAL_clus[in]->getHitCenter().y(), all_neutral_HCAL_clus[in]->getHitCenter().z(), pfo_angle, clus_index, min_angle );
      }

//cout<<"    In Loop "<<loop_count<<": current deltaE = "<<delta_energy<<", closest virtual cluster index "<<clus_index<<", min Angle "<<min_angle<<endl;
      if(clus_index<0) break;  // No neutral Hcal cluster to be merged

      double tmp_delta_E = delta_energy + settings.map_floatPars["HCALNeutralCalib"]*all_neutral_HCAL_clus[clus_index]->getHitsE();
      if (TMath::Abs(tmp_delta_E) >= TMath::Abs(delta_energy)){
        skip_clus.push_back(all_neutral_HCAL_clus[clus_index]);
        continue;  // No need to merge this HCAL cluster
      }

      // Update delta_energy
      delta_energy = tmp_delta_E;
      // Add this HCAL cluster to charged PFO
//cout<<"    Add neutral HCAL cluster #"<<clus_index<<" into charged PFO #"<<ic<<". Cluster E "<<settings.map_floatPars["HCALCalib"]*all_neutral_HCAL_clus[clus_index]->getHitsE()<<", deltaE "<<delta_energy<<endl;
     m_chargedPFOs[ic]->addHCALCluster(all_neutral_HCAL_clus[clus_index]);

      // Remove this HCAL cluster from neutral PFO
      bool is_found = false;
      for(int in=0; in<m_neutralPFOs.size(); in++){
        std::vector<const Cyber::Calo3DCluster*> neutral_cluster = m_neutralPFOs[in]->getHCALClusters();
        int tmp_index=-1;
        for(int ii=0; ii<neutral_cluster.size(); ii++){
          if (all_neutral_HCAL_clus[clus_index]==neutral_cluster[ii]){
            tmp_index = ii;
            break;
          }
        }
        if (tmp_index==-1) continue;
//cout<<"  Remove a neutral cluster: En "<<settings.map_floatPars["HCALCalib"]*neutral_cluster[tmp_index]->getHitsE()<<endl;

        neutral_cluster.erase(neutral_cluster.begin()+tmp_index);
        m_neutralPFOs[in]->setHCALCluster(neutral_cluster);

        if(m_neutralPFOs[in]->getTracks().size() + m_neutralPFOs[in]->getECALClusters().size() + m_neutralPFOs[in]->getHCALClusters().size()==0){
//cout<<"  Remove a neutral PFO: ECAL En "<<m_neutralPFOs[in]->getECALClusterEnergy()<<", HCAL En "<<m_neutralPFOs[in]->getHCALClusterEnergy()<<endl;
          auto iter = find(p_PFObjects->begin(), p_PFObjects->end(), m_neutralPFOs[in]);
          if(iter==p_PFObjects->end()){
            std::cout<<"ERROR: can not find this neutral PFO in p_PFObjects. "<<std::endl;
          }
          else{
            m_neutralPFOs.erase(m_neutralPFOs.begin()+in);
            p_PFObjects->erase(iter);
          }
        }
        is_found = true;
        break;
      }
      if(!is_found){
        cout << "Error! Can not find the HCAL cluster in neutral PFO" << endl;
      }

      // Remove this HCAL cluster from all_neutral_HCAL_clus
      all_neutral_HCAL_clus.erase(all_neutral_HCAL_clus.begin()+clus_index);
    }
  }


  // double totE_Ecal = 0;
  // double totE_Hcal = 0;
  // cout<<"After merge real neu to Ch: charged "<<m_chargedPFOs.size()<<", neutral "<<m_neutralPFOs.size()<<", total "<<p_PFObjects->size()<<endl;
  // for(int i=0; i<m_neutralPFOs.size(); i++){
  //   cout<<"    PFO #"<<i<<": track size "<<m_neutralPFOs[i]->getTracks().size()<<", leading P "<<m_neutralPFOs[i]->getTrackMomentum();
  //   cout<<", ECAL cluster size "<<m_neutralPFOs[i]->getECALClusters().size()<<", totE "<<m_neutralPFOs[i]->getECALClusterEnergy()*settings.map_floatPars["ECALNeutralCalib"];
  //   cout<<", HCAL cluster size "<<m_neutralPFOs[i]->getHCALClusters().size()<<", totE "<<m_neutralPFOs[i]->getHCALClusterEnergy()*settings.map_floatPars["HCALNeutralCalib"]<<endl;
  //   totE_Ecal += m_neutralPFOs[i]->getECALClusterEnergy()*settings.map_floatPars["ECALNeutralCalib"];
  //   totE_Hcal += m_neutralPFOs[i]->getHCALClusterEnergy()*settings.map_floatPars["HCALNeutralCalib"];
  // }
  // cout<<"-----Neutral cluster Ecal total energy: "<<totE_Ecal<<", Hcal total energy: "<<totE_Hcal<<endl;
  // totE_Ecal = 0;
  // totE_Hcal = 0;
  // for(int i=0; i<m_chargedPFOs.size(); i++){
  //   cout<<"    PFO #"<<i<<": track size "<<m_chargedPFOs[i]->getTracks().size()<<", leading P "<<m_chargedPFOs[i]->getTrackMomentum();
  //   cout<<", ECAL cluster size "<<m_chargedPFOs[i]->getECALClusters().size()<<", totE "<<m_chargedPFOs[i]->getECALClusterEnergy()*settings.map_floatPars["ECALChargedCalib"];
  //   cout<<", HCAL cluster size "<<m_chargedPFOs[i]->getHCALClusters().size()<<", totE "<<m_chargedPFOs[i]->getHCALClusterEnergy()*settings.map_floatPars["HCALChargedCalib"]<<endl;
  //   totE_Ecal += m_chargedPFOs[i]->getECALClusterEnergy()*settings.map_floatPars["ECALChargedCalib"];
  //   totE_Hcal += m_chargedPFOs[i]->getHCALClusterEnergy()*settings.map_floatPars["HCALChargedCalib"];
  // }
  // cout<<"-----Charged cluster Ecal total energy: "<<totE_Ecal<<", Hcal total energy: "<<totE_Hcal<<endl;



  //Merge virtual neutral PFOs created from splitting. 
  for(int ic=0; ic<m_chargedPFOs.size(); ic++){
    if(m_chargedPFOs[ic]->getTracks().size()==0) continue;
    if( m_chargedPFOs[ic]->getECALClusters().size()==0 && 
        m_chargedPFOs[ic]->getHCALClusters().size()==0 && 
        m_chargedPFOs[ic]->getTracks()[0]->getTrackStates("Hcal").size()==0) continue;

    double track_energy = m_chargedPFOs[ic]->getTrackMomentum();
    double ECAL_energy = settings.map_floatPars["ECALChargedCalib"]*m_chargedPFOs[ic]->getECALClusterEnergy();
    double HCAL_energy = settings.map_floatPars["HCALChargedCalib"]*m_chargedPFOs[ic]->getHCALClusterEnergy();
    if(track_energy<0 || ECAL_energy<0 || HCAL_energy<0){
      std::cout<<"ERROR: Charged PFO info break. Ptrk "<<track_energy<<", E_ecal "<<ECAL_energy<<", E_hcal "<<HCAL_energy<<endl;
      continue;
    }

    double sigmaE = settings.map_floatPars["EnergyRes"] * sqrt(track_energy);
    double delta_energy = ECAL_energy + HCAL_energy - track_energy;


    if(delta_energy > settings.map_floatPars["VirtualMergeSigma"]*sigmaE) continue;

    // Virtual HCAL clusters in the neutral PFO
    std::vector<const Cyber::Calo3DCluster*> all_neutral_HCAL_clus; all_neutral_HCAL_clus.clear();
    for(int ip=0; ip<m_neutralPFOs.size(); ip++){
      std::vector<const Cyber::Calo3DCluster*> tmp_HCAL_clus = m_neutralPFOs[ip]->getHCALClusters();
      if(tmp_HCAL_clus.size()!=1) continue;
      if(tmp_HCAL_clus[0]->getType()!=-1) continue;
      all_neutral_HCAL_clus.push_back(tmp_HCAL_clus[0]);
    }
//cout<<"    Virtual HCAL cluster size "<<all_neutral_HCAL_clus.size()<<", print them: "<<endl;
//for(int aa=0; aa<all_neutral_HCAL_clus.size(); aa++){
//  printf("      #%d: pos (%.3f, %.3f, %.3f), En %.3f, Nhit %d, type %d \n", aa, all_neutral_HCAL_clus[aa]->getHitCenter().x(), all_neutral_HCAL_clus[aa]->getHitCenter().y(), all_neutral_HCAL_clus[aa]->getHitCenter().z(), all_neutral_HCAL_clus[aa]->getHitsE()*settings.map_floatPars["HCALCalib"], all_neutral_HCAL_clus[aa]->getCaloHits().size(), all_neutral_HCAL_clus[aa]->getType() );
//}

    TVector3 trackclus_pos(0, 0, 0);
    if(m_chargedPFOs[ic]->getTracks()[0]->getTrackStates("Hcal").size()>0){
      std::vector<TrackState> m_extTrkStats = m_chargedPFOs[ic]->getTracks()[0]->getTrackStates("Hcal");

      int min_hit_index = -1;
      double min_distance = 999999;
      for(int i=0; i<m_extTrkStats.size(); i++){
        double hit_distance = m_extTrkStats[i].referencePoint.Perp();
        if(hit_distance<min_distance){
          min_distance = hit_distance;
          min_hit_index = i;
        }
      }
      
      if(min_hit_index>=0) trackclus_pos = m_extTrkStats[min_hit_index].referencePoint;
    }
    else if(m_chargedPFOs[ic]->getECALClusters().size()>0){
      trackclus_pos = m_chargedPFOs[ic]->getECALClusters()[0]->getShowerCenter();
    }
    else{
      for(int jc=0; jc<m_chargedPFOs[ic]->getHCALClusters().size(); jc++)
        trackclus_pos += m_chargedPFOs[ic]->getHCALClusters()[jc]->getHitCenter() * settings.map_floatPars["HCALChargedCalib"]*m_chargedPFOs[ic]->getHCALClusters()[jc]->getHitsE();

      trackclus_pos = trackclus_pos*(1./HCAL_energy);
    }
// printf("    charged cluster trk pos: [%.3f, %.3f, %.3f] \n", trackclus_pos.x(), trackclus_pos.y(), trackclus_pos.z());

    int loop_count = 0;
    std::vector<const Cyber::Calo3DCluster*> skip_clus;
    while(delta_energy < sigmaE*settings.map_floatPars["VirtualMergeSigma"]-(1e-6)){
      loop_count++;
      if(loop_count>all_neutral_HCAL_clus.size()+10){
        break;
      }
      double min_angle = 999.0;
      int clus_index = -1;

      for(int in=0; in<all_neutral_HCAL_clus.size(); in++){
        if(find(skip_clus.begin(), skip_clus.end(), all_neutral_HCAL_clus[in]) != skip_clus.end()) continue;

        TVector3 neutral_clus_pos = all_neutral_HCAL_clus[in]->getHitCenter();
        double pfo_angle = trackclus_pos.Angle(neutral_clus_pos);
        if (pfo_angle<settings.map_floatPars["MinAngleForVirMerge"] && pfo_angle<min_angle){
          min_angle=pfo_angle;
          clus_index = in;
        }
      }
// cout<<"    In Loop "<<loop_count<<": current deltaE = "<<delta_energy<<", closest virtual cluster index "<<clus_index<<endl;

      if(clus_index<0) break;  // No neutral Hcal cluster to be merged

      double tmp_delta_E = delta_energy + settings.map_floatPars["HCALNeutralCalib"]*all_neutral_HCAL_clus[clus_index]->getHitsE();
// cout<<"    If include this cluster: new deltaE "<<tmp_delta_E<<", merge = "<<(tmp_delta_E > sigmaE * settings.map_floatPars["VirtualMergeSigma"])<<endl;

      if(tmp_delta_E > sigmaE * settings.map_floatPars["VirtualMergeSigma"]){
        double absorbed_energy = sigmaE*settings.map_floatPars["VirtualMergeSigma"] - delta_energy;
        delta_energy = delta_energy + absorbed_energy;
// cout<<"      Absorberd energy: "<<absorbed_energy<<", current delta energy "<<delta_energy<<endl;

        //Create a new virtual neutral cluster with energy = absorbed_energy. 

        std::shared_ptr<Cyber::CaloHit> m_hit = all_neutral_HCAL_clus[clus_index]->getCaloHits()[0]->Clone();
        m_hit->setEnergy(absorbed_energy/settings.map_floatPars["HCALNeutralCalib"]);

        std::shared_ptr<Cyber::Calo3DCluster> m_clus = std::make_shared<Cyber::Calo3DCluster>();
        m_clus->addHit(m_hit.get());
        m_clus->setType(-1);

        m_bkCol.map_CaloHit["bkHit"].push_back( m_hit );
        m_bkCol.map_CaloCluster["bk3DCluster"].push_back(m_clus);

        m_chargedPFOs[ic]->addHCALCluster( m_clus.get() );


        //Re-set neutral virtual cluster energy
        bool is_found = false;
        for(int ip=0; ip<m_neutralPFOs.size(); ip++){
          auto tmp_HCAL_clus = m_neutralPFOs[ip]->getHCALClusters();

          if(tmp_HCAL_clus.size()!=1) continue;
          if(tmp_HCAL_clus[0]->getType()!=-1) continue;

          if(tmp_HCAL_clus[0]==all_neutral_HCAL_clus[clus_index]){
            std::shared_ptr<Cyber::CaloHit> m_newhit = all_neutral_HCAL_clus[clus_index]->getCaloHits()[0]->Clone();
            m_newhit->setEnergy(all_neutral_HCAL_clus[clus_index]->getHitsE() - absorbed_energy/settings.map_floatPars["HCALNeutralCalib"] );
            std::shared_ptr<Cyber::Calo3DCluster> m_newclus = std::make_shared<Cyber::Calo3DCluster>();
            m_newclus->addHit(m_newhit.get());
            m_newclus->setType(-1);            

            m_bkCol.map_CaloHit["bkHit"].push_back( m_newhit );
            m_bkCol.map_CaloCluster["bk3DCluster"].push_back(m_newclus);

            std::vector<const Calo3DCluster*> tmp_clusters; tmp_clusters.clear();
            tmp_clusters.push_back(m_newclus.get());
            m_neutralPFOs[ip]->setHCALCluster( tmp_clusters );

            is_found = true;
            break;
          }

        }
        if(!is_found){
          cout << "Error! Can not find the HCAL cluster in neutral PFO to delete part of the energy" << endl;
        }
        //Reset the energy
        all_neutral_HCAL_clus.erase(all_neutral_HCAL_clus.begin()+clus_index);
        all_neutral_HCAL_clus.push_back(m_clus.get());
        skip_clus.push_back(m_clus.get());
      }
      else{

        if (TMath::Abs(tmp_delta_E) >= TMath::Abs(delta_energy)){
          skip_clus.push_back(all_neutral_HCAL_clus[clus_index]);
          continue;  // No need to merge this HCAL cluster
        }
   
        // Update delta_energy
        delta_energy = tmp_delta_E;
        // Add this HCAL cluster to charged PFO
        m_chargedPFOs[ic]->addHCALCluster(all_neutral_HCAL_clus[clus_index]);
   
        // Remove this HCAL cluster from neutral PFO
        bool is_found = false;
        for(int in=0; in<m_neutralPFOs.size(); in++){
          std::vector<const Cyber::Calo3DCluster*> neutral_cluster = m_neutralPFOs[in]->getHCALClusters();
          int tmp_index=-1;
          for(int ii=0; ii<neutral_cluster.size(); ii++){
            if (all_neutral_HCAL_clus[clus_index]==neutral_cluster[ii]){
              tmp_index = ii;
              break;
            }
          }
          if (tmp_index==-1) continue;
// cout<<"  Remove a neutral cluster: En "<<settings.map_floatPars["HCALCalib"]*neutral_cluster[tmp_index]->getHitsE()<<endl;

          neutral_cluster.erase(neutral_cluster.begin()+tmp_index);
          m_neutralPFOs[in]->setHCALCluster(neutral_cluster);
   
          if(m_neutralPFOs[in]->getTracks().size() + m_neutralPFOs[in]->getECALClusters().size() + m_neutralPFOs[in]->getHCALClusters().size()==0){
// cout<<"  Remove a neutral PFO: ECAL En "<<m_neutralPFOs[in]->getECALClusterEnergy()<<", HCAL En "<<m_neutralPFOs[in]->getHCALClusterEnergy()<<endl;
            auto iter = find(p_PFObjects->begin(), p_PFObjects->end(), m_neutralPFOs[in]);
            if(iter==p_PFObjects->end()){
              std::cout<<"ERROR: can not find this neutral PFO in p_PFObjects. "<<std::endl;
            }
            else{
              m_neutralPFOs.erase(m_neutralPFOs.begin()+in);
              p_PFObjects->erase(iter);
            }
          }
          is_found = true;
          break;
        }
        if(!is_found){
          cout << "Error! Can not find the HCAL cluster in neutral PFO" << endl;
        }
   
        // Remove this HCAL cluster from all_neutral_HCAL_clus
        all_neutral_HCAL_clus.erase(all_neutral_HCAL_clus.begin()+clus_index);

      }

    }
  }

  return StatusCode::SUCCESS;
};


StatusCode PFOReclusteringAlg::ReCluster_SplitFromChg( std::vector< std::shared_ptr<Cyber::PFObject> >& m_chargedPFOs,
                                                       std::vector< std::shared_ptr<Cyber::PFObject> >& m_neutralPFOs ){

  for(int ipfo=0; ipfo<m_chargedPFOs.size(); ipfo++){
    if(m_chargedPFOs[ipfo]->getECALClusters().size()==0 && m_chargedPFOs[ipfo]->getHCALClusters().size()==0) continue;

    double track_energy = m_chargedPFOs[ipfo]->getTrackMomentum();
    double ECAL_energy = settings.map_floatPars["ECALChargedEMCalib"]*m_chargedPFOs[ipfo]->getECALClusterEnergy();
    double HCAL_energy = settings.map_floatPars["HCALChargedCalib"]*m_chargedPFOs[ipfo]->getHCALClusterEnergy();

//cout<<"  Charged PFO #"<<ipfo<<": track P "<<track_energy<<", cluster E "<<ECAL_energy+HCAL_energy<<", PID "<<m_chargedPFOs[ipfo]->getPID()<<endl;

    //FIXME: If PFO is identified as electron (E/p > 0.8): use ECALChargedEMCalib(1.0). Else: use ECALChargedCalib(1.26). 
    if( track_energy < 15. || ECAL_energy/track_energy < settings.map_floatPars["ElecEoverP"] ) 
    //if( (track_energy<15. && abs(m_chargedPFOs[ipfo]->getPID())!=11) || 
    //    (track_energy>15. && ECAL_energy/track_energy < settings.map_floatPars["ElecEoverP"]) )
      ECAL_energy = ECAL_energy*settings.map_floatPars["ECALChargedCalib"]/settings.map_floatPars["ECALChargedEMCalib"]; 

    if(track_energy<0 || ECAL_energy<0 || HCAL_energy<0){
      std::cout<<"ERROR: Charged PFO info break. Ptrk "<<track_energy<<", E_ecal "<<ECAL_energy<<", E_hcal "<<HCAL_energy<<endl;
      continue;
    }

    double delta_energy = ECAL_energy + HCAL_energy - track_energy;
    double sigmaE = settings.map_floatPars["EnergyRes"] * sqrt(track_energy);
    
//cout<<"  ReCluster_MergeToChg: In ChPFO #"<<ipfo<<": Ptrk = "<<track_energy<<", Eecal = "<<ECAL_energy<<", Ehcal = "<<HCAL_energy<<", deltaE = "<<delta_energy<<", sigmaE = "<<sigmaE<<endl;

    if(delta_energy <= settings.map_floatPars["SplitSigma"]*sigmaE) continue;
//cout<<"    Do charged PFO splitting. "<<endl;

    //Create a new hit and cluster
    TVector3 tmp_pos(0,0,0);
    for(int ic=0; ic<m_chargedPFOs[ipfo]->getECALClusters().size(); ic++)
      tmp_pos += settings.map_floatPars["ECALChargedCalib"] * m_chargedPFOs[ipfo]->getECALClusters()[ic]->getLongiE() * m_chargedPFOs[ipfo]->getECALClusters()[ic]->getShowerCenter();
    for(int ic=0; ic<m_chargedPFOs[ipfo]->getHCALClusters().size(); ic++)
      tmp_pos += settings.map_floatPars["HCALChargedCalib"] * m_chargedPFOs[ipfo]->getHCALClusters()[ic]->getHitsE() * m_chargedPFOs[ipfo]->getHCALClusters()[ic]->getHitCenter();
    tmp_pos = tmp_pos*(1./(ECAL_energy+HCAL_energy));

    if(ECAL_energy>0 || HCAL_energy>0){
      //Create a new ECAL virtual PFO
      if(ECAL_energy>0){
        double tmp_virtualE = delta_energy * ECAL_energy / (ECAL_energy + HCAL_energy) / settings.map_floatPars["ECALChargedCalib"];
        // -- Create bar
        std::shared_ptr<Cyber::CaloUnit> m_bar_u = std::make_shared<Cyber::CaloUnit>();
        std::shared_ptr<Cyber::CaloUnit> m_bar_v = std::make_shared<Cyber::CaloUnit>();
        m_bar_u->setQ(tmp_virtualE/4., tmp_virtualE/4);
        m_bar_u->setPosition(tmp_pos);
        m_bar_v->setQ(tmp_virtualE/4., tmp_virtualE/4);
        m_bar_v->setPosition(tmp_pos);
   
        // -- Create 1D cluster U/V
        std::shared_ptr<Cyber::Calo1DCluster> m_1dclus_u = std::make_shared<Cyber::Calo1DCluster>();
        std::shared_ptr<Cyber::Calo1DCluster> m_1dclus_v = std::make_shared<Cyber::Calo1DCluster>();
        m_1dclus_u->addUnit(m_bar_u.get());
        m_1dclus_v->addUnit(m_bar_v.get());
   
        // -- Create 2D cluster 
        std::shared_ptr<Cyber::Calo2DCluster> m_2dclus = std::make_shared<Cyber::Calo2DCluster>();
        m_2dclus->addShowerU(m_1dclus_u.get());
        m_2dclus->addShowerV(m_1dclus_v.get());
        m_2dclus->addTowerID(0., 0., 0., 0.);
        m_2dclus->setPos(tmp_pos);
   
        // -- Create Half cluster U/V
        std::shared_ptr<Cyber::CaloHalfCluster> m_hfclus_u = std::make_shared<Cyber::CaloHalfCluster>();
        std::shared_ptr<Cyber::CaloHalfCluster> m_hfclus_v = std::make_shared<Cyber::CaloHalfCluster>();
        m_hfclus_u->addUnit(m_1dclus_u.get());
        m_hfclus_v->addUnit(m_1dclus_v.get());
   
        // -- Create ECAL 3D cluster
        std::shared_ptr<Cyber::Calo3DCluster> m_clus = std::make_shared<Cyber::Calo3DCluster>();
        m_clus->addHalfClusterU("LinkedLongiCluster", m_hfclus_u.get());
        m_clus->addHalfClusterV("LinkedLongiCluster", m_hfclus_v.get());
        m_clus->addUnit(m_2dclus.get());
        m_clus->setType(-1);
   
        std::shared_ptr<Cyber::PFObject> m_pfo = std::make_shared<Cyber::PFObject>();
        m_pfo->addECALCluster( m_clus.get() );
   
        m_neutralPFOs.push_back( m_pfo );
        p_PFObjects->push_back( m_pfo );
        
//cout<<"  In ChPFO #"<<ipfo<<": create a neutral ECAL PFO. En "<<m_pfo->getECALClusterEnergy()*settings.map_floatPars["ECALCalib"]<<endl;
   
        m_bkCol.map_BarCol["bkBar"].push_back( m_bar_u );
        m_bkCol.map_BarCol["bkBar"].push_back( m_bar_v );
        m_bkCol.map_1DCluster["bk1DCluster"].push_back( m_1dclus_u );
        m_bkCol.map_1DCluster["bk1DCluster"].push_back( m_1dclus_v );
        m_bkCol.map_HalfCluster["bkHalfCluster"].push_back( m_hfclus_u );
        m_bkCol.map_HalfCluster["bkHalfCluster"].push_back( m_hfclus_v );
        m_bkCol.map_2DCluster["bk2DCluster"].push_back( m_2dclus );
        m_bkCol.map_CaloCluster["bk3DCluster"].push_back(m_clus);
        m_bkCol.map_PFObjects["bkPFO"].push_back(m_pfo);
      }

      //Create a new virtual HCAL PFO
      if(HCAL_energy>0){
        double tmp_virtualE = delta_energy * HCAL_energy / (ECAL_energy + HCAL_energy) / settings.map_floatPars["HCALChargedCalib"];
        std::shared_ptr<Cyber::CaloHit> m_hit = std::make_shared<Cyber::CaloHit>();
        m_hit->setPosition( tmp_pos );
        m_hit->setEnergy( tmp_virtualE );

        std::shared_ptr<Cyber::Calo3DCluster> m_clus_hcal = std::make_shared<Cyber::Calo3DCluster>();
        m_clus_hcal->addHit(m_hit.get());
        m_clus_hcal->setType(-1);

        std::shared_ptr<Cyber::PFObject> m_pfo_hcal = std::make_shared<Cyber::PFObject>();
        m_pfo_hcal->addHCALCluster( m_clus_hcal.get() );

        m_neutralPFOs.push_back( m_pfo_hcal );
        p_PFObjects->push_back( m_pfo_hcal );

//cout<<"  In ChPFO #"<<ipfo<<": create a neutral HCAL PFO. En "<<m_pfo_hcal->getHCALClusterEnergy()*settings.map_floatPars["HCALCalib"]<<endl;

        m_bkCol.map_CaloHit["bkHit"].push_back( m_hit );
        m_bkCol.map_CaloCluster["bk3DCluster"].push_back(m_clus_hcal);
        m_bkCol.map_PFObjects["bkPFO"].push_back(m_pfo_hcal);
      }

      //Reset origin cluster energy
      double m_EnScale = (ECAL_energy + HCAL_energy - delta_energy)/(ECAL_energy + HCAL_energy);
//cout<<"[FY debug] In PFO "<<ipfo<<": track P "<<m_chargedPFOs[ipfo]->getTrackMomentum()<<", HCAL cluster size "<<m_chargedPFOs[ipfo]->getHCALClusters().size()<<endl;
//", first HCAL has Nhit "<<m_chargedPFOs[ipfo]->getHCALClusters()[0]->getCaloHits().size()<<endl;
//cout<<"[FY debug] Hcal energy scale = "<<m_EnScale<<endl;

      //Create a new PFO
      std::shared_ptr<Cyber::PFObject> m_newpfo = m_chargedPFOs[ipfo]->Clone();
      // -- Reset ECAL cluster
      std::vector<const Cyber::Calo3DCluster*> tmp_clusvec; tmp_clusvec.clear();
      for(int ic=0; ic<m_chargedPFOs[ipfo]->getECALClusters().size(); ic++){
        std::shared_ptr<Cyber::Calo3DCluster> m_newclus = m_chargedPFOs[ipfo]->getECALClusters()[ic]->Clone();
        m_newclus->setEnergyScale( m_EnScale*m_newclus->getEnergyScale() );
        tmp_clusvec.push_back(m_newclus.get()); 
        m_bkCol.map_CaloCluster["bk3DCluster"].push_back(m_newclus);
      }
      m_newpfo->setECALCluster(tmp_clusvec);
    
      // -- Reset HCAL cluster
      std::shared_ptr<Cyber::Calo3DCluster> m_newclus = std::make_shared<Cyber::Calo3DCluster>();
      for(int ic=0; ic<m_chargedPFOs[ipfo]->getHCALClusters().size(); ic++){
        std::vector<const Cyber::CaloHit*> tmp_hits = m_chargedPFOs[ipfo]->getHCALClusters()[ic]->getCaloHits();
        for(int ih=0; ih<tmp_hits.size(); ih++){
          std::shared_ptr<CaloHit> tmp_newhit = tmp_hits[ih]->Clone();
          tmp_newhit->setEnergy( tmp_newhit->getEnergy()*m_EnScale );
          m_newclus->addHit(tmp_newhit.get());
          m_bkCol.map_CaloHit["bkHit"].push_back( tmp_newhit );
        }
      }
      m_bkCol.map_CaloCluster["bk3DCluster"].push_back(m_newclus);
      std::vector<const Cyber::Calo3DCluster*> tmp_hcalclus; tmp_hcalclus.clear();
      tmp_hcalclus.push_back(m_newclus.get());
      m_newpfo->setHCALCluster(tmp_hcalclus);

      m_bkCol.map_PFObjects["bkPFO"].push_back(m_newpfo);


      auto iter = find(p_PFObjects->begin(), p_PFObjects->end(), m_chargedPFOs[ipfo]);
      if(iter!=p_PFObjects->end())
        *iter = m_newpfo;

      m_chargedPFOs[ipfo] = m_newpfo;
//cout<<"[FY debug] The splitted new charged PFO: track size "<<m_newpfo->getTracks().size()<<", leading P "<<m_newpfo->getTrackMomentum();
//cout<<", ECAL cluster size "<<m_newpfo->getECALClusters().size()<<", totE "<<m_newpfo->getECALClusterEnergy();
//cout<<", HCAL cluster size "<<m_newpfo->getHCALClusters().size()<<", totE "<<m_newpfo->getHCALClusterEnergy()<<endl;
//cout<<endl;
    }
  }

  return StatusCode::SUCCESS;
};

#endif
