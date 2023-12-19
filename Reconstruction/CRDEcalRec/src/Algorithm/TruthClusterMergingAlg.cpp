#ifndef _TRUTHCLUSTERMERGING_ALG_C
#define _TRUTHCLUSTERMERGING_ALG_C

#include "Algorithm/TruthClusterMergingAlg.h"

StatusCode TruthClusterMergingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_boolPars.find("DoECALMerge")==settings.map_boolPars.end()) settings.map_boolPars["DoECALMerge"] = 1;
  if(settings.map_boolPars.find("DoHCALMerge")==settings.map_boolPars.end()) settings.map_boolPars["DoHCALMerge"] = 1;
  if(settings.map_boolPars.find("DoECALHCALConnection")==settings.map_boolPars.end()) settings.map_boolPars["DoECALHCALConnection"] = 1;

  if(settings.map_stringPars.find("ReadinECALClusters")==settings.map_stringPars.end()) settings.map_stringPars["ReadinECALClusters"] = "TruthEcalCluster";
  if(settings.map_stringPars.find("ReadinHCALClusters")==settings.map_stringPars.end()) settings.map_stringPars["ReadinHCALClusters"] = "TruthHcalCluster";
  if(settings.map_stringPars.find("OutputECALCluster")==settings.map_stringPars.end()) settings.map_stringPars["OutputECALCluster"] = "TruthMergedEcalCluster";
  if(settings.map_stringPars.find("OutputHCALCluster")==settings.map_stringPars.end()) settings.map_stringPars["OutputHCALCluster"] = "TruthMergedHcalCluster";
  if(settings.map_stringPars.find("OutputCombPFO")==settings.map_stringPars.end()) settings.map_stringPars["OutputCombPFO"] = "TruthCombPFO";

  return StatusCode::SUCCESS;
};

StatusCode TruthClusterMergingAlg::Initialize( PandoraPlusDataCol& m_datacol ){

  for(int icl=0; icl<m_datacol.map_CaloCluster[settings.map_stringPars["ReadinECALClusters"]].size(); icl++)
    m_EcalClusterCol.push_back(m_datacol.map_CaloCluster[settings.map_stringPars["ReadinECALClusters"]][icl].get()); 
  for(int icl=0; icl<m_datacol.map_CaloCluster[settings.map_stringPars["ReadinHCALClusters"]].size(); icl++)
    m_HcalClusterCol.push_back(m_datacol.map_CaloCluster[settings.map_stringPars["ReadinHCALClusters"]][icl].get()); 


cout<<"Print input ECAL cluster"<<endl;
for(int ic=0; ic<m_EcalClusterCol.size(); ic++){
  printf("  Pos+E (%.3f, %.3f, %.3f, %.3f), linked leading MCP [%d, %.3f], linked track size %d \n", m_EcalClusterCol[ic]->getShowerCenter().x(), m_EcalClusterCol[ic]->getShowerCenter().y(), m_EcalClusterCol[ic]->getShowerCenter().z(), m_EcalClusterCol[ic]->getLongiE(), m_EcalClusterCol[ic]->getLeadingMCP().getPDG(), m_EcalClusterCol[ic]->getLeadingMCPweight(), m_EcalClusterCol[ic]->getAssociatedTracks().size() );
  if(m_EcalClusterCol[ic]->getAssociatedTracks().size()>0) printf("    Track MCP: %d \n", m_EcalClusterCol[ic]->getAssociatedTracks()[0]->getLeadingMCP().getPDG());  
}
cout<<"Print input HCAL cluster"<<endl;
for(int ic=0; ic<m_HcalClusterCol.size(); ic++){
  printf("  Pos+E (%.3f, %.3f, %.3f, %.3f), linked leading MCP [%d, %.3f], linked track size %d \n", m_HcalClusterCol[ic]->getHitCenter().x(), m_HcalClusterCol[ic]->getHitCenter().y(), m_HcalClusterCol[ic]->getHitCenter().z(), m_HcalClusterCol[ic]->getHitsE(), m_HcalClusterCol[ic]->getLeadingMCP().getPDG(), m_HcalClusterCol[ic]->getLeadingMCPweight(), m_HcalClusterCol[ic]->getAssociatedTracks().size() );
  if(m_HcalClusterCol[ic]->getAssociatedTracks().size()>0) printf("    Track MCP: %d \n", m_HcalClusterCol[ic]->getAssociatedTracks()[0]->getLeadingMCP().getPDG());
}


  return StatusCode::SUCCESS;
};

StatusCode TruthClusterMergingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

//cout<<"  TruthClusterMergingAlg: input ECAL cluster size "<<m_EcalClusterCol.size()<<", HCAL cluster size "<<m_HcalClusterCol.size()<<endl;
//cout<<"  Print ECAL cluster energy: "<<endl;
//for(int ic=0; ic<m_EcalClusterCol.size(); ic++)
//  cout<<"    #"<<ic<<" En = "<<m_EcalClusterCol[ic]->getLongiE()<<endl;
//cout<<"  Print HCAL cluster energy: "<<endl;
//for(int ic=0; ic<m_HcalClusterCol.size(); ic++)
//  cout<<"    #"<<ic<<" En = "<<m_HcalClusterCol[ic]->getHitsE()<<endl;

  std::map<edm4hep::MCParticle, std::vector<const PandoraPlus::Calo3DCluster*>> map_ClusterCol_Ecal; 
  std::map<edm4hep::MCParticle, std::vector<const PandoraPlus::Calo3DCluster*>> map_ClusterCol_Hcal; 
  for(int ic=0; ic<m_EcalClusterCol.size(); ic++){
    auto mcp = m_EcalClusterCol[ic]->getLeadingMCP();
    map_ClusterCol_Ecal[mcp].push_back( m_EcalClusterCol[ic] );
  }
  for(int ic=0; ic<m_HcalClusterCol.size(); ic++){
    auto mcp = m_HcalClusterCol[ic]->getLeadingMCP();
    map_ClusterCol_Hcal[mcp].push_back( m_HcalClusterCol[ic] );
  }

  if(settings.map_boolPars["DoECALMerge"]){
    for(auto& iter: map_ClusterCol_Ecal){
      if(iter.second.size()==0) continue;
      auto tmp_newcluster = iter.second[0]->Clone();
      for(int icl=1; icl<iter.second.size(); icl++){
        tmp_newcluster->mergeCluster( iter.second[icl] );
      }
      //tmp_newcluster->getLinkedMCPfromUnit();
      merged_EcalClusterCol.push_back(tmp_newcluster);
      m_datacol.map_CaloCluster["bk3DCluster"].push_back(tmp_newcluster);
    }
  }

  if(settings.map_boolPars["DoHCALMerge"]){
    for(auto& iter: map_ClusterCol_Hcal){
      if(iter.second.size()==0) continue;
      auto tmp_newcluster = iter.second[0]->Clone();
      for(int icl=1; icl<iter.second.size(); icl++){
        tmp_newcluster->mergeCluster( iter.second[icl] );
      }
      //tmp_newcluster->getLinkedMCPfromUnit();
      merged_HcalClusterCol.push_back(tmp_newcluster);
      m_datacol.map_CaloCluster["bk3DCluster"].push_back(tmp_newcluster);
    }
  }


  if(settings.map_boolPars["DoECALHCALConnection"]){
    map_ClusterCol_Ecal.clear();
    map_ClusterCol_Hcal.clear();

    if(settings.map_boolPars["DoECALMerge"]){
      for(int ic=0; ic<merged_EcalClusterCol.size(); ic++){
        auto mcp = merged_EcalClusterCol[ic]->getLeadingMCP();
        map_ClusterCol_Ecal[mcp].push_back( merged_EcalClusterCol[ic].get() );
      }
    }
    else{
      for(int ic=0; ic<m_EcalClusterCol.size(); ic++){
        auto mcp = m_EcalClusterCol[ic]->getLeadingMCP();
        map_ClusterCol_Ecal[mcp].push_back( m_EcalClusterCol[ic] );
      }
    }

    if(settings.map_boolPars["DoHCALMerge"]){
      for(int ic=0; ic<merged_HcalClusterCol.size(); ic++){
        auto mcp = merged_HcalClusterCol[ic]->getLeadingMCP();
        map_ClusterCol_Hcal[mcp].push_back( merged_HcalClusterCol[ic].get() );
      }
    }
    else{
      for(int ic=0; ic<m_HcalClusterCol.size(); ic++){
        auto mcp = m_HcalClusterCol[ic]->getLeadingMCP();
        map_ClusterCol_Hcal[mcp].push_back( m_HcalClusterCol[ic] );
      }
    }


    for(auto& iter: map_ClusterCol_Ecal){
      auto mcp = iter.first;
      std::vector<const PandoraPlus::Calo3DCluster*> m_linkedEcalClus = iter.second;
      std::vector<const PandoraPlus::Calo3DCluster*> m_linkedHcalClus = map_ClusterCol_Hcal[mcp];   

      if(m_linkedEcalClus.size()==0 && m_linkedHcalClus.size()==0) continue;

      std::shared_ptr<PandoraPlus::PFObject> tmp_newpfo = std::make_shared<PandoraPlus::PFObject>();
      double Emax = -99;
      int index = -1;
      for(int ic=0; ic<m_linkedEcalClus.size(); ic++){ 
        tmp_newpfo->addECALCluster(m_linkedEcalClus[ic]);
        if(m_linkedEcalClus[ic]->getLongiE()>Emax){
          Emax = m_linkedEcalClus[ic]->getLongiE();
          index = ic;
        }
      }
      const PandoraPlus::Track* p_EcalLeadingTrk = nullptr;
      if(index>=0 && m_linkedEcalClus[index]->getAssociatedTracks().size()>0) p_EcalLeadingTrk = m_linkedEcalClus[index]->getAssociatedTracks()[0];

      Emax = -99; index = -1;
      for(int ic=0; ic<m_linkedHcalClus.size(); ic++){
        tmp_newpfo->addHCALCluster(m_linkedHcalClus[ic]);
        if(m_linkedHcalClus[ic]->getHitsE()>Emax){
          Emax = m_linkedHcalClus[ic]->getHitsE();
          index = ic;
        }
      }
      const PandoraPlus::Track* p_HcalLeadingTrk = nullptr;
      if(index>=0 && m_linkedHcalClus[index]->getAssociatedTracks().size()>0) p_HcalLeadingTrk = m_linkedHcalClus[index]->getAssociatedTracks()[0];

      if(!p_EcalLeadingTrk && p_HcalLeadingTrk) tmp_newpfo->addTrack(p_HcalLeadingTrk);
      else if(p_EcalLeadingTrk) tmp_newpfo->addTrack(p_EcalLeadingTrk);

      merged_CombClusterCol.push_back(tmp_newpfo);
      m_datacol.map_PFObjects["bkPFO"].push_back(tmp_newpfo);
    }

  }
 /* 
cout<<"Print PFO"<<endl;
for(int ip=0; ip<merged_CombClusterCol.size(); ip++){
  printf("  PFO #%d: track size %d, ECAL cluster size %d, HCAL cluster size %d \n", ip, merged_CombClusterCol[ip]->getTracks().size(), merged_CombClusterCol[ip]->getECALClusters().size(), merged_CombClusterCol[ip]->getHCALClusters().size() );
  cout<<"    Track: "<<endl;
  for(int itrk=0; itrk<merged_CombClusterCol[ip]->getTracks().size(); itrk++) 
    printf("      Trk #%d: momentum %.3f, linked MCP %d \n", itrk, merged_CombClusterCol[ip]->getTracks()[itrk]->getMomentum(), merged_CombClusterCol[ip]->getTracks()[itrk]->getLeadingMCP().getPDG() );
  cout<<"    ECAL cluster: "<<endl;
  for(int icl=0; icl<merged_CombClusterCol[ip]->getECALClusters().size(); icl++)
    printf("      ECAL cluster #%d: energy %.3f, linked MCP %d \n", icl, merged_CombClusterCol[ip]->getECALClusters()[icl]->getLongiE(), merged_CombClusterCol[ip]->getECALClusters()[icl]->getLeadingMCP().getPDG() );
  cout<<"    HCAL cluster: "<<endl;
  for(int icl=0; icl<merged_CombClusterCol[ip]->getHCALClusters().size(); icl++)
    printf("      HCAL cluster #%d: energy %.3f, linked MCP %d \n", icl, merged_CombClusterCol[ip]->getHCALClusters()[icl]->getHitsE(), merged_CombClusterCol[ip]->getHCALClusters()[icl]->getLeadingMCP().getPDG() );
  cout<<endl;
}
*/

  m_datacol.map_CaloCluster[settings.map_stringPars["OutputECALCluster"]] = merged_EcalClusterCol;
  m_datacol.map_CaloCluster[settings.map_stringPars["OutputHCALCluster"]] = merged_HcalClusterCol;
  m_datacol.map_PFObjects[settings.map_stringPars["OutputCombPFO"]] = merged_CombClusterCol;

  return StatusCode::SUCCESS;
};

StatusCode TruthClusterMergingAlg::ClearAlgorithm(){

  m_EcalClusterCol.clear();
  m_HcalClusterCol.clear();
  merged_EcalClusterCol.clear();
  merged_HcalClusterCol.clear();
  merged_CombClusterCol.clear();

  return StatusCode::SUCCESS;
};


#endif
