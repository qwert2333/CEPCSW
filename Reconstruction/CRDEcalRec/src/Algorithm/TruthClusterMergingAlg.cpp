#ifndef _TRUTHCLUSTERMERGING_ALG_C
#define _TRUTHCLUSTERMERGING_ALG_C

#include "Algorithm/TruthClusterMergingAlg.h"

StatusCode TruthClusterMergingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_stringPars.find("ReadinECALClusters")==settings.map_stringPars.end()) settings.map_stringPars["ReadinECALClusters"] = "TruthEcalCluster";
  if(settings.map_stringPars.find("ReadinHCALClusters")==settings.map_stringPars.end()) settings.map_stringPars["ReadinHCALClusters"] = "TruthHcalCluster";
  if(settings.map_stringPars.find("OutputECALCluster")==settings.map_stringPars.end()) settings.map_stringPars["OutputECALCluster"] = "TruthMergedEcalCluster";
  if(settings.map_stringPars.find("OutputHCALCluster")==settings.map_stringPars.end()) settings.map_stringPars["OutputHCALCluster"] = "TruthMergedHcalCluster";
  if(settings.map_stringPars.find("OutputCombCluster")==settings.map_stringPars.end()) settings.map_stringPars["OutputCombCluster"] = "TruthCombCluster";

  return StatusCode::SUCCESS;
};

StatusCode TruthClusterMergingAlg::Initialize( PandoraPlusDataCol& m_datacol ){

  for(int icl=0; icl<m_datacol.map_CaloCluster[settings.map_stringPars["ReadinECALClusters"]].size(); icl++)
    m_EcalClusterCol.push_back(m_datacol.map_CaloCluster[settings.map_stringPars["ReadinECALClusters"]][icl].get()); 
  for(int icl=0; icl<m_datacol.map_CaloCluster[settings.map_stringPars["ReadinHCALClusters"]].size(); icl++)
    m_HcalClusterCol.push_back(m_datacol.map_CaloCluster[settings.map_stringPars["ReadinHCALClusters"]][icl].get()); 

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
  std::map<edm4hep::MCParticle, std::vector<const PandoraPlus::Calo3DCluster*>> map_ClusterCol_Comb; 
  for(int ic=0; ic<m_EcalClusterCol.size(); ic++){
    auto mcp = m_EcalClusterCol[ic]->getLeadingMCP();
    map_ClusterCol_Ecal[mcp].push_back( m_EcalClusterCol[ic] );
    map_ClusterCol_Comb[mcp].push_back( m_EcalClusterCol[ic] );
  }
  for(int ic=0; ic<m_HcalClusterCol.size(); ic++){
    auto mcp = m_HcalClusterCol[ic]->getLeadingMCP();
    map_ClusterCol_Hcal[mcp].push_back( m_HcalClusterCol[ic] );
    map_ClusterCol_Comb[mcp].push_back( m_HcalClusterCol[ic] );
  }


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

  for(auto& iter: map_ClusterCol_Comb){
    if(iter.second.size()==0) continue;
    auto tmp_newcluster = iter.second[0]->Clone();
    for(int icl=1; icl<iter.second.size(); icl++){
      tmp_newcluster->mergeCluster( iter.second[icl] );    
    }
    //tmp_newcluster->getLinkedMCPfromUnit();
    merged_CombClusterCol.push_back(tmp_newcluster);
    m_datacol.map_CaloCluster["bk3DCluster"].push_back(tmp_newcluster);
  }

//cout<<"  TruthClusterMergingAlg: after merge ECAL cluster size "<<merged_EcalClusterCol.size()<<", HCAL cluster size "<<merged_HcalClusterCol.size()<<endl;
//cout<<"  Print ECAL cluster energy: "<<endl;
//for(int ic=0; ic<merged_EcalClusterCol.size(); ic++)
//  cout<<"    #"<<ic<<" En = "<<merged_EcalClusterCol[ic]->getLongiE()<<endl;
//cout<<"  Print HCAL cluster energy: "<<endl;
//for(int ic=0; ic<merged_HcalClusterCol.size(); ic++)
//  cout<<"    #"<<ic<<" En = "<<merged_HcalClusterCol[ic]->getHitsE()<<endl;



  m_datacol.map_CaloCluster[settings.map_stringPars["OutputECALCluster"]] = merged_EcalClusterCol;
  m_datacol.map_CaloCluster[settings.map_stringPars["OutputHCALCluster"]] = merged_HcalClusterCol;
  m_datacol.map_CaloCluster[settings.map_stringPars["OutputCombCluster"]] = merged_CombClusterCol;

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
