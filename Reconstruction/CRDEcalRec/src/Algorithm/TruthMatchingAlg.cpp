#ifndef _TRUTHMATCHING_ALG_C
#define _TRUTHMATCHING_ALG_C

#include "Algorithm/TruthMatchingAlg.h"

StatusCode TruthMatchingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_stringPars.find("ReadinHFClusterName")==settings.map_stringPars.end()) settings.map_stringPars["ReadinHFClusterName"] = "ESHalfCluster";

  return StatusCode::SUCCESS;
};

StatusCode TruthMatchingAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  std::cout<<"Initialize TruthMatchingAlg"<<std::endl;

  return StatusCode::SUCCESS;
};

StatusCode TruthMatchingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  for(int it=0; it<m_towerCol.size(); it++){

    std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> tmp_clusters; tmp_clusters.clear();
    m_HFClusUCol = m_towerCol.at(it)->getHalfClusterUCol(settings.map_stringPars["ReadinHFClusterName"]+"U");
    m_HFClusVCol = m_towerCol.at(it)->getHalfClusterVCol(settings.map_stringPars["ReadinHFClusterName"]+"V");

prtf("  In Tower [%d, %d, %d], ", m_towerCol[it]->getTowerID()[0][0], m_towerCol[it]->getTowerID()[0][1], m_towerCol[it]->getTowerID()[0][2] );
prtf(" HalfCluster size (%d, %d), input total energy %.3f \n", m_HFClusUCol.size(), m_HFClusVCol.size(), m_towerCol[it]->getEnergy());

    TruthMatching(m_HFClusUCol, m_HFClusVCol, tmp_clusters);
    m_clusterCol.insert(m_clusterCol.end(), tmp_clusters.begin(), tmp_clusters.end());

co<<"  Reconstructed cluster energy: "<<endl;
foint i=0; i<tmp_clusters.size(); i++) printf("  Pos+E: (%.3f, %.3f, %.3f). En %.3f \n", tmp_clusters[i]->getShowerCenter().x(), tmp_clusters[i]->getShowerCenter().y(), tmp_clusters[i]->getShowerCenter().z(), tmp_clusters[i]->getEnergy() );
  }

  //Merge clusters linked to the same MCP.
  for(int ic=0; ic<m_clusterCol.size() && m_clusterCol.size()>1; ic++){
    edm4hep::MCParticle mcp1 = m_clusterCol[ic].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getLinkedMCP()[0].first;
    for(int jc=ic+1; jc<m_clusterCol.size(); jc++){
      if(ic>m_clusterCol.size()) ic--;
      edm4hep::MCParticle mcp2 = m_clusterCol[jc].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getLinkedMCP()[0].first;
      if(mcp1==mcp2){
        m_clusterCol[ic].get()->mergeCluster( m_clusterCol[jc].get() );
        m_clusterCol.erase(m_clusterCol.begin()+jc);
        jc--;
        if(jc<ic) jc=ic;
      }
    }
  }


cout<<"    After merge: cluster size "<<m_clusterCol.size()<<endl;
   for(int ic=0; ic<m_clusterCol.size(); ic++) m_clusterCol[ic].get()->getLinkedMCPfromHFCluster("LinkedLongiCluster");

cout<<"    After merge: cluster size "<<m_clusterCol.size()<<endl;
cout<<"    Print cluster MC truth info"<<endl;
for(int i=0; i<m_clusterCol.size(); i++){
  printf("  Cluster #%d: En = %.3f, Ntrk %d, Nmcp = %d \n ",i, m_clusterCol[i]->getLongiE(), m_clusterCol[i]->getAssociatedTracks().size(), m_clusterCol[i]->getLinkedMCP().size());
  for(auto imcp : m_clusterCol[i]->getLinkedMCP()){
    printf("    MCP mom (%.3f, %.3f, %.3f), PDG %d, weight %.3f \n", imcp.first.getMomentum().x, imcp.first.getMomentum().y, imcp.first.getMomentum().z, imcp.first.getPDG(), imcp.second);
  }

}

  m_datacol.map_CaloCluster["EcalCluster"] = m_clusterCol;

  m_datacol.map_BarCol["bkBar"].insert( m_datacol.map_BarCol["bkBar"].end(), m_bkCol.map_BarCol["bkBar"].begin(), m_bkCol.map_BarCol["bkBar"].end() );
  m_datacol.map_1DCluster["bk1DCluster"].insert( m_datacol.map_1DCluster["bk1DCluster"].end(), m_bkCol.map_1DCluster["bk1DCluster"].begin(), m_bkCol.map_1DCluster["bk1DCluster"].end() );
  m_datacol.map_2DCluster["bk2DCluster"].insert( m_datacol.map_2DCluster["bk2DCluster"].end(), m_bkCol.map_2DCluster["bk2DCluster"].begin(), m_bkCol.map_2DCluster["bk2DCluster"].end() );
  m_datacol.map_HalfCluster["bkHalfCluster"].insert( m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_bkCol.map_HalfCluster["bkHalfCluster"].begin(), m_bkCol.map_HalfCluster["bkHalfCluster"].end() );



  return StatusCode::SUCCESS;
};

StatusCode TruthMatchingAlg::ClearAlgorithm(){
  std::cout<<"End run TruthMatchingAlg. Clean it."<<std::endl;

  return StatusCode::SUCCESS;
};


#endif
