#ifndef _TRUTHTRKMATCHING_ALG_C
#define _TRUTHTRKMATCHING_ALG_C

#include "Algorithm/TruthTrackMatchingAlg.h"

StatusCode TruthTrackMatchingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_stringPars.find("MatchType")==settings.map_stringPars.end())
    settings.map_stringPars["MatchType"] = "LocalMax";  //LocalMax or Bar. TODO: implatement this option. 
  if(settings.map_stringPars.find("ReadinLocalMaxName")==settings.map_stringPars.end())
    settings.map_stringPars["ReadinLocalMaxName"] = "AllLocalMax";
  if(settings.map_stringPars.find("OutputLongiClusName")==settings.map_stringPars.end())
    settings.map_stringPars["OutputLongiClusName"] = "TruthTrackAxis";

  return StatusCode::SUCCESS;
};

StatusCode TruthTrackMatchingAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  m_TrackCol.clear();
  p_HalfClusterV = nullptr;
  p_HalfClusterU = nullptr;

  for(int itrk=0; itrk<m_datacol.TrackCol.size(); itrk++ ) m_TrackCol.push_back(m_datacol.TrackCol[itrk].get());
  p_HalfClusterU = &(m_datacol.map_HalfCluster["HalfClusterColU"]);
  p_HalfClusterV = &(m_datacol.map_HalfCluster["HalfClusterColV"]);

  return StatusCode::SUCCESS;
};

StatusCode TruthTrackMatchingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){
  if( m_TrackCol.size()==0) return StatusCode::SUCCESS;

  for(int itrk=0; itrk<m_TrackCol.size(); itrk++){

    //Get MCP linked to the track
    edm4hep::MCParticle mcp_trk = m_TrackCol[itrk]->getLeadingMCP();
 
    //Loop in HFCluster U
    for(int ihf=0; ihf<p_HalfClusterU->size(); ihf++){
      //Get Local max of the HalfCluster
      std::vector<const PandoraPlus::Calo1DCluster*> localMaxColU = p_HalfClusterU->at(ihf).get()->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

      //Loop for the localMax: 
      std::shared_ptr<PandoraPlus::CaloHalfCluster> t_axis = std::make_shared<PandoraPlus::CaloHalfCluster>();
      for(int ilm=0; ilm<localMaxColU.size(); ilm++){
        auto truthMap_lm = localMaxColU[ilm]->getLinkedMCP();
        edm4hep::MCParticle mcp_lm = localMaxColU[ilm]->getLeadingMCP(); 
        if(mcp_lm==mcp_trk)
          t_axis->addUnit(localMaxColU[ilm]);
        
      }

      // If the track does not match the Halfcluster, the track axis candidate will have no 1DCluster
      if(t_axis->getCluster().size()==0)
        continue;

      t_axis->addAssociatedTrack(m_TrackCol[itrk]);
      t_axis->setType(10000); //Track-type axis.
      m_TrackCol[itrk]->addAssociatedHalfClusterU( p_HalfClusterU->at(ihf).get() );
      m_datacol.map_HalfCluster["bkHalfCluster"].push_back(t_axis);
      p_HalfClusterU->at(ihf).get()->addHalfCluster(settings.map_stringPars["OutputLongiClusName"], t_axis.get());
    }

    //Loop in HFCluster V
    for(int ihf=0; ihf<p_HalfClusterV->size(); ihf++){
      //Get Local max of the HalfCluster
      std::vector<const PandoraPlus::Calo1DCluster*> localMaxColV = p_HalfClusterV->at(ihf).get()->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

      //Loop for the localMax:
      std::shared_ptr<PandoraPlus::CaloHalfCluster> t_axis = std::make_shared<PandoraPlus::CaloHalfCluster>();
      for(int ilm=0; ilm<localMaxColV.size(); ilm++){
        auto truthMap_lm = localMaxColV[ilm]->getLinkedMCP();
        edm4hep::MCParticle mcp_lm = localMaxColV[ilm]->getLeadingMCP();
        if(mcp_lm==mcp_trk)
          t_axis->addUnit(localMaxColV[ilm]);
        
      }

      // If the track does not match the Halfcluster, the track axis candidate will have no 1DCluster
      if(t_axis->getCluster().size()==0)
        continue;

      t_axis->addAssociatedTrack(m_TrackCol[itrk]);
      t_axis->setType(10000); //Track-type axis.
      m_TrackCol[itrk]->addAssociatedHalfClusterV( p_HalfClusterV->at(ihf).get() );
      m_datacol.map_HalfCluster["bkHalfCluster"].push_back(t_axis);
      p_HalfClusterV->at(ihf).get()->addHalfCluster(settings.map_stringPars["OutputLongiClusName"], t_axis.get());
    }

  }//End loop track


  //Loop track to check the associated cluster: merge clusters if they are associated to the same track.
  std::vector<PandoraPlus::CaloHalfCluster*> tmp_deleteClus; tmp_deleteClus.clear();
  for(auto &itrk : m_TrackCol){
    std::vector<PandoraPlus::CaloHalfCluster*> m_matchedUCol = itrk->getAssociatedHalfClustersU();
    std::vector<PandoraPlus::CaloHalfCluster*> m_matchedVCol = itrk->getAssociatedHalfClustersV();

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
  }//End loop track


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


cout<<"Print HalfClusters and axis"<<endl;
cout<<"  HalfClusterU size: "<<p_HalfClusterU->size()<<endl;
for(int i=0; i<p_HalfClusterU->size(); i++){
  std::vector<const PandoraPlus::Calo1DCluster*> m_1dcluster = p_HalfClusterU->at(i)->getCluster();
  std::vector<const PandoraPlus::Calo1DCluster*> m_localMax = p_HalfClusterU->at(i).get()->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);
  std::vector<const PandoraPlus::CaloHalfCluster*> m_axis = p_HalfClusterU->at(i)->getHalfClusterCol(settings.map_stringPars["OutputLongiClusName"]);
  printf("    In HalfCluster #%d: 1D cluster size %d, localMax size %d, axis size %d \n", i, m_1dcluster.size(), m_localMax.size(), m_axis.size());

  cout<<"    Print 1D cluster"<<endl;
  for(int i1d=0; i1d<m_1dcluster.size(); i1d++)
    printf("      ");
}



  return StatusCode::SUCCESS;
};

StatusCode TruthTrackMatchingAlg::ClearAlgorithm(){
  m_TrackCol.clear();
  p_HalfClusterV = nullptr;
  p_HalfClusterU = nullptr;

  return StatusCode::SUCCESS;
};



#endif
