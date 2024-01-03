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

cout<<"TruthTrackMatchingAlg: Input track size "<<m_TrackCol.size()<<", HFClusterU size "<<p_HalfClusterU->size()<<", HFClusterV size "<<p_HalfClusterV->size()<<endl;

  return StatusCode::SUCCESS;
};

StatusCode TruthTrackMatchingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){
  if( m_TrackCol.size()==0) return StatusCode::SUCCESS;

  for(int itrk=0; itrk<m_TrackCol.size(); itrk++){

    //Get MCP linked to the track
    edm4hep::MCParticle mcp_trk = m_TrackCol[itrk]->getLeadingMCP();
//cout<<"  Track #"<<itrk<<": MC pid "<<mcp_trk.getPDG()<<endl;
 
    //Loop in HFCluster U
    for(int ihf=0; ihf<p_HalfClusterU->size(); ihf++){
      //Get Local max of the HalfCluster
      std::vector<const PandoraPlus::Calo1DCluster*> localMaxColU = p_HalfClusterU->at(ihf).get()->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

      //Loop for the localMax: 
      std::shared_ptr<PandoraPlus::CaloHalfCluster> t_axis = std::make_shared<PandoraPlus::CaloHalfCluster>();
      for(int ilm=0; ilm<localMaxColU.size(); ilm++){
        edm4hep::MCParticle mcp_lm = localMaxColU[ilm]->getLeadingMCP(); 
        if(mcp_lm==mcp_trk)
          t_axis->addUnit(localMaxColU[ilm]);
        
      }

      // If the track does not match the Halfcluster, the track axis candidate will have no 1DCluster
      if(!t_axis || t_axis->getCluster().size()==0)
        continue;

//cout<<"    HFCluster #"<<ihf<<" match with track. New axis localMax size "<<t_axis->getCluster().size()<<endl;
      t_axis->getLinkedMCPfromUnit();
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
        edm4hep::MCParticle mcp_lm = localMaxColV[ilm]->getLeadingMCP();
        if(mcp_lm==mcp_trk)
          t_axis->addUnit(localMaxColV[ilm]);
        
      }

      // If the track does not match the Halfcluster, the track axis candidate will have no 1DCluster
      if(!t_axis || t_axis->getCluster().size()==0)
        continue;

//cout<<"    HFCluster #"<<ihf<<" match with track. New axis localMax size "<<t_axis->getCluster().size()<<endl;
      t_axis->getLinkedMCPfromUnit();
      t_axis->addAssociatedTrack(m_TrackCol[itrk]);
      t_axis->setType(10000); //Track-type axis.
      m_TrackCol[itrk]->addAssociatedHalfClusterV( p_HalfClusterV->at(ihf).get() );
      m_datacol.map_HalfCluster["bkHalfCluster"].push_back(t_axis);
      p_HalfClusterV->at(ihf).get()->addHalfCluster(settings.map_stringPars["OutputLongiClusName"], t_axis.get());
    }

  }//End loop track

//cout<<"Check HFClusters linked to the same track"<<endl;
  //Loop track to check the associated cluster: merge clusters if they are associated to the same track.
  std::vector<PandoraPlus::CaloHalfCluster*> tmp_deleteClus; tmp_deleteClus.clear();
  for(auto &itrk : m_TrackCol){
    std::vector<PandoraPlus::CaloHalfCluster*> m_matchedUCol = itrk->getAssociatedHalfClustersU();
    std::vector<PandoraPlus::CaloHalfCluster*> m_matchedVCol = itrk->getAssociatedHalfClustersV();
//cout<<"  In track "<<itrk<<": linked HFU size "<<m_matchedUCol.size()<<", linked HFV size "<<m_matchedVCol.size()<<endl;

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

  //Clean the duplicated axes
  for(int ihc=0; ihc<p_HalfClusterU->size(); ihc++){
    std::vector<const PandoraPlus::CaloHalfCluster*> m_axis = p_HalfClusterU->at(ihc)->getHalfClusterCol(settings.map_stringPars["OutputLongiClusName"]);
    std::vector<const PandoraPlus::CaloHalfCluster*> m_mergedAxis; m_mergedAxis.clear();

//cout<<"  In HFCluster #"<<ihc<<": track axis size "<<m_axis.size();
    if(m_axis.size()<=1){ 
//cout<<endl;
      continue;
    }

    std::map<const PandoraPlus::Track*, std::vector<const PandoraPlus::CaloHalfCluster*>> map_trk; map_trk.clear();
    for(int iax=0; iax<m_axis.size(); iax++){
      if(m_axis[iax]->getAssociatedTracks().size()==0){ 
        m_mergedAxis.push_back(m_axis[iax]);
        continue;
      }
      if(m_axis[iax]->getAssociatedTracks().size()>1)
        std::cout<<"Warning: track axis has "<<m_axis[iax]->getAssociatedTracks().size()<<" matched tracks! Need to check! "<<endl;

      map_trk[m_axis[iax]->getAssociatedTracks()[0]].push_back(m_axis[iax]);
    }
//cout<<"  map size "<<map_trk.size()<<", neutral axis size "<<m_mergedAxis.size()<<endl;;

    for(auto &itrk: map_trk){
//cout<<"  In track: axis size "<<itrk.second.size()<<endl;
      if(itrk.second.size()==0) continue;
      if(itrk.second.size()==1){ 
        m_mergedAxis.push_back(itrk.second[0]);
        continue;
      }
      PandoraPlus::CaloHalfCluster* p_axis = const_cast<PandoraPlus::CaloHalfCluster*>(itrk.second[0]);
      for(int iax=1; iax<itrk.second.size(); iax++) p_axis->mergeHalfCluster(itrk.second[iax]);
      m_mergedAxis.push_back(p_axis);
    }
    
//cout<<". After merge: axis size "<<m_mergedAxis.size()<<endl;
     p_HalfClusterU->at(ihc)->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], m_mergedAxis);
  }

  for(int ihc=0; ihc<p_HalfClusterV->size(); ihc++){
    std::vector<const PandoraPlus::CaloHalfCluster*> m_axis = p_HalfClusterV->at(ihc)->getHalfClusterCol(settings.map_stringPars["OutputLongiClusName"]);
    std::vector<const PandoraPlus::CaloHalfCluster*> m_mergedAxis; m_mergedAxis.clear();

//cout<<"  In HFCluster #"<<ihc<<": track axis size "<<m_axis.size();
    if(m_axis.size()<=1){
//cout<<endl;
      continue;
    }

    std::map<const PandoraPlus::Track*, std::vector<const PandoraPlus::CaloHalfCluster*>> map_trk; map_trk.clear();
    for(int iax=0; iax<m_axis.size(); iax++){
      if(m_axis[iax]->getAssociatedTracks().size()==0){
        m_mergedAxis.push_back(m_axis[iax]);
        continue;
      }
      if(m_axis[iax]->getAssociatedTracks().size()>1)
        std::cout<<"Warning: track axis has "<<m_axis[iax]->getAssociatedTracks().size()<<" matched tracks! Need to check! "<<endl;

      map_trk[m_axis[iax]->getAssociatedTracks()[0]].push_back(m_axis[iax]);
    }
//cout<<"  map size "<<map_trk.size()<<", neutral axis size "<<m_mergedAxis.size()<<endl;;

    for(auto &itrk: map_trk){
//cout<<"  In track: axis size "<<itrk.second.size()<<endl;
      if(itrk.second.size()==0) continue;
      if(itrk.second.size()==1){
        m_mergedAxis.push_back(itrk.second[0]);
        continue;
      }
      PandoraPlus::CaloHalfCluster* p_axis = const_cast<PandoraPlus::CaloHalfCluster*>(itrk.second[0]);
      for(int iax=1; iax<itrk.second.size(); iax++) p_axis->mergeHalfCluster(itrk.second[iax]);
      m_mergedAxis.push_back(p_axis);
    }

//cout<<". After merge: axis size "<<m_mergedAxis.size()<<endl;
     p_HalfClusterV->at(ihc)->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], m_mergedAxis);
  }



/*
cout<<"Track size: "<<m_TrackCol.size()<<endl;
cout<<"Print HalfClusters and axis"<<endl;
cout<<"  HalfClusterU size: "<<p_HalfClusterU->size()<<endl;
for(int i=0; i<p_HalfClusterU->size(); i++){
  std::vector<const PandoraPlus::Calo1DCluster*> m_1dcluster = p_HalfClusterU->at(i)->getCluster();
  std::vector<const PandoraPlus::Calo1DCluster*> m_localMax = p_HalfClusterU->at(i).get()->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);
  std::vector<const PandoraPlus::CaloHalfCluster*> m_axis = p_HalfClusterU->at(i)->getHalfClusterCol(settings.map_stringPars["OutputLongiClusName"]);
  printf("    In HalfCluster #%d: 1D cluster size %d, localMax size %d, axis size %d \n", i, m_1dcluster.size(), m_localMax.size(), m_axis.size());

  for(int iax=0; iax<m_axis.size(); iax++){
    printf("      In axis #%d: localMax size %d \n", iax, m_axis[iax]->getCluster().size() );
    for(int ilm=0; ilm<m_axis[iax]->getCluster().size(); ilm++){
      printf("        LocalMax #%d: layer %d, towersize %d, towerID [%d, %d, %d], En %.4f, leadingMC pid %d, address %p \n", 
        ilm, m_axis[iax]->getCluster()[ilm]->getDlayer(), m_axis[iax]->getCluster()[ilm]->getTowerID().size(), 
        m_axis[iax]->getCluster()[ilm]->getTowerID()[0][0], m_axis[iax]->getCluster()[ilm]->getTowerID()[0][1],m_axis[iax]->getCluster()[ilm]->getTowerID()[0][2],
        m_axis[iax]->getCluster()[ilm]->getEnergy(), m_axis[iax]->getCluster()[ilm]->getLeadingMCP().getPDG(), m_axis[iax]->getCluster()[ilm] );
    }
  }
}

cout<<"  HalfClusterV size: "<<p_HalfClusterV->size()<<endl;
for(int i=0; i<p_HalfClusterV->size(); i++){
  std::vector<const PandoraPlus::Calo1DCluster*> m_1dcluster = p_HalfClusterV->at(i)->getCluster();
  std::vector<const PandoraPlus::Calo1DCluster*> m_localMax = p_HalfClusterV->at(i).get()->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);
  std::vector<const PandoraPlus::CaloHalfCluster*> m_axis = p_HalfClusterV->at(i)->getHalfClusterCol(settings.map_stringPars["OutputLongiClusName"]);
  printf("    In HalfCluster #%d: 1D cluster size %d, localMax size %d, axis size %d \n", i, m_1dcluster.size(), m_localMax.size(), m_axis.size());

  for(int iax=0; iax<m_axis.size(); iax++){
    printf("      In axis #%d: localMax size %d \n", iax, m_axis[iax]->getCluster().size() );
    for(int ilm=0; ilm<m_axis[iax]->getCluster().size(); ilm++){
      printf("        LocalMax #%d: layer %d, towersize %d, towerID [%d, %d, %d], En %.4f, leadingMC pid %d, address %p \n",
        ilm, m_axis[iax]->getCluster()[ilm]->getDlayer(), m_axis[iax]->getCluster()[ilm]->getTowerID().size(),
        m_axis[iax]->getCluster()[ilm]->getTowerID()[0][0], m_axis[iax]->getCluster()[ilm]->getTowerID()[0][1],m_axis[iax]->getCluster()[ilm]->getTowerID()[0][2],
        m_axis[iax]->getCluster()[ilm]->getEnergy(), m_axis[iax]->getCluster()[ilm]->getLeadingMCP().getPDG(), m_axis[iax]->getCluster()[ilm] );
    }
  }
}
*/

  return StatusCode::SUCCESS;
};

StatusCode TruthTrackMatchingAlg::ClearAlgorithm(){
  m_TrackCol.clear();
  p_HalfClusterV = nullptr;
  p_HalfClusterU = nullptr;

  return StatusCode::SUCCESS;
};



#endif
