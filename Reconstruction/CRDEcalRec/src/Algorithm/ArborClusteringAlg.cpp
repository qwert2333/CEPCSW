#ifndef _ARBORCLUSTERING_ALG_C
#define _ARBORCLUSTERING_ALG_C

#include "Algorithm/ArborClusteringAlg.h"

StatusCode ArborClusteringAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters

  return StatusCode::SUCCESS;
};

StatusCode ArborClusteringAlg::Initialize(){

  return StatusCode::SUCCESS;
};

StatusCode ArborClusteringAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  //Get readin info
  std::vector<PandoraPlus::CaloHit*> m_ECalHitCol; m_ECalHitCol.clear();
  std::vector<PandoraPlus::CaloHit*> m_HCalHitCol; m_HCalHitCol.clear();

  std::vector<PandoraPlus::Calo2DCluster*>* p_Tshowers = &(m_datacol.map_ShowerInLayer[settings.map_stringPars["ReadinHit"]]);
    //Convert transShower to hit
    for(int is=0; is<p_Tshowers->size(); is++){
      PandoraPlus::CaloHit* m_hit = new PandoraPlus::CaloHit();
      m_hit->setEnergy( p_Tshowers->at(is)->getEnergy() );
      m_hit->setPosition( p_Tshowers->at(is)->getPos() );
      m_hit->setLayer( p_Tshowers->at(is)->getDlayer() );
      m_hit->setType(1);
      //m_hit->setParentShower( p_Tshowers->at(is) );
      m_ECalHitCol.push_back(m_hit);
      m_datacol.bk_HitCol.push_back(m_hit);
    }
    p_Tshowers = nullptr;

  m_HCalHitCol = m_datacol.map_CaloHit["HCALBarrel"];


  std::vector<PandoraPlus::Track*>* p_trackCol = &(m_datacol.TrackCol);

  std::vector<PandoraPlus::Calo3DCluster*> m_clusterCol; m_clusterCol.clear(); 
  //Track driven clustering: find all CaloHits nearby a tracks.
  for(int itrk=0; itrk<p_trackCol->size(); itrk++){

    PandoraPlus::Calo3DCluster* m_clus = new PandoraPlus::Calo3DCluster();
    for(int ihit=0; ihit<m_HCalHitCol.size(); ihit++){
      if(m_HCalHitCol[ihit]->getLayer() > settings.map_intPars["TrkClusteringMaxLayer"]) continue; 
      double TrkHitDistance = GetTrkCaloHitDistance( p_trackCol->at(itrk), m_HCalHitCol[ihit] );
      if(TrkHitDistance<settings.map_floatPars["TrkClusteringRThresh"])
        m_clus->addHit( m_HCalHitCol[ihit] );
    }
    if(m_clus->getCaloHits().size()!=0){ 
      m_clus->addTrack(p_trackCol->at(itrk));
      m_clusterCol.push_back(m_clus);
    }
    m_datacol.bk_Cluster3DCol.push_back(m_clus);
  }
  p_trackCol = nullptr;

  //Convert hits to ArborNodes

  //Convert existed clusters to ArborTree

/*
  InitArborTree(m_orderedNodes, m_ArborTreeCol, m_isoNodes);

  std::vector<CRDEcalEDM::CRDArborTree> tmpTrees; tmpTrees.clear();
  MergeConnectedTrees( m_ArborTreeCol, tmpTrees );
  m_ArborTreeCol.clear(); m_ArborTreeCol = tmpTrees;

  //Clean connection
  for(int it=0; it<m_ArborTreeCol.size(); it++) CleanConnection(m_ArborTreeCol[it]);

  //Depart trees after connection cleaning.
  std::vector<CRDEcalEDM::CRDArborTree> m_departedTrees; m_departedTrees.clear();
  for(int it=0; it<m_ArborTreeCol.size(); it++){
    tmpTrees.clear();
    DepartArborTree(m_ArborTreeCol[it], tmpTrees, m_isoNodes);
    m_departedTrees.insert(m_departedTrees.end(), tmpTrees.begin(), tmpTrees.end());
    //m_ArborTreeCol.clear(); m_ArborTreeCol = tmpTrees;
  }
  m_ArborTreeCol.clear(); m_ArborTreeCol = m_departedTrees; m_departedTrees.clear();




  //Transform ArborTree to Clusters. 



  m_datacol.map_CaloCluster["ArborHCALCluster"] = 
*/
  return StatusCode::SUCCESS;
};

StatusCode ArborClusteringAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
};

double ArborClusteringAlg::GetTrkCaloHitDistance(const PandoraPlus::Track* trk, const PandoraPlus::CaloHit* hit){
  double dis = 10e5;
  for(int i=0; i<trk->trackStates_size(); i++){
    TVector3 pos_state = trk->getTrackStates(i).referencePoint; 
    double m_dis = (pos_state-hit->getPosition()).Mag();
    if(m_dis<dis) dis = m_dis;
  }

  return dis; 
};


StatusCode ArborClusteringAlg::InitArborTree( std::map<int, std::vector<ArborNode*> >& m_orderedNodes,
                                              std::vector<ArborTree>& m_treeCol,
                                              std::vector<ArborNode*>& m_isoNodes ){

  return StatusCode::SUCCESS;
}

StatusCode ArborClusteringAlg::MergeConnectedTrees( std::vector<ArborTree>& m_inTreeCol, std::vector<ArborTree>& m_outTreeCol ){

  return StatusCode::SUCCESS;
}

StatusCode ArborClusteringAlg::CleanConnection( ArborTree& m_tree ){

  return StatusCode::SUCCESS;
}

StatusCode ArborClusteringAlg::DepartArborTree( ArborTree& m_tree, std::vector<ArborTree>& m_departedTrees, std::vector<ArborNode*>& m_isoNodes ){

  return StatusCode::SUCCESS;
}

StatusCode ArborClusteringAlg::MergeNeighborTree( std::vector<ArborTree>& m_inTreeCol, std::vector<ArborTree>& m_outTreeCol ){

  return StatusCode::SUCCESS;
}


#endif
