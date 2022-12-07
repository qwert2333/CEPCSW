#ifndef _ARBORCLUSTERING_ALG_C
#define _ARBORCLUSTERING_ALG_C

#include "Algorithm/ArborClusteringAlg.h"

StatusCode ArborClusteringAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //if(settings.map_stringPars.find("ReadinHit")==settings.map_stringPars.end()) settings.map_stringPars["ReadinHit"] = "";
  if(settings.map_floatPars.find("Rth")==settings.map_floatPars.end()) settings.map_floatPars["Rth"] = 50;
  if(settings.map_floatPars.find("TrkClusteringMaxLayer")==settings.map_floatPars.end()) settings.map_floatPars["TrkClusteringMaxLayer"] = 5;
  if(settings.map_floatPars.find("TrkClusteringRThresh")==settings.map_floatPars.end()) settings.map_floatPars["TrkClusteringRThresh"] = 20.;

  //Initialize parameters

  return StatusCode::SUCCESS;
};

StatusCode ArborClusteringAlg::Initialize(){

  return StatusCode::SUCCESS;
};

StatusCode ArborClusteringAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  //Get readin info
  std::vector<PandoraPlus::CaloHit*> m_EcalHitCol; m_EcalHitCol.clear();
  std::vector<PandoraPlus::CaloHit*> m_HcalHitCol; m_HcalHitCol.clear();

  std::vector<PandoraPlus::Calo2DCluster*>* p_Tshowers = &(m_datacol.map_ShowerInLayer["EcalShowerInLayer"]);
    //Convert transShower to hit
    for(int is=0; is<p_Tshowers->size(); is++){
      PandoraPlus::CaloHit* m_hit = new PandoraPlus::CaloHit();
      m_hit->setEnergy( p_Tshowers->at(is)->getEnergy() );
      m_hit->setPosition( p_Tshowers->at(is)->getPos() );
      m_hit->setLayer( p_Tshowers->at(is)->getDlayer() );
      m_hit->setType(1);
      //m_hit->setParentShower( p_Tshowers->at(is) );
      m_EcalHitCol.push_back(m_hit);
      m_datacol.bk_HitCol.push_back(m_hit);
    }
    p_Tshowers = nullptr;

  m_HcalHitCol = m_datacol.map_CaloHit["HCALBarrel"];


  //Track driven clustering: find all CaloHits nearby a tracks.
  std::vector<PandoraPlus::Calo3DCluster*> m_clusterCol; m_clusterCol.clear(); 
  std::vector<PandoraPlus::Track*>* p_trackCol = &(m_datacol.TrackCol);
  for(int itrk=0; itrk<p_trackCol->size(); itrk++){

    PandoraPlus::Calo3DCluster* m_clus = new PandoraPlus::Calo3DCluster();
    for(int ihit=0; ihit<m_HcalHitCol.size(); ihit++){
      if(m_HcalHitCol[ihit]->getLayer() > settings.map_floatPars["TrkClusteringMaxLayer"]) continue; 
      double TrkHitDistance = GetTrkCaloHitDistance( p_trackCol->at(itrk), m_HcalHitCol[ihit] );
      if(TrkHitDistance<settings.map_floatPars["TrkClusteringRThresh"])
        m_clus->addHit( m_HcalHitCol[ihit] );
    }
    if(m_clus->getCaloHits().size()!=0){ 
      m_clus->addTrack(p_trackCol->at(itrk));
      m_clusterCol.push_back(m_clus);
    }
    m_datacol.bk_Cluster3DCol.push_back(m_clus);
  }
  p_trackCol = nullptr;
  m_datacol.map_CaloCluster["TrackHCALCluster"] = m_clusterCol; 
  
  //Convert hits to ArborNodes
  std::vector<PandoraPlus::ArborNode> m_EcalNode; m_EcalNode.clear();
  std::vector<PandoraPlus::ArborNode> m_HcalNode; m_HcalNode.clear();
  for(int ih=0; ih<m_EcalHitCol.size(); ih++){
    ArborNode m_node(m_EcalNode[ih]);
    m_EcalNode.push_back(m_node);
  }
  for(int ih=0; ih<m_HcalHitCol.size(); ih++){
    ArborNode m_node(m_HcalNode[ih]);
    m_HcalNode.push_back(m_node);
  }

  //Convert existed clusters to ArborTree
  std::vector<ArborTree> m_ArborTreeCol; m_ArborTreeCol.clear();
  for(int icl=0; icl<m_clusterCol.size(); icl++){
    ArborTree m_tree; m_tree.clear();
    std::vector<PandoraPlus::ArborNode*> tmp_nodes; tmp_nodes.clear();
    for(int ih=0; ih<m_clusterCol[icl]->getCaloHits().size(); ih++){
      ArborNode m_node(m_clusterCol[icl]->getCaloHits()[ih]);
      auto iter = find( m_HcalNode.begin(), m_HcalNode.end(), m_node );
      if(iter!=m_HcalNode.end()) tmp_nodes.push_back( &(*iter) );
    }

    if(tmp_nodes.size()<2) continue; 
    //sort(tmp_nodes.begin(), tmp_nodes.end(), compR);
    for(int inode=0; inode<tmp_nodes.size()-1; inode++){
      tmp_nodes[inode]->ConnectDaughter(tmp_nodes[inode+1]);
      ArborLink m_link = make_pair(tmp_nodes[inode], tmp_nodes[inode+1]);
      m_tree.push_back(m_link);
    }
  }

  //Initialize full links 
  std::vector<PandoraPlus::ArborNode> isoNodes; 
  //InitArborTree(m_HcalNode, m_ArborTreeCol, isoNodes );
  //InitArborTree(m_EcalNode, m_ArborTreeCol, isoNodes );

  //Link cleaning
  for(int itree=0; itree<m_ArborTreeCol.size(); itree++) CleanConnection( m_ArborTreeCol[itree] );


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


StatusCode ArborClusteringAlg::InitArborTree( std::vector<ArborNode>* p_nodeCol,
                                              std::vector<ArborTree>* p_treeCol,
                                              std::vector<ArborNode>* p_isoNodes ){
  //Build full connections within radius R. 
/*
  //Extand existing trees
  std::vector<ArborNode*> m_leftNodes; m_leftNodes.clear();
  for(int it=0; it<p_treeCol->size(); it++){
    ArborNode* m_nodeA = p_treeCol->at(it).first; 
    for(int in=0; in<p_nodeCol->size(); in++){
      TVector3 pos = p_nodeCol->at(in).getPosition();
      double R_AB = (pos-m_nodeA->getPosition()).Mag();
      if(pos.Mag()<m_nodeA->getPosition().Mag()) continue; //In front of node A. 
      if(R_AB>settings.map_floatPars["Rth"]) continue;

      ArborLink m_link = make_pair(m_nodeA, &(p_nodeCol->at(i)));
      auto iter = find(p_nodeCol->begin(), p_nodeCol->end(), m_link );
      if(iter==p_nodeCol->end()) p_treeCol->at(it).push_back(m_link);
      else{
        auto iter_node = find(m_leftNodes.begin(), m_leftNodes.end(), &(p_nodeCol->at(i)) );
        if( iter_node==m_leftNodes.end() ) m_leftNodes.push_back( &(p_nodeCol->at(i)) );
      }
    }
  }


  //Build new trees 
  for(int in=0; in<m_leftNodes.size(); in++){
  for(int jn=0; jn!=in && jn<m_leftNodes.size(); jn++){
    TVector3 posA = m_leftNodes[in]->getPosition();
    TVector3 posB = m_leftNodes[jn]->getPosition();
    if(posA.Mag()>posB.Mag()) continue;
    if( (posA-posB).Mag()>settings.map_floatPars["Rth"] ) continue;

    ArborLink m_link = make_pair( m_leftNodes[in], m_leftNodes[jn] );
    
  }

  }
*/


  return StatusCode::SUCCESS;
}

//StatusCode ArborClusteringAlg::MergeConnectedTrees( std::vector<ArborTree>& m_inTreeCol, std::vector<ArborTree>& m_outTreeCol ){

//  return StatusCode::SUCCESS;
//}

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
