#ifndef _ARBORTREEMERGING_ALG_C
#define _ARBORTREEMERGING_ALG_C
#include "Algorithm/ArborTreeMergingAlg.h"
void ArborTreeMergingAlg::Settings::SetInitialValue(){

    th_Nsigma = 2; 
    th_daughterR = 40.; 
    th_GoodTreeLevel = 3;
    th_GoodTreeNodes = 10;
    th_MergeR = 100; 
    th_MergeTheta = PI/10.;
    fl_MergeTrees = true; 
}

StatusCode ArborTreeMergingAlg::Initialize(){

  return StatusCode::SUCCESS;
}

StatusCode ArborTreeMergingAlg::RunAlgorithm( ArborTreeMergingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){
  settings = m_settings; 

  std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_2DshowerCol = m_datacol.Shower2DCol;
  std::vector<CRDEcalEDM::CRDArborTree> m_ArborTreeCol = m_datacol.ArborTreeCol; 
  std::vector<CRDEcalEDM::CRDArborNode*> m_isoNodes = m_datacol.IsoNodes; 

  //Classify good and bad Tree.
  std::vector<CRDEcalEDM::CRDArborTree*> m_goodTreeCol; m_goodTreeCol.clear();
  std::vector<CRDEcalEDM::CRDArborTree*> m_badTreeCol; m_badTreeCol.clear();
  for(int it=0; it<m_ArborTreeCol.size(); it++){
    if( (m_ArborTreeCol[it].GetMaxDlayer()-m_ArborTreeCol[it].GetMinDlayer()) >= settings.th_GoodTreeLevel &&
         m_ArborTreeCol[it].GetNodes().size()>=settings.th_GoodTreeNodes)
      m_goodTreeCol.push_back( &m_ArborTreeCol[it] );
    else m_badTreeCol.push_back( &m_ArborTreeCol[it] );
  }

  if(settings.fl_MergeTrees){
  std::map<CRDEcalEDM::CRDArborNode*, CRDEcalEDM::CRDArborTree*> m_svtx; m_svtx.clear();
  //Type 1 secondary vertex: large variation in variance. 
  GetSecondaryVtxT1(m_ArborTreeCol, m_2DshowerCol, m_svtx);
  //Type 2 secondary vertex: branch node in a tree & large distance between daughter nodes.
  GetSecondaryVtxT2(m_ArborTreeCol, m_svtx);

  //Merge clusters pointing to the vertex. 
  std::map<CRDEcalEDM::CRDArborNode*, CRDEcalEDM::CRDArborTree*>::iterator iter = m_svtx.begin(); 
  for(iter; iter!=m_svtx.end(); iter++){

    //Get vertex
    CRDEcalEDM::CRDArborNode* m_node = iter->first; 
    int m_layer = m_node->GetDlayer(); 
    TVector3 m_vtxPos = m_node->GetPosition(); 

//printf("vtx: (%.2f, %.2f, %.2f) (Layer %d) \n", m_vtxPos.x(), m_vtxPos.y(), m_vtxPos.z(), m_layer);

    for(int it=0; it<m_goodTreeCol.size(); it++){
      if( m_goodTreeCol[it]->GetMinDlayer()<0      || m_goodTreeCol[it]->GetMaxDlayer()<0 ||  
          (m_layer>m_goodTreeCol[it]->GetMinDlayer() && m_layer<m_goodTreeCol[it]->GetMaxDlayer()) ) continue; 

      m_goodTreeCol[it]->SortNodes();  //Sort from outer layer to inner layer. 
      //std::vector<CRDEcalEDM::CRDArborNode*> m_nodesInTree; m_nodesInTree.clear(); 
      //m_nodesInTree = m_goodTreeCol[it]->GetNodes(); 

      TVector3 m_nodePos; 
      if( m_layer<m_goodTreeCol[it]->GetMinDlayer() )   m_nodePos = m_goodTreeCol[it]->GetNodes().back()->GetPosition(); 
      if( m_layer>m_goodTreeCol[it]->GetMaxDlayer() )   m_nodePos = m_goodTreeCol[it]->GetNodes().front()->GetPosition(); 

      CRDEcalEDM::CRDCaloHit3DCluster m_clus = m_goodTreeCol[it]->ConvertTreeToCluster(); 
      m_clus.FitAxis(); 
      TVector3 m_treeAxis = m_clus.getAxis(); 

//printf("  Root/leaf node: (%.2f, %.2f, %.2f) \t", m_nodePos.x(), m_nodePos.y(), m_nodePos.z());
//printf("  cluster axis: (%.2f, %.2f, %.2f) \n", m_treeAxis.x(), m_treeAxis.y(), m_treeAxis.z()); 

      //Merge clusters
      if( (m_vtxPos-m_nodePos).Mag()<settings.th_MergeR && 
          ( m_treeAxis.Angle(m_vtxPos-m_nodePos)<settings.th_MergeTheta || (PI-m_treeAxis.Angle(m_vtxPos-m_nodePos))<settings.th_MergeTheta ) ){
        iter->second->AddNode( m_goodTreeCol[it]->GetNodes() );
        m_goodTreeCol.erase( m_goodTreeCol.begin()+it );
        it--; 
      }
    }
  }
  }

  //Save Trees into CRDCaloHit3DCluster
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_goodClusterCol;  m_goodClusterCol.clear();
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_badClusterCol;  m_badClusterCol.clear();
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_ClusterCol;  m_ClusterCol.clear();
  for(int it=0; it<m_goodTreeCol.size(); it++)
    m_goodClusterCol.push_back( m_goodTreeCol[it]->ConvertTreeToCluster() );
  for(int it=0; it<m_badTreeCol.size(); it++)
    m_badClusterCol.push_back( m_badTreeCol[it]->ConvertTreeToCluster() );

  for(int in=0; in<m_isoNodes.size(); in++){
    CRDEcalEDM::CRDCaloHit3DCluster clus_isonodes;
    CRDEcalEDM::CRDCaloHit2DShower m_shower = m_isoNodes[in]->GetOriginShower();
    clus_isonodes.AddShower( m_shower );
    m_badClusterCol.push_back( clus_isonodes );
  }

  for(int it=0; it<m_ArborTreeCol.size(); it++)
    m_ClusterCol.push_back( m_ArborTreeCol[it].ConvertTreeToCluster() );

  m_datacol.GoodClus3DCol = m_goodClusterCol;
  m_datacol.BadClus3DCol = m_badClusterCol;
  m_datacol.Clus3DCol = m_ClusterCol;

//  std::cout<<"Print Good Tree: "<<std::endl;
//  for(int it=0; it<m_goodTreeCol.size(); it++) m_goodTreeCol[it]->PrintTree();
//  std::cout<<"Print Bad Tree: "<<std::endl;
//  for(int it=0; it<m_badTreeCol.size(); it++) m_badTreeCol[it]->PrintTree();
//  std::cout<<"Print Isolated Nodes: "<<std::endl;
//  for(int i=0; i<m_isoNodes.size(); i++) 
//    printf("(%.2f, %.2f, %.2f, %d) \n", m_isoNodes[i]->GetPosition().x(), m_isoNodes[i]->GetPosition().y(), m_isoNodes[i]->GetPosition().z(), m_isoNodes[i]->GetType());

  //Clear nodes pointer. 
  for(int i=0; i<m_goodTreeCol.size(); i++) { m_goodTreeCol[i]->Clear(); }
  for(int i=0; i<m_badTreeCol.size(); i++) { m_badTreeCol[i]->Clear(); }
  for(int i=0; i<m_isoNodes.size(); i++)     { delete m_isoNodes[i]; }

  return StatusCode::SUCCESS;
}

StatusCode ArborTreeMergingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}

StatusCode ArborTreeMergingAlg::GetSecondaryVtxT1(  std::vector<CRDEcalEDM::CRDArborTree>& m_ArborTreeCol,  
                                                    std::vector<CRDEcalEDM::CRDCaloHit2DShower>& m_2DshowerCol,
                                                    std::map<CRDEcalEDM::CRDArborNode*, CRDEcalEDM::CRDArborTree*>& m_svtx )
{

  if( m_ArborTreeCol.size()==0 || m_2DshowerCol.size()==0 )  return StatusCode::SUCCESS;

  //Get the variance in each layer:
  std::vector<float> m_scndM; m_scndM.clear(); m_scndM.resize(14);
  std::map<int, std::vector<CRDEcalEDM::CRDCaloHit2DShower> > m_orderedShower; m_orderedShower.clear();
  for(int is=0;is<m_2DshowerCol.size();is++){
    m_orderedShower[m_2DshowerCol[is].getDlayer()].push_back(m_2DshowerCol[is]);
  }
  for(int il=0; il<14; il++){
    std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_showers = m_orderedShower[il+1];
    if(m_showers.size()==0){ m_scndM[il]=0; continue; }

    double sumE=0;
    TVector3 cent(0,0,0);
    for(int is=0; is<m_showers.size(); is++){
      sumE+=m_showers[is].getShowerE();
      TVector3 pos(m_showers[is].getPos().x(), m_showers[is].getPos().y(), m_showers[is].getPos().z());
      cent += pos;
    }
    cent = cent*(1./(double)m_showers.size());

    double scndM=0;
    for(int is=0; is<m_showers.size(); is++){
      TVector3 pos(m_showers[is].getPos().x(), m_showers[is].getPos().y(), m_showers[is].getPos().z());
      double rpos = (pos-cent).Mag2();
      double wi = m_showers[is].getShowerE()/sumE;
      scndM += wi*rpos;
    }
    m_scndM[il] = scndM;
//std::cout<<"  Layer: "<<il<<"  variance: "<<m_scndM[il]<<std::endl;
  }
//std::cout<<std::endl; 

  //Get secondary vertex layer and position
  vector<int> Layer_svtx; Layer_svtx.clear();
  for(int il=1; il<14; il++)
    if( m_scndM[il-1]!=0 && (fabs(m_scndM[il]-m_scndM[il-1]) > settings.th_Nsigma * m_scndM[il-1]) ){ Layer_svtx.push_back(il-1); }

//std::cout<<std::endl; 
//std::cout<<"  Layer with svtx: "; 
//for(int i=0; i<Layer_svtx.size(); i++) std::cout<<Layer_svtx[i]<<'\t'; 
//std::cout<<std::endl; 


  for(int iv=0; iv<Layer_svtx.size(); iv++){
    for(int it=0; it<m_ArborTreeCol.size(); it++){
      if( Layer_svtx[iv]>m_ArborTreeCol[it].GetMaxDlayer() || Layer_svtx[iv]<m_ArborTreeCol[it].GetMinDlayer() ||
          m_ArborTreeCol[it].GetMinDlayer()==-1 || m_ArborTreeCol[it].GetMaxDlayer()==-1 ) continue;

      std::vector<CRDEcalEDM::CRDArborNode*> m_nodes = m_ArborTreeCol[it].GetNodes(Layer_svtx[iv]);
      for(int in=0; in<m_nodes.size(); in++){
//printf("  vtx node: (%.2f, %.2f, %.2f) \n", m_nodes[in]->GetPosition().x(), m_nodes[in]->GetPosition().y(), m_nodes[in]->GetPosition().z());
        m_svtx[m_nodes[in]] = &m_ArborTreeCol[it];
      }
    }
  }

  return StatusCode::SUCCESS;
}


StatusCode ArborTreeMergingAlg::GetSecondaryVtxT2(  std::vector<CRDEcalEDM::CRDArborTree>& m_ArborTreeCol,
                                                    std::map<CRDEcalEDM::CRDArborNode*, CRDEcalEDM::CRDArborTree*>& m_svtx )
{

  if(m_ArborTreeCol.size()==0) return StatusCode::SUCCESS;

  //Get first branch node
  for(int it=0; it<m_ArborTreeCol.size(); it++){
    m_ArborTreeCol[it].SortNodes(); 
    std::vector<CRDEcalEDM::CRDArborNode*> m_nodes = m_ArborTreeCol[it].GetNodes(); 

    std::vector<CRDEcalEDM::CRDArborNode*> m_BrNodes; m_BrNodes.clear(); 
    for(int in=0; in<m_nodes.size(); in++){
      int nodeType = m_nodes[in]->GetType(); 
      if(nodeType!=3) continue; 

      double maxR_daughter = -999;
      std::vector<CRDEcalEDM::CRDArborNode*> m_dtrNodes = m_nodes[in]->GetDaughterNodes(); 
      for(int id=0; id<m_dtrNodes.size(); id++){
      for(int jd=id+1; jd<m_dtrNodes.size(); jd++){
        double m_Dis = (m_dtrNodes[id]->GetPosition() - m_dtrNodes[jd]->GetPosition()).Mag(); 
        if( m_Dis>maxR_daughter ) maxR_daughter = m_Dis; 
      }}
      if(maxR_daughter<0){ std::cout<<"Error in GetSecondaryVtxT2: failed to get daughter node distance! Check it!"<<std::endl; continue; }

      if( nodeType==3 && maxR_daughter>settings.th_daughterR ) m_BrNodes.push_back(m_nodes[in]); 
    }

    for(int in=0; in<m_BrNodes.size(); in++){
      std::map<CRDEcalEDM::CRDArborNode*, CRDEcalEDM::CRDArborTree*>::iterator it_find = m_svtx.find(m_BrNodes[in]);
      if( it_find==m_svtx.end() ) m_svtx[m_BrNodes[in]] = &m_ArborTreeCol[it]; 
    }

  }

  return StatusCode::SUCCESS;
}

#endif
