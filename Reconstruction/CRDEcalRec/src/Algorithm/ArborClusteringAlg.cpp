#ifndef _ARBORCLUSTERING_ALG_C
#define _ARBORCLUSTERING_ALG_C

#include "Algorithm/ArborClusteringAlg.h"
#include <map>

void ArborClusteringAlg::Settings::SetInitialValue(){

    clusType="";

    Lth_start = 2; 
    Rth_start = 30; 
    Rth_value = 40; 
    Rth_slope = 0.5; 

    wiB = 1;
    wiF = 1;
    pTheta = 1;
    pR = 1;
    pE = 1; 
    Rth_nbrRoot = 40; // Root node distance threshold when merging neighbor trees.
    Debug = 0; 
}

void ArborClusteringAlg::Settings::PrintSettings() const{
  std::cout<<"  ArborClusteringAlg Settings: "<<std::endl;
  std::cout<<"Rth_value: "<<'\t'<<Rth_value<<std::endl; 
  std::cout<<"Rth_slope: "<<'\t'<<Rth_slope<<std::endl; 
  std::cout<<"wiB: "<<'\t'<<wiB<<std::endl; 
  std::cout<<"wiF: "<<'\t'<<wiF<<std::endl; 
  std::cout<<"pTheta: "<<'\t'<<pTheta<<std::endl; 
  std::cout<<"pR: "<<'\t'<<pR<<std::endl; 
  std::cout<<"pE: "<<'\t'<<pE<<std::endl; 
  std::cout<<"Rth_nbrRoot: "<<'\t'<<Rth_nbrRoot<<std::endl; 
  std::cout<<"Debug: "<<'\t'<<Debug<<std::endl; 

}

StatusCode ArborClusteringAlg::Initialize(){

  return StatusCode::SUCCESS;
}

StatusCode ArborClusteringAlg::RunAlgorithm( ArborClusteringAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol ){
  settings = m_settings;

  std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_2DshowerCol; m_2DshowerCol.clear();
  if( settings.clusType=="MIP" )     m_2DshowerCol = m_datacol.MIPShower2DCol;
  else if( settings.clusType=="EM" ) m_2DshowerCol = m_datacol.EMShower2DCol;
  else m_2DshowerCol = m_datacol.Shower2DCol;


  if(m_2DshowerCol.size()==0){ 
    std::cout<<"Warning: Empty input in ArborClusteringAlg. Please check previous algorithm!"<<endl;  
    m_datacol.ClearArbor(); 
    return StatusCode::SUCCESS; 
  }
//m_datacol.PrintShower(); 

  std::map<int, std::vector<CRDEcalEDM::CRDArborNode*> > m_orderedNodes;  m_orderedNodes.clear(); //map<layer, showers>

  //Transform showers to nodes, order the nodes by layer. 
  for(int is=0; is<m_2DshowerCol.size(); is++ ){
    int m_dlayer = m_2DshowerCol[is].getDlayer();
    //CRDEcalEDM::CRDArborNode m_node(m_2DshowerCol[is]);
    CRDEcalEDM::CRDArborNode* m_node = new CRDEcalEDM::CRDArborNode(m_2DshowerCol[is]);
    m_orderedNodes[m_dlayer].push_back(m_node);
  }

  int NtotNodes=0; 
  std::map<int, std::vector<CRDEcalEDM::CRDArborNode*> >::iterator iter_count = m_orderedNodes.begin(); 
  for(iter_count; iter_count!=m_orderedNodes.end(); iter_count++){
    NtotNodes += iter_count->second.size(); 
  }
  if(settings.Debug>=1) std::cout<<"Total node size: "<<NtotNodes<<std::endl; 

/*
std::cout<<"Print ordered nodes: "<<std::endl;
for(auto iter=m_orderedNodes.begin(); iter!=m_orderedNodes.end(); iter++){
  std::cout<<"#Layer: "<<iter->first<<'\t';
  for(int i=0; i<iter->second.size(); i++) printf("(%.2f, %.2f, %.2f, %d) \t", iter->second[i].GetPosition().x(), iter->second[i].GetPosition().y(), iter->second[i].GetPosition().z(), iter->second[i].GetType());
  std::cout<<std::endl;
}
*/
  std::vector<CRDEcalEDM::CRDArborTree> m_ArborTreeCol; m_ArborTreeCol.clear();
  std::vector<CRDEcalEDM::CRDArborNode*> m_isoNodes; m_isoNodes.clear();
  

  //Build tree
  InitArborTree(m_orderedNodes, m_ArborTreeCol, m_isoNodes);

  if(settings.Debug>=1) std::cout<<"Initial Ntrees = "<<m_ArborTreeCol.size()<<std::endl;
  if(settings.Debug>=2){
    for(int it=0; it<m_ArborTreeCol.size(); it++) m_ArborTreeCol[it].PrintTree();
    std::cout<<"Print Isolated Nodes: "<<std::endl;
    for(int i=0; i<m_isoNodes.size(); i++)
      printf("(%.2f, %.2f, %.2f, %d) \n", m_isoNodes[i]->GetPosition().x(), m_isoNodes[i]->GetPosition().y(), m_isoNodes[i]->GetPosition().z(), m_isoNodes[i]->GetType());
    std::cout<<std::endl;
  }  


  std::vector<CRDEcalEDM::CRDArborTree> tmpTrees; tmpTrees.clear();
  MergeConnectedTrees( m_ArborTreeCol, tmpTrees );
  m_ArborTreeCol.clear(); m_ArborTreeCol = tmpTrees;

  if(settings.Debug>=1) std::cout<<"After ConnectedTree merging: Ntree = "<<m_ArborTreeCol.size()<<std::endl;
  if(settings.Debug>=2){
    for(int it=0; it<m_ArborTreeCol.size(); it++) m_ArborTreeCol[it].PrintTree();
    std::cout<<"Print Isolated Nodes: "<<std::endl;
    for(int i=0; i<m_isoNodes.size(); i++)
      printf("(%.2f, %.2f, %.2f, %d) \n", m_isoNodes[i]->GetPosition().x(), m_isoNodes[i]->GetPosition().y(), m_isoNodes[i]->GetPosition().z(), m_isoNodes[i]->GetType());
    std::cout<<std::endl;
  }

  //Clean connection
  for(int it=0; it<m_ArborTreeCol.size(); it++) CleanConnection(m_ArborTreeCol[it]);

  if(settings.Debug>=1) std::cout<<"After connection cleaning: Ntree = "<<m_ArborTreeCol.size()<<std::endl;
  if(settings.Debug>=2){ 
    for(int it=0; it<m_ArborTreeCol.size(); it++) m_ArborTreeCol[it].PrintTree();
    std::cout<<"Print Isolated Nodes: "<<std::endl;
    for(int i=0; i<m_isoNodes.size(); i++)
      printf("(%.2f, %.2f, %.2f, %d) \n", m_isoNodes[i]->GetPosition().x(), m_isoNodes[i]->GetPosition().y(), m_isoNodes[i]->GetPosition().z(), m_isoNodes[i]->GetType());
    std::cout<<std::endl;
  }

  //Depart trees after connection cleaning. 
  std::vector<CRDEcalEDM::CRDArborTree> m_departedTrees; m_departedTrees.clear(); 
  for(int it=0; it<m_ArborTreeCol.size(); it++){ 
    tmpTrees.clear(); 
    DepartArborTree(m_ArborTreeCol[it], tmpTrees, m_isoNodes);
    m_departedTrees.insert(m_departedTrees.end(), tmpTrees.begin(), tmpTrees.end());
    //m_ArborTreeCol.clear(); m_ArborTreeCol = tmpTrees; 
  }
  m_ArborTreeCol.clear(); m_ArborTreeCol = m_departedTrees; m_departedTrees.clear(); 

  if(settings.Debug>=1) std::cout<<"After DepartArborTree: Ntree = "<<m_ArborTreeCol.size()<<std::endl;
  if(settings.Debug>=2){ 
    for(int it=0; it<m_ArborTreeCol.size(); it++) m_ArborTreeCol[it].PrintTree();
    std::cout<<"Print Isolated Nodes: "<<std::endl;
    for(int i=0; i<m_isoNodes.size(); i++)
      printf("(%.2f, %.2f, %.2f, %d) \n", m_isoNodes[i]->GetPosition().x(), m_isoNodes[i]->GetPosition().y(), m_isoNodes[i]->GetPosition().z(), m_isoNodes[i]->GetType());
    std::cout<<std::endl;
  }

  //Merge Neighbor tree (roots are close within the same layer)
  //tmpTrees.clear();
  //MergeNeighborTree(m_ArborTreeCol, tmpTrees);
  //m_ArborTreeCol.clear(); m_ArborTreeCol = tmpTrees;

  m_datacol.ArborTreeCol = m_ArborTreeCol; 
  m_datacol.IsoNodes = m_isoNodes; 

  int Nnodes_end = 0; 
  for(int it=0; it<m_ArborTreeCol.size(); it++) Nnodes_end += m_ArborTreeCol[it].GetNodes().size(); 
  Nnodes_end += m_isoNodes.size(); 
  if(settings.Debug>=1) std::cout<<"At the end: Node size = "<<NtotNodes<<std::endl;
  if(Nnodes_end != NtotNodes) 
    std::cout<<"WARNING!  ArborClusteringAlg: Initial node size("<<NtotNodes<<") is not equal to final node size("<<Nnodes_end<<")! May have memory leakage! "<<std::endl; 

  return StatusCode::SUCCESS;
}

StatusCode ArborClusteringAlg::ClearAlgorithm() {

  return StatusCode::SUCCESS;
}


StatusCode ArborClusteringAlg::InitArborTree( std::map<int, std::vector<CRDEcalEDM::CRDArborNode*> >& m_orderedNodes, 
                                              std::vector<CRDEcalEDM::CRDArborTree>& m_treeCol,
                                              std::vector<CRDEcalEDM::CRDArborNode*>& m_isoNodes )
{

//std::cout<<"Ordered node collection size: "<<m_orderedNodes.size()<<std::endl;
//std::cout<<"  Layer(Nnode): ";
//for(auto iter = m_orderedNodes.begin(); iter!=m_orderedNodes.end(); iter++) std::cout<<iter->first<<"("<<iter->second.size()<<")  ";
//std::cout<<std::endl;

  if( m_orderedNodes.size()==0 ) return StatusCode::SUCCESS;
  m_treeCol.clear();
  m_isoNodes.clear();

  std::map<int, std::vector<CRDEcalEDM::CRDArborNode*> >::iterator iter = m_orderedNodes.begin();

  //Nodes in First layer: initialize the trees. 
  for(int is=0; is<iter->second.size(); is++){
//std::cout<<"#Layer: "<<iter->first<<", Nnode: "<<iter->second.size()<<endl;
    int m_dlayer = iter->first; 
    double Rth; 
    if(m_dlayer<=settings.Lth_start) Rth = settings.Rth_start; 
    else Rth = settings.Rth_value + (double)m_dlayer * settings.Rth_slope; 

    CRDEcalEDM::CRDArborTree m_tree; 
    std::vector<CRDEcalEDM::CRDArborNode*>* m_nextLayer = &m_orderedNodes[iter->first+1]; 
    for(int js=0; js<m_nextLayer->size(); js++)
      if( (iter->second[is]->GetPosition()-m_nextLayer->at(js)->GetPosition()).Mag()<Rth ) 
        iter->second[is]->ConnectDaughter( (m_nextLayer->at(js)) ); 

    m_tree.AddNode( (iter->second[is]) );
    m_tree.AddNode( iter->second[is]->GetDaughterNodes() );
    m_treeCol.push_back( m_tree );
  }
  iter++;
//std::cout<<"End tree seeding. Tree number: "<< m_treeCol.size()<<", Next layer: "<<iter->first<<std::endl;
//for(int it=0; it<m_treeCol.size(); it++) m_treeCol[it].PrintTree();

  //In later layers
  for(iter; iter!=m_orderedNodes.end(); iter++){

//std::cout<<"  #Layer: "<<iter->first<<", Nnode: "<<iter->second.size()<<", Present Ntrees: "<<m_treeCol.size()<<std::endl;
//for(int it=0; it<m_treeCol.size(); it++) m_treeCol[it].PrintTree();

    int m_dlayer = iter->first;
    double Rth;
    if(m_dlayer<=settings.Lth_start) Rth = settings.Rth_start;
    else Rth = settings.Rth_value + (double)m_dlayer * settings.Rth_slope;

    for(int in=0; in<iter->second.size(); in++){
      CRDEcalEDM::CRDArborNode* m_node = iter->second[in];
//printf("    #Node: (%.2f, %.2f, %.2f) \n", m_node->GetPosition().x(),  m_node->GetPosition().y(),  m_node->GetPosition().z());
      //Loop all trees:
      for(int it=0; it<m_treeCol.size(); it++){
//std::cout<<"    Loop in tree: "<<std::endl;
//m_treeCol[it].PrintTree();
//std::cout<<"    Node in this tree: "<<m_node->isInTree(m_treeCol[it])<<std::endl; 
        //If node in this tree: 
        if( m_node->isInTree(m_treeCol[it]) ){          
          //Get daughter nodes and add them in tree
          std::vector<CRDEcalEDM::CRDArborNode*>* m_nextLayer = &m_orderedNodes[m_dlayer+1];
          for(int jn=0; jn<m_nextLayer->size(); jn++)
            if( (m_node->GetPosition()-m_nextLayer->at(jn)->GetPosition()).Mag()<Rth ) 
              m_node->ConnectDaughter( (m_nextLayer->at(jn)) );
          m_treeCol[it].AddNode( m_node->GetDaughterNodes() );
          break;
        }
      }//end loop trees
//std::cout<<"End loop nodes in this layer! "<<std::endl; 
//std::cout<<std::endl;
    } //end loop nodes in this layer

    //if collection is not empty, create a new tree.

    for(int in=0; in<iter->second.size(); in++){
      CRDEcalEDM::CRDArborNode* m_node = iter->second[in];
      if( m_node->GetDaughterNodes().size()==0 && m_node->GetParentNodes().size()==0){
        CRDEcalEDM::CRDArborTree m_tree;
        m_tree.AddNode( m_node );
        std::vector<CRDEcalEDM::CRDArborNode*>* m_nextLayer = &m_orderedNodes[m_dlayer+1];
        for(int jn=0; jn<m_nextLayer->size(); jn++)
          if( (m_node->GetPosition() - m_nextLayer->at(jn)->GetPosition()).Mag()<Rth ) 
            m_node->ConnectDaughter( (m_nextLayer->at(jn)) );        
        m_tree.AddNode( m_node->GetDaughterNodes() );
        m_treeCol.push_back( m_tree );
      }
    } //end new tree.

  }//end loop all nodes.

//std::cout<<"End loop all nodes! Print all trees: "<<std::endl;
//for(int it=0; it<m_treeCol.size(); it++) m_treeCol[it].PrintTree();

  int NtotNodes=0;
//  for(int it=0; it<m_treeCol.size(); it++){
//    NtotNodes += m_treeCol[it].GetNodes().size();
//  }
//  std::cout<<"Total nodes before picking out isoNodes: "<<NtotNodes<<std::endl;
//  std::cout<<"Print all trees: "<<std::endl;
//  for(int it=0; it<m_treeCol.size(); it++) m_treeCol[it].PrintTree();

  //Loop treeCol, pickout isoNodes, update the orderedNodes.
  for(int it=0; it<m_treeCol.size(); it++){
    if(m_treeCol[it].GetNodes().size()==1){ 
      CRDEcalEDM::CRDArborNode* m_node = m_treeCol[it].GetNodes()[0];
      m_node->SetType(0);
      m_isoNodes.push_back( m_node );
      m_treeCol.erase( m_treeCol.begin()+it );
      it--;
    }
  }

  NtotNodes=0;
  for(int it=0; it<m_treeCol.size(); it++){
    NtotNodes += m_treeCol[it].GetNodes().size();
  }
//std::cout<<"Final tree size: "<<m_treeCol.size(); 
//std::cout<<"  Total nodes in tree: "<<NtotNodes; 
//std::cout<<"  Isolated note size: "<< m_isoNodes.size()<<std::endl;

  return StatusCode::SUCCESS;
}


StatusCode ArborClusteringAlg::MergeConnectedTrees( std::vector<CRDEcalEDM::CRDArborTree>& m_inTreeCol,
                                                    std::vector<CRDEcalEDM::CRDArborTree>& m_outTreeCol )
{
  m_outTreeCol.clear(); 
  if(m_inTreeCol.size()==0){ return StatusCode::SUCCESS;}
  if(m_inTreeCol.size()==1){ m_outTreeCol=m_inTreeCol; return StatusCode::SUCCESS;}

  //Create a tree with topological link, with one node as seed. 
  for(int it=0; it<m_inTreeCol.size(); it++){
    m_inTreeCol[it].SortNodes(); 
    std::vector<CRDEcalEDM::CRDArborNode*> m_nodeCol = m_inTreeCol[it].GetNodes(); 
    if(m_nodeCol.size()==0) continue; 

    CRDEcalEDM::CRDArborTree m_tree; 
    m_tree.CreateWithTopo( m_nodeCol[0] );
    m_outTreeCol.push_back(m_tree);
    break;
  }

  //Check if other nodes in this tree. If not, create a new tree. 
  for(int it=0; it<m_inTreeCol.size(); it++){
    std::vector<CRDEcalEDM::CRDArborNode*> m_nodeCol = m_inTreeCol[it].GetNodes();
    for(int in=0; in<m_nodeCol.size(); in++){
      bool f_inTree = false; 
      for(int jt=0; jt<m_outTreeCol.size(); jt++)
        if( m_nodeCol[in]->isInTree(m_outTreeCol[jt]) ) {f_inTree=true; break;}
      if(!f_inTree){
        CRDEcalEDM::CRDArborTree m_tree;
        m_tree.CreateWithTopo( m_nodeCol[in] );
        m_outTreeCol.push_back(m_tree);
      }      
    }

  }
//std::cout<<"Ntrees before: "<<m_inTreeCol.size()<<", Ntrees after: "<<m_outTreeCol.size()<<std::endl;

  return StatusCode::SUCCESS;
}


StatusCode ArborClusteringAlg::CleanConnection( CRDEcalEDM::CRDArborTree& m_tree ){

  m_tree.SortNodes(); 
//std::cout<<"First node: Layer="<<m_tree.GetNodes()[0]->GetDlayer()<<" Type="<<m_tree.GetNodes()[0]->GetType()<<std::endl;

  std::vector<CRDEcalEDM::CRDArborNode*> m_nodeCol = m_tree.GetNodes();
  std::vector< std::pair<CRDEcalEDM::CRDArborNode*, CRDEcalEDM::CRDArborNode*> > m_connectorCol; m_connectorCol.clear();  
  for(int in=0; in<m_nodeCol.size(); in++){
    if(in==0 && m_nodeCol[in]->GetDaughterNodes().size()!=0 ){ std::cout<<"WARNING in ArborClusteringAlg: Last node still have daughters! Check it!"<<std::endl; continue; }

//std::cout<<"#Layer: "<<m_nodeCol[in]->GetDlayer()<<" Parent size:"<<m_nodeCol[in]->GetParentNodes().size()<<", Daughter size: "<<m_nodeCol[in]->GetDaughterNodes().size()<<std::endl;

    //Get the connection with minimum kappa order
    TVector3 Cref = m_nodeCol[in]->GetRefDir(settings.wiF, settings.wiB);
    double kappaMin = 9999.;
    double deltaMax = -999.; 
    double deltaEmax = -999.; 
    CRDEcalEDM::CRDArborNode* m_foundnode = nullptr;
    std::vector<CRDEcalEDM::CRDArborNode*> m_parentNodes = m_nodeCol[in]->GetParentNodes();
    if(m_parentNodes.size()<=1) continue; 

    for(int ip=0; ip<m_parentNodes.size(); ip++){
      TVector3 relP = m_nodeCol[in]->GetPosition() - m_parentNodes[ip]->GetPosition();
      if(relP.Mag()>deltaMax) deltaMax = relP.Mag();
      double deltaE = fabs( m_nodeCol[in]->GetEnergy()-m_parentNodes[ip]->GetEnergy() ); 
      if(deltaE>deltaEmax) deltaEmax = deltaE; 
    }

    for(int ip=0; ip<m_parentNodes.size(); ip++){
      TVector3 relP = m_nodeCol[in]->GetPosition() - m_parentNodes[ip]->GetPosition();
      double theta  = relP.Angle(Cref)/(PI) ;
      double deltaR = relP.Mag();
      double deltaE = fabs( m_nodeCol[in]->GetEnergy()-m_parentNodes[ip]->GetEnergy() ); 

      double kappa = pow(theta, settings.pTheta) * pow(deltaR/deltaMax, settings.pR ) * pow(deltaE/deltaEmax, settings.pE ); 
//std::cout<<"   Node kappa order: "<<theta<<"  "<<deltaR<<"  "<<kappa<<std::endl;
      if(kappa<kappaMin) { kappaMin=kappa; m_foundnode=m_parentNodes[ip]; }
    }
//std::cout<<"  Min Kappa order: "<<kappaMin<<std::endl;
    if( m_foundnode==nullptr ){ std::cout<<"WARNING in ArborClusteringAlg: Did not find minimum kappa order!"<<std::endl; continue; }

    //Save two nodes into connector collection
    for(int ip=0; ip<m_parentNodes.size(); ip++){
      if(m_parentNodes[ip]==m_foundnode) continue; 
      std::pair<CRDEcalEDM::CRDArborNode*, CRDEcalEDM::CRDArborNode*> m_connector;
      m_connector.first  = m_nodeCol[in]; 
      m_connector.second = m_parentNodes[ip]; 
      m_connectorCol.push_back( m_connector );
    }
  }

  //Clean the connections in tree
//std::cout<<" Connection size: "<<m_connectorCol.size()<<std::endl;

  m_tree.CleanConnection( m_connectorCol );
  m_tree.NodeClassification(); 

  return StatusCode::SUCCESS;
}


StatusCode ArborClusteringAlg::DepartArborTree( CRDEcalEDM::CRDArborTree& m_tree, 
                                                std::vector<CRDEcalEDM::CRDArborTree>& m_departedTrees,
                                                std::vector<CRDEcalEDM::CRDArborNode*>& m_isoNodes ) 
{
//std::cout<<std::endl;
//std::cout<<"ArborClusteringAlg::DepartArborTree  Input tree node size: "<<m_tree.GetNodes().size()<<std::endl;
//std::cout<<"ArborClusteringAlg::DepartArborTree  Input isonode size: "<<m_isoNodes.size()<<std::endl;

  std::vector<CRDEcalEDM::CRDArborNode*> m_rootNodes; m_rootNodes.clear(); 
  for(int in=0; in<m_tree.GetNodes().size(); in++){
    CRDEcalEDM::CRDArborNode* m_node = m_tree.GetNodes()[in];
//printf("  Debug: Loop in node: (%.2f, %.2f, %.2f, %d) \n", m_node->GetPosition().x(), m_node->GetPosition().y(), m_node->GetPosition().z(), m_node->GetType() );
    if(m_node->GetType()==0){ m_isoNodes.push_back(m_node); m_tree.CleanNode(m_node); in--; }
    if(m_node->GetType()==4 || m_node->GetType()==5) m_rootNodes.push_back(m_node); 
//std::cout<<"  Debug: Iso node size: "<<m_isoNodes.size()<<"  Root node size: "<<m_rootNodes.size()<<std::endl;
  }
//std::cout<<"ArborClusteringAlg::DepartArborTree  isonode size after classification: "<<m_isoNodes.size()<<std::endl;
//std::cout<<"ArborClusteringAlg::DepartArborTree  tree node size after classification: "<<m_tree.GetNodes().size()<<std::endl;
//std::cout<<"ArborClusteringAlg::DepartArborTree  Root node size: "<<m_rootNodes.size()<<std::endl;
//std::cout<<std::endl;

  //Only one tree: 
  if( m_rootNodes.size()<=1 ){ 
    if(m_rootNodes.size()==0) std::cout<<"ERROR in ArborClusteringAlg::DepartArborTree: Did not find root node! Check it!"<<std::endl; 
    m_departedTrees.push_back(m_tree);
    return StatusCode::SUCCESS;
  }

  //Multiple trees: 
  else{
  for(int ir=0; ir<m_rootNodes.size(); ir++){
    CRDEcalEDM::CRDArborTree m_newTree; m_newTree.Clear(); 
    m_newTree.CreateWithRoot( m_rootNodes[ir] ); //Create a tree with recursion 
    m_departedTrees.push_back(m_newTree);
  }}

  return StatusCode::SUCCESS;
}


StatusCode ArborClusteringAlg::MergeNeighborTree( std::vector<CRDEcalEDM::CRDArborTree>& m_inTreeCol, 
                                                  std::vector<CRDEcalEDM::CRDArborTree>& m_outTreeCol )
{
  if(m_inTreeCol.size()<=1){ m_outTreeCol=m_inTreeCol; return StatusCode::SUCCESS; }

  std::vector<CRDEcalEDM::CRDArborNode*> m_rootCol; m_rootCol.clear();
  m_rootCol.resize(m_inTreeCol.size());
  for(int it=0; it<m_inTreeCol.size(); it++)
    m_rootCol[it] = m_inTreeCol[it].GetRootNode();

  //Merge trees whose root nodes are in the same layer and close to each other. 
  for(int in=0; in<m_rootCol.size(); in++){
  for(int jn=in+1; jn<m_rootCol.size(); jn++){
    if( m_rootCol[in]->GetDlayer() != m_rootCol[jn]->GetDlayer() ) continue; 
    if( (m_rootCol[in]->GetPosition()-m_rootCol[jn]->GetPosition()).Mag()<settings.Rth_nbrRoot ){
      CRDEcalEDM::CRDArborTree m_newTree;
      m_newTree.AddNode( m_inTreeCol[in].GetNodes() );
      m_newTree.AddNode( m_inTreeCol[jn].GetNodes() );
      m_outTreeCol.push_back( m_newTree );
    }
  }}

  return StatusCode::SUCCESS;
}

/*
StatusCode ArborClusteringAlg::MergeBranches( std::vector<CRDEcalEDM::CRDArborTree>& m_inTreeCol, 
                                              std::vector<CRDEcalEDM::CRDArborTree>& m_goodTreeCol, 
                                              std::vector<CRDEcalEDM::CRDArborTree>& m_badTreeCol)
{
  if(m_inTreeCol.size()==0) return StatusCode::SUCCESS; 

  m_goodTreeCol.clear(); m_badTreeCol.clear(); 

  for(int it=0; it<m_inTreeCol.size(); it++){
    if( (m_inTreeCol[it].GetMaxDlayer()-m_inTreeCol[it].GetMinDlayer()) >= settings.fl_GoodTreeLevel && 
         m_inTreeCol[it].GetNodes().size()>=settings.th_GoodTreeNodes) 
      m_goodTreeCol.push_back( m_inTreeCol[it] ); 
    else m_badTreeCol.push_back( m_inTreeCol[it] );
  }

//std::cout<<"  ArborClusteringAlg::MergeBranches  Ngoodtrees: "<<m_goodTreeCol.size()<<"  Nbadtrees: "<<m_badTreeCol.size()<<std::endl;

//  for(int it=0; it<m_badTreeCol.size(); it++){
//    CRDEcalEDM::CRDArborTree m_tree = GetClosestTree( m_badTreeCol[it], m_goodTreeCol );
//    std::vector<CRDEcalEDM::CRDArborTree>::iterator iter = find(m_goodTreeCol.begin(), m_goodTreeCol.end(), m_tree);
//    if( iter==m_goodTreeCol.end() ) { cout<<"Warning in MergeBranches: Tree Merging Fail!"<<std::endl; continue; }

//    iter->AddNode( m_badTreeCol[it].GetNodes() );
//  }


  return StatusCode::SUCCESS;
}


CRDEcalEDM::CRDArborTree ArborClusteringAlg::GetClosestTree( CRDEcalEDM::CRDArborTree m_badTree, std::vector<CRDEcalEDM::CRDArborTree> m_goodTreeCol ){

  CRDEcalEDM::CRDArborTree m_tree; m_tree.Clear(); 

  TVector3 cent_bad = m_badTree.GetBarycenter();
  double minR = 9999; 
  int index = -1; 
  for(int it=0; it<m_goodTreeCol.size(); it++){
    TVector3 cent_good = m_goodTreeCol[it].GetBarycenter();
    double dis = (cent_bad-cent_good).Mag();
    if(dis<minR) { minR = dis; index = it;}
  }

  if(index<0){ std::cout<<"Warning in ArborClusteringAlg::GetClosestTree: Can not find closest good tree!"<<std::endl; return m_tree; }
  m_tree = m_goodTreeCol[index];
  return m_tree; 
}
*/
#endif
