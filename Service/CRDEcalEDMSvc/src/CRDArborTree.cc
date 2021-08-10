#ifndef _CRD_ARBORTREE_C
#define _CRD_ARBORTREE_C
#include "CRDEcalEDMSvc/CRDArborTree.h"

namespace CRDEcalEDM{

  CRDEcalEDM::CRDArborNode* CRDArborTree::GetRootNode() const{
    std::vector<CRDEcalEDM::CRDArborNode*> _rootCol; _rootCol.clear(); 
    for(int i=0; i<Nodes.size(); i++)
      if( Nodes[i]!=NULL && (Nodes[i]->GetType()==4 || Nodes[i]->GetType()==5) ) _rootCol.push_back(Nodes[i]);

    if(_rootCol.size()==0){ std::cout<<"Error in GetRootNode: No Root node."<<std::endl; return nullptr; }
    if(_rootCol.size()>=2){ std::cout<<"Error in GetRootNode: >=2 root node in tree, need splitting."<<std::endl; return nullptr; }
    if(_rootCol.size()==1) return _rootCol[0]; 
  };


  std::vector<CRDEcalEDM::CRDArborNode*> CRDArborTree::GetNodes(int Layer) const{
    std::vector<CRDEcalEDM::CRDArborNode*> m_nodes; m_nodes.clear(); 
    for(int in=0; in<Nodes.size(); in++) 
      if(Nodes[in]->GetDlayer() == Layer) m_nodes.push_back(Nodes[in]);
    return m_nodes; 
  }


  int CRDArborTree::GetMinDlayer() const{
    int min=999; 
    for(int i=0; i<Nodes.size(); i++)
      if( Nodes[i]->GetDlayer()<min && Nodes[i]->GetDlayer()>0 ) min = Nodes[i]->GetDlayer(); 
    if(min==999){ std::cout<<"Error in GetMinDlayer: Did not get minimum Dlayer. "<<std::endl; return -1; }
    return min; 
  };


  int CRDArborTree::GetMaxDlayer() const{
    int max=-999;
    for(int i=0; i<Nodes.size(); i++)
      if( Nodes[i]->GetDlayer()>max && Nodes[i]->GetDlayer()>0 ) max = Nodes[i]->GetDlayer();
    if(max<=0){ std::cout<<"Error in GetMaxDlayer: Did not get maximum Dlayer. "<<std::endl; return -1;}
    return max;
  };


  TVector3 CRDArborTree::GetBarycenter() {
    TVector3 vec(0., 0., 0.);
    double totE = 0;
    for(int i=0; i<Nodes.size(); i++) totE += Nodes[i]->GetEnergy(); 
    for(int i=0; i<Nodes.size(); i++) vec += (Nodes[i]->GetEnergy()/totE)*Nodes[i]->GetPosition(); 
    pBarycent = vec; 

    return pBarycent; 
  };
  
  void CRDArborTree::AddNode( CRDEcalEDM::CRDArborNode* _node ) {
    std::vector<CRDEcalEDM::CRDArborNode*>::iterator iter = find(Nodes.begin(), Nodes.end(), _node);
    if( iter==Nodes.end() ) Nodes.push_back( _node );
  }

  void CRDArborTree::AddNode( std::vector<CRDEcalEDM::CRDArborNode*> _nodes ) {
    for(int i=0; i<_nodes.size(); i++) AddNode( _nodes[i] );
  };


  void CRDArborTree::CreateWithRoot(CRDEcalEDM::CRDArborNode* _root) {
    Nodes.clear(); 
    Nodes.push_back( _root );
    for(int i=0; i< _root->GetDaughterNodes().size(); i++) this->CreateWithNode( _root->GetDaughterNodes()[i], true );
  };

  void CRDArborTree::CreateWithTopo( CRDEcalEDM::CRDArborNode* _node ){
    Nodes.clear();
    this->CreateWithNode( _node, false);
  }

  void CRDArborTree::CreateWithNode( CRDEcalEDM::CRDArborNode* _node, bool f_t2b ){
    std::vector<CRDEcalEDM::CRDArborNode*>::iterator iter = find(Nodes.begin(), Nodes.end(), _node);
    if( iter!=Nodes.end() ) return; 

    this->Nodes.push_back(_node);
    if(f_t2b) 
      for(int i=0; i< _node->GetDaughterNodes().size(); i++) this->CreateWithNode( _node->GetDaughterNodes()[i], true );
    else{
      for(int i=0; i< _node->GetDaughterNodes().size(); i++) this->CreateWithNode( _node->GetDaughterNodes()[i], false );
      for(int i=0; i< _node->GetParentNodes().size(); i++) this->CreateWithNode( _node->GetParentNodes()[i], false );
    }
  }


  void CRDArborTree::CleanNode( CRDEcalEDM::CRDArborNode* _node ){
    std::vector<CRDEcalEDM::CRDArborNode*>::iterator iter = find( Nodes.begin(), Nodes.end(), _node );
    if(iter!=Nodes.end()) Nodes.erase(iter);
  };


  void CRDArborTree::CleanConnection( std::vector< std::pair<CRDEcalEDM::CRDArborNode*, CRDEcalEDM::CRDArborNode*> >& _connectionCol){
    for(int i=0; i<_connectionCol.size(); i++){
//std::cout<<"CleanConnection: #"<<i<<", daughter and parent node: "<<std::endl;
//printf("(%.2f, %.2f, %.2f, %d) \t", _connectionCol[i].first->GetPosition().x(), _connectionCol[i].first->GetPosition().y(), _connectionCol[i].first->GetPosition().z(), _connectionCol[i].first->GetType());
//printf("(%.2f, %.2f, %.2f, %d) \n", _connectionCol[i].second->GetPosition().x(), _connectionCol[i].second->GetPosition().y(), _connectionCol[i].second->GetPosition().z(), _connectionCol[i].second->GetType());

//      _connectionCol[i].first->DisconnectDaughter( _connectionCol[i].second );
      _connectionCol[i].first->DisconnectParent  ( _connectionCol[i].second );
    }
  };


  void CRDArborTree::NodeClassification() {
    for(int i=0; i<Nodes.size(); i++){
      int Nparent = Nodes[i]->GetParentNodes().size(); 
      int Ndaughter = Nodes[i]->GetDaughterNodes().size(); 
      if(Ndaughter==0 && Nparent==1)  Nodes[i]->SetType(1);     // Leaf node
      else if(Ndaughter==1 && Nparent==1) Nodes[i]->SetType(2); // Joint node
      else if(Ndaughter>1  && Nparent==1) Nodes[i]->SetType(3); // Branch node
      else if(Ndaughter==1 && Nparent==0) Nodes[i]->SetType(4); // Root node
      else if(Ndaughter>1  && Nparent==0) Nodes[i]->SetType(5); // Star root node
      else if(Ndaughter==0 && Nparent==0) Nodes[i]->SetType(0); // Isolated node
      else std::cout<<"Warning: Unknown node type appear!  Ndaughter="<<Ndaughter<<"  Nparent="<<Nparent<<" Dlayer="<<Nodes[i]->GetDlayer()<<std::endl;
    }
  };

  void CRDArborTree::PrintTree() const{
    std::cout<<"N node: "<<Nodes.size()<<std::endl;
    std::cout<<"Nodes in this tree: (x, y, z, Type)"<<std::endl;
    for(int i=0; i<Nodes.size(); i++)
      printf("(%.2f, %.2f, %.2f, %d) \n", Nodes[i]->GetPosition().x(), Nodes[i]->GetPosition().y(), Nodes[i]->GetPosition().z(), Nodes[i]->GetType());
    std::cout<<std::endl;
  };

  CRDCaloHit3DShower CRDArborTree::ConvertTreeToCluster() const{
    CRDEcalEDM::CRDCaloHit3DShower m_clus; m_clus.Clear();
    for(int in=0; in<Nodes.size(); in++){
      CRDEcalEDM::CRDCaloHit2DShower m_shower = Nodes[in]->GetOriginShower();
      m_clus.AddShower( m_shower );
    }
    return m_clus; 
  };

  bool CRDArborTree::compLayer( CRDEcalEDM::CRDArborNode* node1, CRDEcalEDM::CRDArborNode* node2 ) { return node1->GetDlayer()>node2->GetDlayer(); }

};
#endif
