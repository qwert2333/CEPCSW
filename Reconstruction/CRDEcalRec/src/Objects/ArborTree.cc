#ifndef _ARBORTREE_C
#define _ARBORTREE_C
#include "Objects/ArborTree.h"

namespace PandoraPlus{

  ArborNode* ArborTree::getRootNode() const{
    std::vector<ArborNode*> _rootCol; _rootCol.clear(); 
    for(int i=0; i<Nodes.size(); i++)
      if( Nodes[i]!=NULL && (Nodes[i]->getType()==4 || Nodes[i]->getType()==5) ) _rootCol.push_back(Nodes[i]);

    if(_rootCol.size()==0){ std::cout<<"Error in getRootNode: No Root node."<<std::endl; return nullptr; }
    if(_rootCol.size()>=2){ std::cout<<"Error in getRootNode: >=2 root node in tree, need splitting."<<std::endl; return nullptr; }
    if(_rootCol.size()==1) return _rootCol[0]; 
  };


  std::vector<PandoraPlus::ArborNode*> ArborTree::getNodes(int Layer) const{
    std::vector<PandoraPlus::ArborNode*> m_nodes; m_nodes.clear(); 
    for(int in=0; in<Nodes.size(); in++) 
      if(Nodes[in]->getLayer() == Layer) m_nodes.push_back(Nodes[in]);
    return m_nodes; 
  }


  int ArborTree::getMinDlayer() const{
    int min=999; 
    for(int i=0; i<Nodes.size(); i++)
      if( Nodes[i]->getLayer()<min && Nodes[i]->getLayer()>0 ) min = Nodes[i]->getLayer(); 
    if(min==999){ std::cout<<"Error in getMinDlayer: Did not get minimum Dlayer. "<<std::endl; return -1; }
    return min; 
  };


  int ArborTree::getMaxDlayer() const{
    int max=-999;
    for(int i=0; i<Nodes.size(); i++)
      if( Nodes[i]->getLayer()>max && Nodes[i]->getLayer()>0 ) max = Nodes[i]->getLayer();
    if(max<=0){ std::cout<<"Error in getMaxDlayer: Did not get maximum Dlayer. "<<std::endl; return -1;}
    return max;
  };


  TVector3 ArborTree::getBarycenter() {
    TVector3 vec(0., 0., 0.);
    double totE = 0;
    for(int i=0; i<Nodes.size(); i++) totE += Nodes[i]->getEnergy(); 
    for(int i=0; i<Nodes.size(); i++) vec += (Nodes[i]->getEnergy()/totE)*Nodes[i]->getPosition(); 
    pBarycent = vec; 

    return pBarycent; 
  };
  
  void ArborTree::addNode( PandoraPlus::ArborNode* _node ) {
    std::vector<PandoraPlus::ArborNode*>::iterator iter = find(Nodes.begin(), Nodes.end(), _node);
    if( iter==Nodes.end() ) Nodes.push_back( _node );
  }

  void ArborTree::addNode( std::vector<PandoraPlus::ArborNode*> _nodes ) {
    for(int i=0; i<_nodes.size(); i++) addNode( _nodes[i] );
  };


  void ArborTree::createWithRoot(PandoraPlus::ArborNode* _root) {
    Nodes.clear(); 
    Nodes.push_back( _root );
    for(int i=0; i< _root->getDaughterNodes().size(); i++) this->createWithNode( _root->getDaughterNodes()[i], true );
  };

  void ArborTree::createWithTopo( PandoraPlus::ArborNode* _node ){
    Nodes.clear();
    this->createWithNode( _node, false);
  }

  void ArborTree::createWithNode( PandoraPlus::ArborNode* _node, bool f_t2b ){
    std::vector<PandoraPlus::ArborNode*>::iterator iter = find(Nodes.begin(), Nodes.end(), _node);
    if( iter!=Nodes.end() ) return; 

    this->Nodes.push_back(_node);
    if(f_t2b) 
      for(int i=0; i< _node->getDaughterNodes().size(); i++) this->createWithNode( _node->getDaughterNodes()[i], true );
    else{
      for(int i=0; i< _node->getDaughterNodes().size(); i++) this->createWithNode( _node->getDaughterNodes()[i], false );
      for(int i=0; i< _node->getParentNodes().size(); i++) this->createWithNode( _node->getParentNodes()[i], false );
    }
  }


  void ArborTree::cleanNode( PandoraPlus::ArborNode* _node ){
    std::vector<PandoraPlus::ArborNode*>::iterator iter = find( Nodes.begin(), Nodes.end(), _node );
    if(iter!=Nodes.end()) Nodes.erase(iter);
  };


  void ArborTree::cleanConnection( std::vector< std::pair<PandoraPlus::ArborNode*, PandoraPlus::ArborNode*> >& _connectionCol){
    for(int i=0; i<_connectionCol.size(); i++){
//std::cout<<"CleanConnection: #"<<i<<", daughter and parent node: "<<std::endl;
//printf("(%.2f, %.2f, %.2f, %d) \t", _connectionCol[i].first->getPosition().x(), _connectionCol[i].first->getPosition().y(), _connectionCol[i].first->getPosition().z(), _connectionCol[i].first->getType());
//printf("(%.2f, %.2f, %.2f, %d) \n", _connectionCol[i].second->getPosition().x(), _connectionCol[i].second->getPosition().y(), _connectionCol[i].second->getPosition().z(), _connectionCol[i].second->getType());

//      _connectionCol[i].first->DisconnectDaughter( _connectionCol[i].second );
      _connectionCol[i].first->DisconnectParent  ( _connectionCol[i].second );
    }
  };


  void ArborTree::nodeClassification() {
    for(int i=0; i<Nodes.size(); i++){
      int Nparent = Nodes[i]->getParentNodes().size(); 
      int Ndaughter = Nodes[i]->getDaughterNodes().size(); 
      if(Ndaughter==0 && Nparent==1)  Nodes[i]->setType(1);     // Leaf node
      else if(Ndaughter==1 && Nparent==1) Nodes[i]->setType(2); // Joint node
      else if(Ndaughter>1  && Nparent==1) Nodes[i]->setType(3); // Branch node
      else if(Ndaughter==1 && Nparent==0) Nodes[i]->setType(4); // Root node
      else if(Ndaughter>1  && Nparent==0) Nodes[i]->setType(5); // Star root node
      else if(Ndaughter==0 && Nparent==0) Nodes[i]->setType(0); // Isolated node
      else std::cout<<"Warning: Unknown node type appear!  Ndaughter="<<Ndaughter<<"  Nparent="<<Nparent<<" Dlayer="<<Nodes[i]->getLayer()<<std::endl;
    }
  };

  void ArborTree::printTree() const{
    std::cout<<"N node: "<<Nodes.size()<<std::endl;
    std::cout<<"Nodes in this tree: (x, y, z, Type)"<<std::endl;
    for(int i=0; i<Nodes.size(); i++)
      printf("(%.2f, %.2f, %.2f, %d) \n", Nodes[i]->getPosition().x(), Nodes[i]->getPosition().y(), Nodes[i]->getPosition().z(), Nodes[i]->getType());
    std::cout<<std::endl;
  };

  Calo3DCluster* ArborTree::convertTreeToCluster() const{
    PandoraPlus::Calo3DCluster* m_clus = new Calo3DCluster();
    for(int in=0; in<Nodes.size(); in++){
      CaloHit* m_hit = Nodes[in]->getOriginCaloHit();
      m_clus->addHit( m_hit );
    }
    m_clus->FitAxisHit();
    return m_clus; 
  };

  bool ArborTree::compLayer( PandoraPlus::ArborNode* node1, PandoraPlus::ArborNode* node2 ) { return node1->getLayer()>node2->getLayer(); }

};
#endif
