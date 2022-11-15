#ifndef _ARBORNODE_C
#define _ARBORNODE_C

#include "Objects/ArborNode.h"

namespace PandoraPlus{

  ArborNode::ArborNode( PandoraPlus::CaloHit* _calohit ){
    calohit = _calohit; 
    En = _calohit->getEnergy();
    pos.SetXYZ( _calohit->getPosition().x(), _calohit->getPosition().y(), _calohit->getPosition().z() );
    Type = 0;
    parentNodes.clear(); 
    daughterNodes.clear(); 
  };

  TVector3 ArborNode::getRefDir(double wf, double wb) const{
    TVector3 sumVf(0.,0.,0.);
    TVector3 sumVb(0.,0.,0.);
    for(int i=0; i<parentNodes.size(); i++)
      sumVf += wf*( pos - parentNodes[i]->getPosition() );
    
    for(int i=0; i<daughterNodes.size(); i++)
      sumVb += wb*( daughterNodes[i]->getPosition() - pos );
    
    return sumVf + sumVb; 
  };

  
  void ArborNode::ConnectDaughter( PandoraPlus::ArborNode* _Dnode ){
    daughterNodes.push_back(_Dnode);

    std::vector<PandoraPlus::ArborNode*> _Dnode_parents = _Dnode->getParentNodes(); 
    std::vector<PandoraPlus::ArborNode*>::iterator iter = find(_Dnode_parents.begin(), _Dnode_parents.end(), this);
    if( iter == _Dnode_parents.end() )  _Dnode->ConnectParent(this);
  };

  
  void ArborNode::ConnectParent( PandoraPlus::ArborNode* _Pnode ){
    parentNodes.push_back(_Pnode);

    std::vector<PandoraPlus::ArborNode*> _Pnode_daughter = _Pnode->getDaughterNodes();
    std::vector<PandoraPlus::ArborNode*>::iterator iter = find(_Pnode_daughter.begin(), _Pnode_daughter.end(), this);
    if( iter == _Pnode_daughter.end() )  _Pnode->ConnectDaughter(this);
  };


  void ArborNode::DisconnectDaughter( PandoraPlus::ArborNode* _Dnode ){
    std::vector<PandoraPlus::ArborNode*>::iterator iter = find( daughterNodes.begin(), daughterNodes.end(), _Dnode );
    if( iter!=daughterNodes.end() )  daughterNodes.erase( iter );
    
    std::vector<PandoraPlus::ArborNode*> m_nodes = _Dnode->getParentNodes();
    iter = find(m_nodes.begin(), m_nodes.end(), this);
    if( iter!=m_nodes.end() ) _Dnode->DisconnectParent(this);
  };


  void ArborNode::DisconnectParent( PandoraPlus::ArborNode* _Pnode ){
    std::vector<PandoraPlus::ArborNode*>::iterator iter = find( parentNodes.begin(), parentNodes.end(), _Pnode );
    if( iter!=parentNodes.end() ) parentNodes.erase( iter );

    std::vector<PandoraPlus::ArborNode*> m_nodes = _Pnode->getDaughterNodes();
    iter = find(m_nodes.begin(), m_nodes.end(), this);
    if( iter!=m_nodes.end() ) _Pnode->DisconnectDaughter(this);
  };

  bool ArborNode::isInTree( const PandoraPlus::ArborTree* _tree ) const{
    std::vector<PandoraPlus::ArborNode*> nodes = _tree->getNodes();
    std::vector<PandoraPlus::ArborNode*>::iterator iter = find( nodes.begin(), nodes.end(), this );
    if( iter == nodes.end() ) return false;
    else return true; 
  };


};
#endif
