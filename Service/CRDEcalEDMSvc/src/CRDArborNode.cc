#ifndef _CRD_ARBORNODE_C
#define _CRD_ARBORNODE_C

#include "CRDEcalEDMSvc/CRDArborNode.h"

namespace CRDEcalEDM{

  CRDArborNode::CRDArborNode( CRDEcalEDM::CRDCaloHit2DShower _shower ){
    shower = _shower; 
    En = _shower.getShowerE();
    pos.SetXYZ( _shower.getPos().x(), _shower.getPos().y(), _shower.getPos().z() );
    Type = 0;
    parentNodes.clear(); 
    daughterNodes.clear(); 
  };

  TVector3 CRDArborNode::GetRefDir(double wf, double wb) const{
    TVector3 sumVf(0.,0.,0.);
    TVector3 sumVb(0.,0.,0.);
    for(int i=0; i<parentNodes.size(); i++)
      sumVf += wf*( pos - parentNodes[i]->GetPosition() );
    
    for(int i=0; i<daughterNodes.size(); i++)
      sumVb += wb*( daughterNodes[i]->GetPosition() - pos );
    
    return sumVf + sumVb; 
  };

  
  void CRDArborNode::ConnectDaughter( CRDEcalEDM::CRDArborNode* _Dnode ){
    daughterNodes.push_back(_Dnode);

    std::vector<CRDEcalEDM::CRDArborNode*> _Dnode_parents = _Dnode->GetParentNodes(); 
    std::vector<CRDEcalEDM::CRDArborNode*>::iterator iter = find(_Dnode_parents.begin(), _Dnode_parents.end(), this);
    if( iter == _Dnode_parents.end() )  _Dnode->ConnectParent(this);
  };

  
  void CRDArborNode::ConnectParent( CRDEcalEDM::CRDArborNode* _Pnode ){
    parentNodes.push_back(_Pnode);

    std::vector<CRDEcalEDM::CRDArborNode*> _Pnode_daughter = _Pnode->GetDaughterNodes();
    std::vector<CRDEcalEDM::CRDArborNode*>::iterator iter = find(_Pnode_daughter.begin(), _Pnode_daughter.end(), this);
    if( iter == _Pnode_daughter.end() )  _Pnode->ConnectDaughter(this);
  };


  void CRDArborNode::DisconnectDaughter( CRDEcalEDM::CRDArborNode* _Dnode ){
    std::vector<CRDEcalEDM::CRDArborNode*>::iterator iter = find( daughterNodes.begin(), daughterNodes.end(), _Dnode );
    if( iter!=daughterNodes.end() )  daughterNodes.erase( iter );
    
    std::vector<CRDEcalEDM::CRDArborNode*> m_nodes = _Dnode->GetParentNodes();
    iter = find(m_nodes.begin(), m_nodes.end(), this);
    if( iter!=m_nodes.end() ) _Dnode->DisconnectParent(this);
  };


  void CRDArborNode::DisconnectParent( CRDEcalEDM::CRDArborNode* _Pnode ){
    std::vector<CRDEcalEDM::CRDArborNode*>::iterator iter = find( parentNodes.begin(), parentNodes.end(), _Pnode );
    if( iter!=parentNodes.end() ) parentNodes.erase( iter );

    std::vector<CRDEcalEDM::CRDArborNode*> m_nodes = _Pnode->GetDaughterNodes();
    iter = find(m_nodes.begin(), m_nodes.end(), this);
    if( iter!=m_nodes.end() ) _Pnode->DisconnectDaughter(this);
  };

  bool CRDArborNode::isInTree( CRDEcalEDM::CRDArborTree& _tree ) const{
    std::vector<CRDEcalEDM::CRDArborNode*> nodes = _tree.GetNodes();
    std::vector<CRDEcalEDM::CRDArborNode*>::iterator iter = find( nodes.begin(), nodes.end(), this );
    if( iter == nodes.end() ) return false;
    else return true; 
  };


};
#endif
