#ifndef _CRD_ARBORTREE_
#define _CRD_ARBORTREE_
#include "CRDEcalEDMSvc/CRDArborNode.h"
#include "CRDEcalEDMSvc/CRDArborTree.h"
#include "TVector3.h"

namespace CRDEcalEDM{
  class CRDArborNode; 

  class CRDArborTreeData{
  
  public:

    CRDArborTreeData(){};
    CRDArborTreeData(std::vector<CRDEcalEDM::CRDArborNode> _nodeCol) : Nodes(_nodeCol) {};

    inline bool operator == (const CRDArborTreeData &x) const{
      return Nodes == x.GetNodes(); 
    }
    void Clear() { Nodes.clear(); pBarycent.SetXYZ(0.,0.,0.); } 

    std::vector<CRDEcalEDM::CRDArborNode> GetNodes() const { return Nodes; }
    CRDEcalEDM::CRDArborNode GetRootNode() const; 
    int GetMinDlayer() const; 
    int GetMaxDlayer() const;
    TVector3 GetBarycenter(); 

    void AddNode( CRDEcalEDM::CRDArborNode _node );
    void AddNode( std::vector<CRDEcalEDM::CRDArborNode> _nodes ); 
    void SortNodes() { std::sort(Nodes.begin(), Nodes.end(), compLayer); } 
    void NodeClassification(); 

    void PrintTree() const; 

  private: 
    std::vector<CRDEcalEDM::CRDArborNode> Nodes; 
    TVector3 pBarycent; 

    static bool compLayer( CRDEcalEDM::CRDArborNode node1, CRDEcalEDM::CRDArborNode node2 ); 

  };

};
#endif
