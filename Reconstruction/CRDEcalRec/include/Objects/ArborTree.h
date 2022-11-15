#ifndef _ARBORTREE_H
#define _ARBORTREE_H
#include "Objects/ArborNode.h"
#include "Objects/Calo3DCluster.h"
#include "TVector3.h"

namespace PandoraPlus{
  class ArborNode; 

  class ArborTree{
  
  public:

    ArborTree(){};
    ArborTree(std::vector<PandoraPlus::ArborNode*> _nodeCol) : Nodes(_nodeCol) {};

    inline bool operator == (const ArborTree &x) const{
      return Nodes == x.getNodes(); 
    }
    void Clear() { 
      for(int i=0; i<Nodes.size(); i++) delete Nodes[i];
      Nodes.clear(); 
      pBarycent.SetXYZ(0.,0.,0.); 
    } 

    std::vector<PandoraPlus::ArborNode*> getNodes() const { return Nodes; }
    std::vector<PandoraPlus::ArborNode*> getNodes(int Layer) const; 
    PandoraPlus::ArborNode* getRootNode() const; 
    int getMinDlayer() const; 
    int getMaxDlayer() const;
    TVector3 getBarycenter(); 

    void addNode( PandoraPlus::ArborNode* _node );
    void addNode( std::vector<PandoraPlus::ArborNode*> _nodes ); 
    void createWithRoot( PandoraPlus::ArborNode* _root ); 
    void createWithTopo( PandoraPlus::ArborNode* _node );
    void createWithNode( PandoraPlus::ArborNode* _node, bool f_t2b ); 
    void cleanNode( PandoraPlus::ArborNode* _node ) ;
    void sortNodes() { std::sort(Nodes.begin(), Nodes.end(), compLayer); } 
    void cleanConnection( std::vector< std::pair<PandoraPlus::ArborNode*, PandoraPlus::ArborNode*> >& _connectionCol ) ;
    void nodeClassification(); 

    void printTree() const; 
    Calo3DCluster* convertTreeToCluster() const; 

  private: 
    std::vector<PandoraPlus::ArborNode*> Nodes; 
    TVector3 pBarycent; 

    static bool compLayer( PandoraPlus::ArborNode* node1, PandoraPlus::ArborNode* node2 ); 
    //static bool ArborTree::compLayer( PandoraPlus::ArborNode& node1, PandoraPlus::ArborNode& node2 ) { return node1.GetDlayer()<node2.GetDlayer(); }

  };

};
#endif
