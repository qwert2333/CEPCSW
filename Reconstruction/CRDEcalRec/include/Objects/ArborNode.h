#ifndef _ARBORNODE_
#define _ARBORNODE_
#include "Objects/CaloHit.h"
#include "TVector3.h"
#include <vector>

namespace PandoraPlus{


  class ArborNode{
  
  public:

    ArborNode(const PandoraPlus::CaloHit* _calohit ) ;  
    ~ArborNode() { Clear(); }

    inline bool operator == (const ArborNode &x) const{
      return (pos==x.getPosition() && En==x.getEnergy() && Type==x.getType() && calohit==x.calohit) ;
    }
   
    void Clear() {
      pos.SetXYZ(0.,0.,0.); Rref.SetXYZ(0.,0.,0.); En=0; calohit=nullptr; Type=0;
      parentNodes.clear(); daughterNodes.clear(); 
    } 

    TVector3 getPosition() const { return pos; }
    TVector3 getRefDir(double wf, double wb) const; 
    double getEnergy() const { return En; }
    int getLayer() const { return calohit->getLayer();  }
    int getType() const { return Type; }
    const PandoraPlus::CaloHit* getOriginCaloHit() const { return calohit; }
    std::vector<PandoraPlus::ArborNode*> getParentNodes() const { return parentNodes; }
    std::vector<PandoraPlus::ArborNode*> getDaughterNodes() const { return daughterNodes; } 

    void setPosition(TVector3& _vecP) { pos=_vecP; }
    void setRefDir(TVector3& _vecRef) { Rref=_vecRef; }
    void setEnergy( double _en ) { En=_en; }
    void setType( int _type ) { Type=_type; }

    void ConnectDaughter( PandoraPlus::ArborNode* _Dnode);
    void ConnectParent  ( PandoraPlus::ArborNode* _Pnode);
    void DisconnectDaughter( PandoraPlus::ArborNode* _Dnode);
    void DisconnectParent  ( PandoraPlus::ArborNode* _Pnode);

  private: 
    TVector3 pos; 
    TVector3 Rref;
    double En;
    int Type; //0: isolate node. 
              //1: leaf node. 
              //2: joint node. 
              //3: branch node.
              //4: root node. 
              //5: star root node. 
    const PandoraPlus::CaloHit* calohit;  

    std::vector<PandoraPlus::ArborNode*> parentNodes; 
    std::vector<PandoraPlus::ArborNode*> daughterNodes; 

    static bool compR( const PandoraPlus::ArborNode* node1, const PandoraPlus::ArborNode* node2 )
      { return node1->getPosition().Mag() < node2->getPosition().Mag(); }

  };

  typedef std::pair<ArborNode*, ArborNode*> ArborLink;
  typedef std::vector<ArborLink> ArborTree; 

};
#endif
