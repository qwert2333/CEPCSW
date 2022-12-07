#ifndef _ARBORNODE_H
#define _ARBORNODE_H
#include "TVector3.h"
#include <vector>


  class ArborNode{
  
  public:

    ArborNode(const edm4hep::CalorimeterHit* _calohit ) ;  
    ~ArborNode() { Clear(); }

    inline bool operator == (const ArborNode &x) const{
      return (pos==x.getPosition() && En==x.getEnergy() && Type==x.getType() && calohit==x.calohit) ;
    }
    inline bool operator < (const ArborNode &x) const{
      return (pos.Mag()<x.getPosition().Mag()) ;
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
    const edm4hep::CalorimeterHit* getOriginCalorimeterHit() const { return calohit; }
    std::vector<ArborNode*> getParentNodes() const { return parentNodes; }
    std::vector<ArborNode*> getDaughterNodes() const { return daughterNodes; } 

    void setPosition(TVector3& _vecP) { pos=_vecP; }
    void setRefDir(TVector3& _vecRef) { Rref=_vecRef; }
    void setEnergy( double _en ) { En=_en; }
    void setType( int _type ) { Type=_type; }

    void ConnectDaughter( ArborNode* _Dnode);
    void ConnectParent  ( ArborNode* _Pnode);
    void DisconnectDaughter( ArborNode* _Dnode);
    void DisconnectParent  ( ArborNode* _Pnode);

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
    const edm4hep::CalorimeterHit* calohit;  

    std::vector<ArborNode*> parentNodes; 
    std::vector<ArborNode*> daughterNodes; 

  };

  typedef std::pair<ArborNode*, ArborNode*> ArborLink;
  typedef std::vector<ArborLink> ArborTree; 

};
#endif
