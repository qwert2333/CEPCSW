#ifndef _CRD_ARBORNODE_
#define _CRD_ARBORNODE_
#include "CRDEcalEDMSvc/CRDCaloHit2DShower.h"
#include "CRDEcalEDMSvc/CRDArborTree.h"
#include "TVector3.h"
#include <vector>

namespace CRDEcalEDM{

  class CRDArborTree; 
  class CRDArborNode{
  
  public:

    CRDArborNode(){};
    CRDArborNode( CRDEcalEDM::CRDCaloHit2DShower _shower ) ;  
    ~CRDArborNode() { Clear(); }

    inline bool operator == (const CRDArborNode &x) const{
      return (pos==x.GetPosition() && En==x.GetEnergy() && Type==x.GetType()) ;
    }
   
    void Clear() {
      pos.SetXYZ(0.,0.,0.); Rref.SetXYZ(0.,0.,0.); En=0; shower.Clear(); Type=0;
      parentNodes.clear(); daughterNodes.clear(); 
    } 

    TVector3 GetPosition() const { return pos; }
    TVector3 GetRefDir(double wf, double wb) const; 
    double GetEnergy() const { return En; }
    int GetDlayer() const { return shower.getDlayer();  }
    int GetType() const { return Type; }
    CRDEcalEDM::CRDCaloHit2DShower GetOriginShower() const { return shower; }
    std::vector<CRDEcalEDM::CRDArborNode*> GetParentNodes() const { return parentNodes; }
    std::vector<CRDEcalEDM::CRDArborNode*> GetDaughterNodes() const { return daughterNodes; } 

    void SetPosition(TVector3& _vecP) { pos=_vecP; }
    void SetRefDir(TVector3& _vecRef) { Rref=_vecRef; }
    void SetOriginBarShower( CRDEcalEDM::CRDCaloHit2DShower& _sh ) { shower=_sh; }
    void SetEnergy( double _en ) { En=_en; }
    void SetType( int _type ) { Type=_type; }

    void ConnectDaughter( CRDEcalEDM::CRDArborNode* _Dnode);
    void ConnectParent  ( CRDEcalEDM::CRDArborNode* _Pnode);
    void DisconnectDaughter( CRDEcalEDM::CRDArborNode* _Dnode);
    void DisconnectParent  ( CRDEcalEDM::CRDArborNode* _Pnode);

    bool isInTree( CRDEcalEDM::CRDArborTree& _tree ) const;

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
    CRDEcalEDM::CRDCaloHit2DShower shower;  

    std::vector<CRDEcalEDM::CRDArborNode*> parentNodes; 
    std::vector<CRDEcalEDM::CRDArborNode*> daughterNodes; 


  };

};
#endif
