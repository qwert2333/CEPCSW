#ifndef _CRD_SHADOWCLUS_H
#define _CRD_SHADOWCLUS_H

namespace CRDEcalEDM {

  //Expected shower in one layer. 
  class CRDShadowCluster{
  public: 

    CRDShadowCluster(){}; 
    void Clear(){  Dlayer=0; ExpEshower=0; ExpEseed=0; ExpDepth=0; ExpPos.SetXYZ(0,0,0); slayer=-1; }

    inline bool operator == (const CRDShadowCluster &x) const{
      return ( Dlayer==x.Dlayer  &&
               ExpEshower == x.ExpEshower  &&
               ExpEseed == x.ExpEseed  &&
               ExpDepth == x.ExpDepth &&
               ExpPos == x.ExpPos &&
               Type == x.Type
             );
    }

    int Dlayer; 
    double ExpEshower; 
    double ExpEseed; 
    double ExpDepth;
    TVector3 ExpPos;
    int Type;  //1 for track-type, 0 for neutral-type
    int slayer; 

  };

};
#endif
