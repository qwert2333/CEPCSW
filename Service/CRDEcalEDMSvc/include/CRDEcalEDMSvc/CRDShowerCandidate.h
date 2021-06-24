#ifndef _CRD_EXPEMSHOWER_
#define _CRD_EXPEMSHOWER_

namespace CRDEcalEDM {

  //Expected shower in one layer. 
  class CRDShowerCandidate{
  public: 

    CRDShowerCandidate(){}; 
    void Clear(){  Dlayer=0; ExpEshower=0; ExpEseed=0; ExpDepth=0; ExpPos.SetXYZ(0,0,0);  }

    inline bool operator == (const CRDShowerCandidate &x) const{
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

  };

};
#endif
