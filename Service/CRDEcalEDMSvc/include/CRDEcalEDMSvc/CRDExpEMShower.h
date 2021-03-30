#ifndef _CRD_EXPEMSHOWER_
#define _CRD_EXPEMSHOWER_

namespace CRDEcalEDM {

  //Expected shower in one layer. 
  class CRDExpEMShower{
  public: 

    CRDExpEMShower(){}; 
    void Clear(){  Dlayer=0; ExpEshower=0; ExpEseed=0; ExpDepth=0; ExpPos.SetXYZ(0,0,0);  }

    int Dlayer; 
    double ExpEshower; 
    double ExpEseed; 
    double ExpDepth;
    TVector3 ExpPos;

  };

};
#endif
