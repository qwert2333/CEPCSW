#ifndef _CRD_CALOHOUGHOBJECT_H
#define _CRD_CALOHOUGHOBJECT_H
#include "Objects/CRDCaloBarShower.h"
#include "TVector2.h"
#include "TF1.h"
namespace CRDEcalEDM {

  class CRDHoughObject{

  public: 
    CRDHoughObject() {};
    ~CRDHoughObject() {};

    inline bool operator == (const CRDHoughObject &x) const{
      return  originLocalMax == x.originLocalMax;
    }

		CRDEcalEDM::CRDCaloBarShower* getLocalMax() const { return originLocalMax; }
		TVector2 getConformPointU() const { return ConformalPoint_u; }
		TVector2 getConformPointD() const { return ConformalPoint_d; }
    TF1 getHoughLineU() const { return HoughLine_u; }
    TF1 getHoughLineD() const { return HoughLine_d; }

    void Clear() { originLocalMax=nullptr; ConformalPoint_u.SetX(0.); ConformalPoint_u.SetY(0.); ConformalPoint_d.SetX(0.); ConformalPoint_d.SetY(0.); }
    void SetLocalMax( CRDEcalEDM::CRDCaloBarShower* _localmax ) { originLocalMax=_localmax; }
    void SetConformalPoint(TVector2& _uvec, TVector2& _dvec) { ConformalPoint_u=_uvec; ConformalPoint_d=_dvec; }
		void SetHoughLine(TF1& _func_u, TF1& _func_d) { HoughLine_u=_func_u; HoughLine_d=_func_d; }
		void SetSlayer(int _slayer) { Slayer=_slayer; }

  private: 
	  int Slayer;
    CRDEcalEDM::CRDCaloBarShower* originLocalMax;  //Local max
		TVector2 ConformalPoint_u; 
		TVector2 ConformalPoint_d; 
		TF1 HoughLine_u;
		TF1 HoughLine_d;

};

};
#endif
