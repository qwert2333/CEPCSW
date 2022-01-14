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

    //inline bool operator == (const CRDHoughObject &x) const{
    //  return  originLocalMax == x.originLocalMax;
    //}

		CRDEcalEDM::CRDCaloBarShower getLocalMax() const { return originLocalMax; }
		TVector2 getConformPointUR() const { return ConformalPoint+TVector2(cellSize/2., cellSize/2.);   }
		TVector2 getConformPointUL() const { return ConformalPoint+TVector2(-cellSize/2., cellSize/2.);  }
		TVector2 getConformPointDR() const { return ConformalPoint+TVector2(cellSize/2., -cellSize/2.);  }
		TVector2 getConformPointDL() const { return ConformalPoint+TVector2(-cellSize/2., -cellSize/2.); }
    TF1 getHoughLineUR() const { return HoughLine_ur; }
    TF1 getHoughLineUL() const { return HoughLine_ul; }
    TF1 getHoughLineDR() const { return HoughLine_dr; }
    TF1 getHoughLineDL() const { return HoughLine_dl; }

    void Clear() { originLocalMax.Clear(); ConformalPoint.SetX(0.); ConformalPoint.SetY(0.);}
    void SetCellSize(double _cell) { cellSize=_cell; }
    void SetLocalMax( CRDEcalEDM::CRDCaloBarShower _localmax ) { originLocalMax=_localmax; }
    void SetConformalPoint(TVector2& _vec) { ConformalPoint=_vec;}
		void SetHoughLine(TF1& _func_ur, TF1& _func_ul, TF1& _func_dr, TF1& _func_dl) 
         { HoughLine_ur=_func_ur; HoughLine_ul=_func_ul; HoughLine_dr=_func_dr; HoughLine_dl=_func_dl; }
		void SetSlayer(int _slayer) { Slayer=_slayer; }

  private: 
	  int Slayer;
    double cellSize; //crystal cell size, unit: mm. 
    CRDEcalEDM::CRDCaloBarShower originLocalMax;  //Local max
		TVector2 ConformalPoint; //Center position. 
		TF1 HoughLine_ur;
		TF1 HoughLine_ul;
		TF1 HoughLine_dr;
		TF1 HoughLine_dl;

};

};
#endif
