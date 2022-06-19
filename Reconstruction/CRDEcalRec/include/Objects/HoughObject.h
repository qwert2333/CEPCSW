#ifndef CALOHOUGHOBJECT_H
#define CALOHOUGHOBJECT_H
#include "Objects/CaloBarShower.h"
#include "TVector2.h"
#include "TF1.h"
namespace PandoraPlus {

  class HoughObject{

  public: 
    HoughObject() {};
    ~HoughObject() { Clear(); };

    void Clean() { delete originLocalMax; }
    void Clear() { originLocalMax=nullptr; ConformalPoint.SetX(0.); ConformalPoint.SetY(0.);}

    //inline bool operator == (const HoughObject &x) const{
    //  return  originLocalMax == x.originLocalMax;
    //}

		const PandoraPlus::CaloBarShower* getLocalMax() const { return originLocalMax; }
		TVector2 getConformPointUR() const { return ConformalPoint+TVector2(cellSize/2., cellSize/2.);   }
		TVector2 getConformPointUL() const { return ConformalPoint+TVector2(-cellSize/2., cellSize/2.);  }
		TVector2 getConformPointDR() const { return ConformalPoint+TVector2(cellSize/2., -cellSize/2.);  }
		TVector2 getConformPointDL() const { return ConformalPoint+TVector2(-cellSize/2., -cellSize/2.); }
    TF1 getHoughLineUR() const { return HoughLine_ur; }
    TF1 getHoughLineUL() const { return HoughLine_ul; }
    TF1 getHoughLineDR() const { return HoughLine_dr; }
    TF1 getHoughLineDL() const { return HoughLine_dl; }

    void SetCellSize(double _cell) { cellSize=_cell; }
    void SetLocalMax( const PandoraPlus::CaloBarShower* _localmax ) { originLocalMax=_localmax; }
    void SetConformalPoint(TVector2& _vec) { ConformalPoint=_vec;}
		void SetHoughLine(TF1& _func_ur, TF1& _func_ul, TF1& _func_dr, TF1& _func_dl) 
         { HoughLine_ur=_func_ur; HoughLine_ul=_func_ul; HoughLine_dr=_func_dr; HoughLine_dl=_func_dl; }
		void SetSlayer(int _slayer) { Slayer=_slayer; }

  private: 
	  int Slayer;
    double cellSize; //crystal cell size, unit: mm. 
    const PandoraPlus::CaloBarShower* originLocalMax;  //Local max
		TVector2 ConformalPoint; //Center position. 
		TF1 HoughLine_ur;
		TF1 HoughLine_ul;
		TF1 HoughLine_dr;
		TF1 HoughLine_dl;

};

};
#endif
