#ifndef CALOHOUGHOBJECT_H
#define CALOHOUGHOBJECT_H
#include "Objects/Calo1DCluster.h"
#include "TVector2.h"
#include "TF1.h"
namespace PandoraPlus {

  class HoughObject{

  public: 
    HoughObject() {};
    ~HoughObject() { Clear(); };

    void Clean() { for(auto iter:originLocalMax) delete iter;  }
    void Clear() { originLocalMax.clear(); ConformalPoint.SetX(0.); ConformalPoint.SetY(0.);}

    //inline bool operator == (const HoughObject &x) const{
    //  return  originLocalMax == x.originLocalMax;
    //}

		std::vector<const PandoraPlus::Calo1DCluster*> getLocalMax() const { return originLocalMax; }
		TVector2 getConformPointUR() const { return ConformalPoint+TVector2(cellSize/2., cellSize/2.);   }
		TVector2 getConformPointUL() const { return ConformalPoint+TVector2(-cellSize/2., cellSize/2.);  }
		TVector2 getConformPointDR() const { return ConformalPoint+TVector2(cellSize/2., -cellSize/2.);  }
		TVector2 getConformPointDL() const { return ConformalPoint+TVector2(-cellSize/2., -cellSize/2.); }
    TF1 getHoughLineUR() const { return HoughLine_ur; }
    TF1 getHoughLineUL() const { return HoughLine_ul; }
    TF1 getHoughLineDR() const { return HoughLine_dr; }
    TF1 getHoughLineDL() const { return HoughLine_dl; }

    void setCellSize(double _cell) { cellSize=_cell; }
    void addLocalMax( const PandoraPlus::Calo1DCluster* _localmax ) { originLocalMax.push_back(_localmax); }
    void setConformalPoint(TVector2& _vec) { ConformalPoint=_vec;}
		void setHoughLine(TF1& _func_ur, TF1& _func_ul, TF1& _func_dr, TF1& _func_dl) 
         { HoughLine_ur=_func_ur; HoughLine_ul=_func_ul; HoughLine_dr=_func_dr; HoughLine_dl=_func_dl; }
		void setSlayer(int _slayer) { Slayer=_slayer; }

  private: 
	  int Slayer;
    double cellSize; //crystal cell size, unit: mm. 
    std::vector<const PandoraPlus::Calo1DCluster*> originLocalMax;  //Local max
		TVector2 ConformalPoint; //Center position. 
		TF1 HoughLine_ur;
		TF1 HoughLine_ul;
		TF1 HoughLine_dr;
		TF1 HoughLine_dl;

};

};
#endif
