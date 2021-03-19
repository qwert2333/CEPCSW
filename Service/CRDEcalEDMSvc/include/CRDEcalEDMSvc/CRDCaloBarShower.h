#ifndef _CRD_CALOBARSHOWER_
#define _CRD_CALOBARSHOWER_

#include <vector>
#include "CRDEcalEDMSvc/CRDCaloBar.h"

namespace CRDEcalEDM {

  class CRDCaloBarShower{

  public: 
    CRDCaloBarShower() {}; 

    void Clear() { Energy=0; Bars.clear(); }
    bool isNeighbor(CRDEcalEDM::CRDCaloBar iBar); 
    bool inShower(CRDEcalEDM::CRDCaloBar iBar); 

    dd4hep::Position getPos() const;
    double getE()  const;
    double getT1() const;
    double getT2() const;
    std::vector<CRDEcalEDM::CRDCaloBar> getBars() const { return Bars; }
    CRDEcalEDM::CRDCaloBar getSeed() const { return Seed; };

    void setBars(std::vector<CRDEcalEDM::CRDCaloBar> _bars){ Bars = _bars; }
    void setSeed(CRDEcalEDM::CRDCaloBar _seed ) { Seed = _seed; }
    void setPos(dd4hep::Position _pos) { pos = _pos; }

  private:
		double Energy;
      dd4hep::Position pos;
      std::vector<CRDEcalEDM::CRDCaloBar> Bars;
      CRDEcalEDM::CRDCaloBar Seed;
  };


};
#endif
