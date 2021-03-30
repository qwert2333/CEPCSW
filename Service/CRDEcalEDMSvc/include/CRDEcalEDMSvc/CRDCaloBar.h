#ifndef _CRD_CALOBAR_
#define _CRD_CALOBAR_

#include <DD4hep/Objects.h>
#include "TVector3.h"
using namespace dd4hep;

namespace CRDEcalEDM{

  class CRDCaloBar{

  public:
    CRDCaloBar(unsigned long long _cellID, int _system, int _module, int _dlayer, int _part, int _block, int _slayer, int _bar, dd4hep::Position _pos, double _Q1, double _Q2, double _T1, double _T2)
    : cellID(_cellID), system(_system), module(_module), dlayer(_dlayer), part(_part), block(_block), slayer(_slayer), bar(_bar), position(_pos), Q1(_Q1), Q2(_Q2), T1(_T1), T2(_T2) {}; 
    CRDCaloBar() {};

    inline bool operator < (const CRDCaloBar &x) const {
      return bar<x.bar ;
    }
    inline bool operator == (const CRDCaloBar &x) const{
      return cellID == x.cellID;
    }
    unsigned long long getcellID() const { return cellID; }
    int getSystem() const { return system; }
    int getModule() const { return module; }
    int getDlayer() const { return dlayer; }
    int getPart()   const { return part;   }
    int getBlock()  const { return block;  }
    int getSlayer() const { return slayer; }
    int getBar()    const { return bar;    }
    double getQ1()  const { return Q1;     }
    double getQ2()  const { return Q2;     }
    double getT1()  const { return T1;     }
    double getT2()  const { return T2;     }

    dd4hep::Position getPosition() const { return position; }
    double getEnergy() const { return (Q1+Q2)/2.; }
    bool isNeighbor(CRDCaloBar &x);

    void setcellID(unsigned long long _cellid) { cellID = _cellid; }
    void setcellID(int _system, int _module, int _dlayer, int _part, int _block, int _slayer, int _bar) { system=_system; module=_module; dlayer=_dlayer; part=_part; block=_block; slayer=_slayer; bar=_bar; }
    void setPosition( dd4hep::Position pos) { position = pos; }
    void setPosition( TVector3 posv3) { position.SetXYZ( posv3.x(), posv3.y(), posv3.z() ); }
    void setQ(double _q1, double _q2) { Q1=_q1; Q2=_q2; }
    void setT(double _t1, double _t2) { T1=_t1; T2=_t2; }

  private:
		unsigned long long cellID;
		int system;
		int module;
		int dlayer;
		int part;
		int block;
		int slayer;
		int bar;
		dd4hep::Position position;
		double Q1;      // Q in left readout
		double Q2;      // Q in right readout;
		double T1;    // T in left readout;
		double T2;    // T in right readout;

  };
  typedef  std::vector<CRDEcalEDM::CRDCaloBar> DigiBlock;
  
}
#endif
