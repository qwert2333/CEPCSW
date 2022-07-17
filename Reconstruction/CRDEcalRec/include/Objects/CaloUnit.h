#ifndef CALOBAR_H
#define CALOBAR_H

#include <DD4hep/Objects.h>
#include "TVector3.h"

namespace PandoraPlus{

  class CaloUnit{

  public:
    CaloUnit(unsigned long long _cellID, int _system, int _module, int _slayer, int _dlayer, int _part, int _stave, int _bar, TVector3 _pos, double _Q1, double _Q2, double _T1, double _T2)
    : cellID(_cellID), system(_system), module(_module), stave(_stave), dlayer(_dlayer), part(_part), slayer(_slayer), bar(_bar), position(_pos), Q1(_Q1), Q2(_Q2), T1(_T1), T2(_T2) {}; 
    CaloUnit() {};
    ~CaloUnit() { Clear(); };
    void Clear() { position.SetXYZ(0.,0.,0.); Q1=-1; Q2=-1; T1=-99; T2=-99; }

    inline bool operator < (const CaloUnit &x) const {
      return ((part<x.part) || (part==x.part && stave<x.stave) ||  (part==x.part && stave==x.stave && bar<x.bar));
    }
    inline bool operator == (const CaloUnit &x) const{
      return ( (cellID == x.cellID) && getEnergy()==x.getEnergy() );
    }
    unsigned long long getcellID() const { return cellID; }
    int getSystem() const { return system; }
    int getModule() const { return module; }
    int getStave()  const { return stave;  }
    int getDlayer() const { return dlayer; }
    int getPart()   const { return part;   }
    int getSlayer() const { return slayer; }
    int getBar()    const { return bar;    }
    double getQ1()  const { return Q1;     }
    double getQ2()  const { return Q2;     }
    double getT1()  const { return T1;     }
    double getT2()  const { return T2;     }

    TVector3 getPosition() const { return position; }
    double getEnergy() const { return (Q1+Q2)/2.; }
    bool isAtLowerEdgePhi() const; 
    bool isAtUpperEdgePhi() const; 
    bool isAtLowerEdgeZ() const; 
    bool isAtUpperEdgeZ() const; 
    bool isNeighbor(const CaloUnit* x) const;
    bool isModuleAdjacent( const CaloUnit* x ) const;

    void setGeomInfo(int _m, int _p, int _s, int _l, int _phi, int _z) { Nmodule=_m; Npart=_p; Nstave=_s; Nlayer=_l; NbarPhi=_phi; NbarZ=_z; }
    void setcellID(unsigned long long _cellid) { cellID = _cellid; }
    void setcellID(int _system, int _module, int _stave, int _dlayer, int _part, int _slayer, int _bar) { system=_system; module=_module; stave=_stave; dlayer=_dlayer; part=_part; slayer=_slayer; bar=_bar; }
    void setPosition( TVector3 posv3) { position.SetXYZ( posv3.x(), posv3.y(), posv3.z() ); }
    void setQ(double _q1, double _q2) { Q1=_q1; Q2=_q2; }
    void setT(double _t1, double _t2) { T1=_t1; T2=_t2; }
    CaloUnit* Clone() const;

  private:
		unsigned long long cellID;
		int system;
		int module;
		int stave;
		int dlayer;
		int part;
		int slayer;
		int bar;
		TVector3 position;
		double Q1;      // Q in left readout
		double Q2;      // Q in right readout;
		double T1;    // T in left readout;
		double T2;    // T in right readout;


    static int Nmodule;
    static int Npart;
    static int Nstave;
    static int Nlayer;
    static int NbarPhi;
    static int NbarZ;

  };
  
}
#endif
