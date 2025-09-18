#ifndef CALOBAR_H
#define CALOBAR_H

#include <DD4hep/Objects.h>
#include "edm4hep/MCParticle.h"
#include "TVector3.h"

namespace Cyber{

  class CaloUnit{

  public:
    CaloUnit(unsigned long long _cellID, int _system, int _module, int _slayer, int _dlayer, int _stave, int _part, int _bar, TVector3 _pos, double _Q1, double _Q2, double _T1, double _T2)
    : cellID(_cellID), system(_system), module(_module), stave(_stave), part(_part), dlayer(_dlayer), slayer(_slayer), bar(_bar), position(_pos), Q1(_Q1), Q2(_Q2), T1(_T1), T2(_T2) {}; 
    CaloUnit() {};
    ~CaloUnit() { Clear(); };
    void Clear() { position.SetXYZ(0.,0.,0.); Q1=-1; Q2=-1; T1=-99; T2=-99; }

    inline bool operator < (const CaloUnit &x) const {
      if(x.cellID==cellID) return false; 

      if(system==System_Barrel){
        if(slayer==0) return ( (stave<x.stave) || (stave==x.stave && bar<x.bar) );
        else{
          if( module==Nmodule && x.module==0 ) return true;
          else if( module==0 && x.module==Nmodule ) return false; 
          else{
            return ( (module<x.module) || (module==x.module && bar<x.bar) );
          }
        }
      }
      else{
        if(slayer==0) return (position.x() < x.getPosition().x());
        else return (position.y() < x.getPosition().y());
      }
    }

    inline bool operator == (const CaloUnit &x) const{
      return ( (cellID == x.cellID) && getEnergy()==x.getEnergy() );
    }
    unsigned long long getcellID() const { return cellID; }
    int getSystem() const { return system; }
    int getModule() const { return module; }
    int getStave()  const { return stave;  }
    int getPart()   const { return part;   }
    int getDlayer() const { return dlayer; }
    int getSlayer() const { return slayer; }
    int getBar()    const { return bar;    }
    double getQ1()  const { return Q1;     }
    double getQ2()  const { return Q2;     }
    double getT1()  const { return T1;     }
    double getT2()  const { return T2;     }
    double getBarLength() const {return barLength; }
    int getNBarInLayer() const { return nBarInLayer; }

    TVector3 getPosition() const { return position; }
    double getEnergy() const { return (Q1+Q2); }
    std::vector< std::pair<edm4hep::MCParticle, float> > getLinkedMCP() const { return MCParticleWeight; }
    edm4hep::MCParticle getLeadingMCP() const;
    float getLeadingMCPweight() const;
    bool isAtLowerEdgePhi() const; 
    bool isAtUpperEdgePhi() const; 
    bool isAtLowerEdgeZ() const; 
    bool isAtUpperEdgeZ() const; 
    bool isNeighbor(const CaloUnit* x) const;
    //bool isModuleAdjacent( const CaloUnit* x ) const;
    bool isLongiNeighbor(const CaloUnit* x) const;
    //bool isLongiModuleAdjacent( const CaloUnit* x ) const;

    void addLinkedMCP( std::pair<edm4hep::MCParticle, float> _pair ) { MCParticleWeight.push_back(_pair); }
    void setLinkedMCP( std::vector<std::pair<edm4hep::MCParticle, float>> _pairVec ) { MCParticleWeight.clear(); MCParticleWeight = _pairVec; }
    void setcellID(unsigned long long _cellid) { cellID = _cellid; }
    void setcellID(int _system, int _module, int _stave, int _part, int _dlayer, int _slayer, int _bar) { system=_system; module=_module; stave=_stave; part=_part; dlayer=_dlayer; slayer=_slayer; bar=_bar; }
    void setPosition( TVector3 posv3) { position.SetXYZ( posv3.x(), posv3.y(), posv3.z() ); }
    void setQ(double _q1, double _q2) { Q1=_q1; Q2=_q2; }
    void setT(double _t1, double _t2) { T1=_t1; T2=_t2; }
    void setBarLength(double _barLength) { barLength=_barLength; }
    void setNBarInLayer(int _nBarInLayer) { nBarInLayer=_nBarInLayer; }
    std::shared_ptr<CaloUnit> Clone() const;

    static int System_Barrel;
    static int System_Endcap;
    static int Nmodule;
    static int Nstave;
    static int Nlayer;
    static int NbarPhi_odd[7];
    static int NbarPhi_even[7];
    static int NbarZ;
    static float barsize;
    static float ecal_innerR;
    static float ecal_endcap_deadarea;
    static float ecal_endcap_barsize;

  private:
		unsigned long long cellID;
		int system;
		int module;
    int part;
		int stave;
		int dlayer;
		int slayer;
		int bar;
		TVector3 position;
		double Q1;      // Q in left readout
		double Q2;      // Q in right readout;
		double T1;    // T in left readout;
		double T2;    // T in right readout;

    double barLength;
    std::vector< std::pair<edm4hep::MCParticle, float> > MCParticleWeight; 

    int nBarInLayer;
  };
  
}
#endif
