#ifndef _TRACK_H
#define _TRACK_H
#include "TVector3.h"

namespace CRDEcalEDM {

  class Track{
  public:
    Track() {};
    ~Track() {};
  void Clear(); 

  double getD0()         const { return m_d0; }
  double getZ0()         const { return m_z0; }
  double getPhi0()       const { return m_phi0; }
  double getKappa()      const { return m_kappa; }
  double getTanLambda()  const { return m_tanLambda; }
  int    getPID()        const { return m_pdgID; }
  int    getCharge()     const { return m_charge; }
  TVector3 getMomentum() const {return m_momentum; }
  TVector3 getVertex()   const { return m_vtx; }
  TVector3 getProspectPos(double _var) const; 
  TVector3 getProspectPos(int _dlayer) const; 

  void setHelix(double _d0, double _z0, double _phi0, double _kappa, double _tanL ){
    m_phi0=_phi0; m_d0=_d0; m_z0=_z0; m_kappa=_kappa; m_tanLambda=_tanL;
  }
  void setVertex( TVector3 _vtx ) { m_vtx=_vtx; }
  void setVertex( double _x, double _y, double _z) { m_vtx.SetXYZ(_x, _y, _z); }
  void setMomentum( TVector3 _p) { m_momentum = _p; }
  void setPID(int _pid) { m_pdgID=_pid; }
  void setCharge(int _ch) { m_charge=_ch; }

  private:
    double m_phi0;
    double m_d0;
    double m_z0;
    double m_kappa;
    double m_tanLambda;
    int    m_charge; 
    TVector3 m_momentum;
    TVector3 m_vtx; 

    int m_pdgID;
  };

};
#endif
