#ifndef _TRACK_H
#define _TRACK_H

namespace CRDEcalEDM {

  class Track{
  public:
    Track() {};
    ~Track() {};

  double getD0() const        { return m_d0; }
  double getZ0() const        { return m_z0; }
  double getPhi0() const      { return m_phi0; }
  double getKappa() const     { return m_kappa; }
  double getTanLambda() const { return m_tanLambda; }
  int    getPID() const { return m_pdgID; }
  TVector3 getMomentum(); 


  void setHelix(double _d0, double _z0, double _phi0, double _kappa, double _tanL ){
    m_phi0=_phi0; m_d0=_d0; m_z0=_z0; m_kappa=_kappa; m_tanLambda=_tanL;
  }
  void setPID(int _pid) { m_pdgID=_pid; }

  private:
    double m_phi0;
    double m_d0;
    double m_z0;
    double m_kappa;
    double m_tanLambda;

    int m_pdgID;
  };

};
#endif
