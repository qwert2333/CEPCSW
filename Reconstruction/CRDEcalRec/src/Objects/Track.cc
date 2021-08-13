#ifndef _TRACK_C
#define _TRACK_C

#include "Objects/Track.h"
#include <cmath>

namespace CRDEcalEDM{

  void Track::Clear(){
    m_phi0  = -1;
    m_d0    = -1;
    m_z0    = -1; 
    m_kappa = -1;
    m_tanLambda = -1; 
    m_momentum.SetXYZ(0, 0, 0);
    m_vtx.SetXYZ(0,0,0);
    m_pdgID = 0;
  }

  TVector3 Track::getProspectPos(double _var) const{
    TVector3 vec(0, 0, 0);
    double x, y, z;
    //For helix

    //For line
    x = m_d0*cos(m_phi0) - _var*sin(m_phi0);
    y = m_d0*sin(m_phi0) + _var*cos(m_phi0);
    z = m_z0 + _var*m_tanLambda;
    vec.SetXYZ(x, y, z);

    return vec;
  }

  TVector3 Track::getProspectPos(int _dlayer) const{
    //FIXME: Only for test
    double var = (1800.+_dlayer*20.-10. - (m_d0*cos(m_phi0)))/(-sin(m_phi0));
    return getProspectPos(var);
  }

};
#endif
