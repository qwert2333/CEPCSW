#ifndef _TRACK_C
#define _TRACK_C

#include "Objects/Track.h"
#include <cmath>
#include "TGraph.h"

namespace CRDEcalEDM{

  const double Track::B = 3.;  //NOTE: hard-coding 3T B-field here!
  const double Track::pos[3] = {0., 0., 0.};
  const double Track::sigma = 5.;
  const double Track::front_face_ecal = 1800.;
  const double Track::PI = 3.141592653;

  void Track::Clear(){ 
    m_phi0  = -1;
    m_d0    = -1;
    m_z0    = -1; 
    m_kappa = -1;
    m_tanLambda = -1; 
    m_momentum.SetXYZ(0, 0, 0);
    m_vtx.SetXYZ(0,0,0);
    m_pdgID = 0;
    extrapolation_point_x[7] = {0};
    extrapolation_point_y[7] = {0};
    extrapolation_point_z[7] = {0};
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

  void Track::setHelix(double p, double thetaini, double phiini)
  {
    double pt = p*sin(thetaini);
    double kap = 1.0/pt;

    double tanl = tan(PI/2. - thetaini);

    double rc = 1000./(0.3*B*kap);
    double phi0Now = phiini - PI/2.;
    double xc = pos[0] - rc*cos(phi0Now);
    double yc = pos[1] - rc*sin(phi0Now);
    double rho = sqrt(xc*xc + yc*yc);
    double dr;
    if(kap<0) dr = rc + rho;
    else dr = rc - rho;

    double rw = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
    double aa = rw*rw - (dr-rc)*(dr-rc) - rc*rc;
    double bb = 2*rc*(dr - rc);
    double cc = aa/bb;
    double dphi = acos(cc);
    double phi0 = phi0Now - dphi;
    if(kap < 0.) phi0 = phi0Now + dphi;
    if(phi0>(2.*PI)) phi0 -= (2.*PI);
    else if(phi0<0.0) phi0 += (2.*PI);

    double dz = rc*tanl*dphi + pos[2];
    if(kap > 0.) dz = pos[2] - rc*tanl*dphi;

    m_phi0 = phi0; 
    m_d0 = dr; 
    m_z0 = dz; 
    m_kappa = kap; 
    m_tanLambda = tanl;
  }

  void Track::setExtrapolation_points()
  {
    extrapolation_point_x[0] = pos[0];
    extrapolation_point_y[0] = pos[1];
    extrapolation_point_z[0] = pos[2];

    double phi0 = m_phi0; 
    double dr = m_d0; 
    double dz = m_z0; 
    double kap = m_kappa; 
    double tanl = m_tanLambda;

    for(int i=1; i<=6; i++)
    {
      double radius = 300*i;
      double g = 1000./(0.3*B*kap);
      double fi0 = phi0 + PI/2.;
      if(fi0 > (PI*2.)) fi0 -= PI*2.; //

      double aaa = radius*radius - (dr-g)*(dr-g) - g*g;
      double bbb = 2*g*(dr - g);
      double ccc = aaa/bbb;
      double ddd = acos(ccc);

      double phi;
      if(kap > 0.) phi = fi0 + ddd;
      else phi = fi0 - ddd;    
      if(phi > (PI*2.)) phi -= PI*2.;  //
      if(phi < 0.) phi += PI*2.;  

      double x0 = dr * cos(phi0);
      double y0 = dr * sin(phi0);
      double x = x0 + g * (sin(phi) - sin(fi0));
      double y = y0 + g * (-cos(phi) + cos(fi0));
      double z;
      if(kap > 0.)  z = dz + g * tanl * ddd;
      else z = dz - g * tanl * ddd;

      extrapolation_point_x[i] = x;
      extrapolation_point_y[i] = y;
      extrapolation_point_z[i] = z;
    }
  }

  void Track::setExtrapolation_front_face()
  {
    TGraph *grxy = new TGraph (7, extrapolation_point_x, extrapolation_point_y);
    TGraph *grxz = new TGraph (7, extrapolation_point_x, extrapolation_point_z);

    extrapolation_front_x = front_face_ecal;
    extrapolation_front_y = grxy->Eval(front_face_ecal);
    extrapolation_front_z = grxz->Eval(front_face_ecal);

    delete grxy;
    delete grxz;
  }
};
#endif
