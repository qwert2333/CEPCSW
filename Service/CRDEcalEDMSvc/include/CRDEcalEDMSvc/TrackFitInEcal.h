#ifndef TRACKFITINECAL_H
#define TRACKFITINECAL_H

#include <iostream>
#include <vector>

#include "Rtypes.h"

class TrackFitInEcal{
public:
	 TrackFitInEcal(double barAngle=0.);
	 ~TrackFitInEcal();

    void setBarAngle(double barAngle) { m_barAngle=barAngle; }
	 void setPoint(int flagUZ, double uzPos, double uzPosErr, double depth, double depthErr);
    void setGlobalPoint(int flagUZ, double x, double xerr, double y, double yerr, double z, double zerr);
	 bool fitTrack();
	 void clear();

	 bool fit2D();				/* fit in x-y and w-z respectively */
	 bool mnFit3D();			/* fit with minuit */

	 double getChisquare() const{return m_chisq;}
	 double getTrkPar(int i) const;
	 double getTrkParErr(int i) const{return m_trkParErr[i];}
	 double getCovariance(int i, int k) const{return m_covariance[i][k];}

	 double getDr() const{return m_trkPar[0];}
	 double getDz() const{return m_trkPar[1];}
	 double getPhi() const{return m_trkPar[2]+m_barAngle;}
	 double getTheta() const{return m_trkPar[3];}

	 static const int NTRKPAR = 4;

	 static void fcnTrk(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

	 static std::vector<int> m_flagUZ; /* direction in transverse plane, 0 for u, 1 for z */
	 static std::vector<double> m_uzPos; /* transverse position */
	 static std::vector<double> m_uzPosErr; /* transverse position error */
	 static std::vector<double> m_depth; /* longitudinal position */
	 static std::vector<double> m_depthErr; /* longitudinal position error */

private:
	 double m_barAngle;
	 double m_chisq;
	 double m_trkPar[NTRKPAR];		/* dr, dz, phi, theta */
	 double m_trkParErr[NTRKPAR];
	 double m_covariance[NTRKPAR][NTRKPAR];
};

inline double TrackFitInEcal::getTrkPar(int i) const{
	 if(2==i) return m_trkPar[i]+m_barAngle;
	 else return m_trkPar[i];
}

#endif
