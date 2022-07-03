#include "Tools/TrackFitInEcal.h"

#include <fstream>
#include <iomanip>
#include <cmath>

#include "TMinuit.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TVector3.h"

using namespace std;

vector<int> TrackFitInEcal::m_flagUZ;
vector<double> TrackFitInEcal::m_uzPos;
vector<double> TrackFitInEcal::m_uzPosErr;
vector<double> TrackFitInEcal::m_depth;
vector<double> TrackFitInEcal::m_depthErr;

TrackFitInEcal::TrackFitInEcal(double barAngle){
	 m_barAngle = barAngle;
    m_IPvar[0]=0; m_IPvar[1]=0; 
    m_FixIP = false; 
}

TrackFitInEcal::~TrackFitInEcal(){

}

void TrackFitInEcal::setPoint(int flagUZ, double uzPos, double uzPosErr, double depth, double depthErr){
	 m_flagUZ.push_back(flagUZ);
	 m_uzPos.push_back(uzPos);
	 m_uzPosErr.push_back(uzPosErr);
	 m_depth.push_back(depth);
	 m_depthErr.push_back(depthErr);
}

void TrackFitInEcal::setGlobalPoint(int flagUZ, double x, double xerr, double y, double yerr, double z, double zerr){
  TVector3 vecp(x, y, z);
  TVector3 vecerr(xerr, yerr, zerr);
  vecp.RotateZ(-1.*m_barAngle);
  vecerr.RotateZ(-1.*m_barAngle);

  if(flagUZ==0) setPoint(flagUZ, vecp.y(), vecerr.y(), vecp.x(), vecerr.x());  //U direction
  else          setPoint(flagUZ, vecp.z(), vecerr.z(), vecp.x(), vecerr.x());  //Z direction
}

bool TrackFitInEcal::fitTrack(){
	 int iXY = 0;
	 int iWZ = 0;
	 for(unsigned int i=0; i<m_flagUZ.size(); i++){
		  int flag = m_flagUZ[i];
		  if(0==flag) iXY++;
		  else if(1==flag) iWZ++;
	 }
	 if((iXY<2) || (iWZ<2)){
		  cout << "WARNING: too few hits for track fit" << endl;
		  return false;
	 }

	 if( ! fit2D() ){
		  cout << "ERROR in 2-D fit" << endl;
		  clear();
		  return false;
	 }

	 if( ! mnFit3D() ){
	 	  clear();
	 	  return false;
	 }

	 clear();
	 return true;
}

bool TrackFitInEcal::fit2D(){
	 int iPointXY = 0;
    // cout<<"Fit IP: "<<m_FixIP<<'\t'<<m_IPvar[0]<<'\t'<<m_IPvar[1]<<endl;

	 TGraphErrors* grXY = new TGraphErrors();
	 for(unsigned int i=0; i<m_flagUZ.size(); i++){
		  int flag = m_flagUZ[i];
		  if(0 != flag) continue;
		  grXY->SetPoint(iPointXY, m_depth[i], m_uzPos[i]);
		  grXY->SetPointError(iPointXY, m_depthErr[i], m_uzPosErr[i]);
		  iPointXY++;
	 }
    TF1 *func1 = new TF1("func1", "pol1", 1800, 2100);
	 grXY->Fit("func1", "Q");
	 double a0 = func1->GetParameter(0);
	 double a1 = func1->GetParameter(1);
    if(m_FixIP){
      func1->FixParameter(0, m_IPvar[0]/cos(atan(a1))); //fix dr. 
      grXY->Fit("func1", "Q");
      a0 = func1->GetParameter(0);
      a1 = func1->GetParameter(1);
    }

	 m_trkPar[2] = atan(a1); // phi
	 m_trkPar[0] = a0 * cos(m_trkPar[2]); // dr
	 double x0 = m_trkPar[0] * cos(m_trkPar[2] + TMath::PiOver2());
	 // cout << "fit in x-y: " << setw(15) << a0 << setw(15) << a1 << setw(15) << m_trkPar[0] << setw(15) << m_trkPar[2] << endl;


	 int iPointWZ = 0;
	 TGraphErrors* grWZ = new TGraphErrors();
	 for(unsigned int i=0; i<m_flagUZ.size(); i++){
		  int flag = m_flagUZ[i];
		  if(1 != flag) continue;

		  double w = (m_depth[i] - x0)/cos(m_trkPar[2]);
		  grWZ->SetPoint(iPointWZ, w, m_uzPos[i]);
		  grWZ->SetPointError(iPointWZ, m_depthErr[i]/cos(m_trkPar[2]), m_uzPosErr[i]);
		  // cout << "iPointWZ " << setw(5) << iPointWZ << setw(15) << w << setw(15) << m_uzPos[i] << endl;
		  iPointWZ++;
	 }
    TF1 *func2 = new TF1("func2", "pol1", 1800, 2100);
    if(m_FixIP) func2->FixParameter(0, m_IPvar[1]); //fix dz.
	 grWZ->Fit("func2", "Q");
	 a0 = func2->GetParameter(0);
	 a1 = func2->GetParameter(1);

	 double lamda = atan(a1);
	 m_trkPar[3] = TMath::PiOver2() - lamda; // theta
	 m_trkPar[1] = a0;						 // dz
	 //cout << "fit in w-z: " << setw(15) << a0 << setw(15) << a1 << setw(15) << m_trkPar[1]
	 // 	  << setw(15) << lamda << setw(15) << m_trkPar[3] << endl;

	 delete grXY;
	 delete grWZ;
    delete func1;
    delete func2;
	 return true;
}

bool TrackFitInEcal::mnFit3D(){
     Int_t ierflg;
     Int_t istat;
     Int_t nvpar;
     Int_t nparx;
     Double_t fmin;
     Double_t edm;
     Double_t errdef;
     Double_t arglist[10];

	 TMinuit* mnTrk = new TMinuit(NTRKPAR);
	 mnTrk->SetPrintLevel(-1);
	 mnTrk->SetFCN(fcnTrk);
	 mnTrk->SetErrorDef(1.0);
	 mnTrk->mnparm(0, "dr", m_IPvar[0], 0.1, 0, 0, ierflg);
	 mnTrk->mnparm(1, "dz", m_IPvar[1], 0.1, 0, 0, ierflg);
	 mnTrk->mnparm(2, "phi", 0, 0.1, 0, 0, ierflg);
	 mnTrk->mnparm(3, "theta", 0, 0.1, 0, 0, ierflg);
	 arglist[0] = 0;
    if(m_FixIP){ mnTrk->FixParameter(0); mnTrk->FixParameter(1); }

	 mnTrk->mnexcm("SET NOW", arglist, 0, ierflg);

	 for(int i=0; i<NTRKPAR; i++){
		  arglist[0] = i + 1;
		  arglist[1] = m_trkPar[i];
		  mnTrk->mnexcm("SET PARameter", arglist, 2, ierflg);
	 }

	 arglist[0] = 1000;
	 arglist[1] = 0.1;
	 mnTrk->mnexcm("MIGRAD", arglist, 2, ierflg);
	 mnTrk->mnstat(fmin, edm, errdef, nvpar, nparx, istat);

     if( (0==ierflg) && (3==istat) ){
          for(int i=0; i<5; i++){
               mnTrk->GetParameter(i, m_trkPar[i], m_trkParErr[i]);
          }

          mnTrk->mnemat(*m_covariance, 4);
          m_chisq = fmin;

          return true;
     } else{
//	 	  cout << "ERROR in fit with minuit! ";
//      cout <<"ierflg = "<<ierflg<<" , stat = "<<istat<<endl;
		  return false;
	 }
}

void TrackFitInEcal::clear(){
	 m_flagUZ.clear();
	 m_uzPos.clear();
	 m_uzPosErr.clear();
	 m_depth.clear();
	 m_depthErr.clear();
    m_IPvar[0]=0; m_IPvar[1]=0;
    m_FixIP = false; 
}

void TrackFitInEcal::fcnTrk(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
     Double_t chisq = 0.0;

	 double fi0 = par[2] + TMath::PiOver2();
	 double x0 = par[0] * cos(fi0);
	 double y0 = par[0] * sin(fi0);
	 for(unsigned int i=0; i<m_flagUZ.size(); i++){
		  if(0==m_flagUZ[i]){
			   double uzPosFit = (m_depth[i]-x0)*tan(par[2]) + y0;
			   Double_t deta = (uzPosFit - m_uzPos[i])/m_uzPosErr[i];
			   chisq += deta*deta;
			   // cout << "fcn x-y: " << setw(15) << uzPosFit << setw(15) << m_uzPos[i] << endl;
		  } else if(1==m_flagUZ[i]){
			   double lamda = TMath::PiOver2() - par[3];
			   double uzPosFit = (m_depth[i]-x0)*tan(lamda)/cos(par[2]) + par[1];
			   Double_t deta = (uzPosFit - m_uzPos[i])/m_uzPosErr[i];
			   chisq += deta*deta;
			   // cout << "fcn w-z: " << setw(15) << uzPosFit << setw(15) << m_uzPos[i] << endl;
		  }
	 }
	 f = chisq;
}

