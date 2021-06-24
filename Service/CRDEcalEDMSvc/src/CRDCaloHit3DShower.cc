#ifndef _CRD_ECAL_3DSHOWER_C
#define _CRD_ECAL_3DSHOWER_C

#include "CRDEcalEDMSvc/CRDCaloHit3DShower.h"
#include "CRDEcalEDMSvc/TrackFitInEcal.h"
#include "TH1.h"
#include "TF1.h"
#include "TVector3.h"

namespace CRDEcalEDM{


  //Get the initial hit in this shower. defined as closest to IP. 
  edm4hep::ConstCalorimeterHit CRDCaloHit3DShower::getClusterInitialHit() const{
    double minr=10e9;
    edm4hep::ConstCalorimeterHit initHit;
    for(int i=0;i<CaloHits.size();i++){
       edm4hep::Vector3f pos = CaloHits[i].getPosition();
       double rhit = sqrt(pos.x*pos.x+pos.y*pos.y+pos.z*pos.z);
       if(rhit<minr){ minr = rhit; initHit = CaloHits[i];}
    }
    return initHit;
  }
  
  //Get shower center by weighted hit position. 
  TVector3 CRDCaloHit3DShower::getHitCenter() const{
    TVector3 vec(0,0,0);
    double totE = getHitsE();
    for(int i=0;i<CaloHits.size(); i++){
       edm4hep::Vector3f v_cent = CaloHits[i].getPosition();
       TVector3 m_vec(v_cent.x, v_cent.y, v_cent.z);
       vec += m_vec* (CaloHits[i].getEnergy()/totE);
    }
    return vec;
  }

  //Get shower center by weighted 2Dshower(shower in layers) position. 
  TVector3 CRDCaloHit3DShower::getShowerCenter() const{
    dd4hep::Position spos(0,0,0);
    double totE = getShowerE();
    for(int i=0;i<ShowerinLayer.size(); i++) spos += ShowerinLayer[i].getPos()*ShowerinLayer[i].getShowerE()/totE;
    TVector3 vec(spos.x(), spos.y(), spos.z());
    return vec;
  }

  //Sum over the hit energy
  double CRDCaloHit3DShower::getHitsE() const{
    double en=0;
    for(int i=0;i<CaloHits.size(); i++) en+=CaloHits[i].getEnergy();
    return en;
  }

  //Sum over the 2Dshower energy
  double CRDCaloHit3DShower::getShowerE() const{
    double en=0;
    for(int i=0;i<ShowerinLayer.size(); i++) en+=ShowerinLayer[i].getShowerE();
    return en;
  }


  //Get the expected deposited energy in layer i. Now is the average of i+1 and i-1. Need to update with profile. 
  double CRDCaloHit3DShower::getExpEnergy(int& dlayer) const{ 
    if(dlayer<getBeginningDlayer()) return -99; 
    else if(dlayer>getEndDlayer() ) return 1; 

    double Einlayer[2] = {0}; 
    for(int i=0;i<ShowerinLayer.size(); i++){
      if(ShowerinLayer[i].getDlayer()==dlayer-1) Einlayer[0] = ShowerinLayer[i].getShowerE();
      if(ShowerinLayer[i].getDlayer()==dlayer+1) Einlayer[1] = ShowerinLayer[i].getShowerE();
    }
    return (Einlayer[0]+Einlayer[1])/2.; 
  }


  //Get the expected shower position in layer i, with maxE layer as base.  
  TVector3 CRDCaloHit3DShower::getExpPos(int& dlayer) const{
    double maxE = -99; 
    int dlayer_maxE = 0;
    TVector3 pos_maxE(0,0,0);
    if( dlayer<getBeginningDlayer() ) return pos_maxE; 

    for(int i=0;i<ShowerinLayer.size();i++){
      if(ShowerinLayer[i].getShowerE()>maxE){
        maxE = ShowerinLayer[i].getShowerE(); 
        dlayer_maxE = ShowerinLayer[i].getDlayer(); 
        pos_maxE.SetXYZ(ShowerinLayer[i].getPos().x(), ShowerinLayer[i].getPos().y(), ShowerinLayer[i].getPos().z());
    }}
    return pos_maxE + axis*20*(dlayer-dlayer_maxE);
  }

  int CRDCaloHit3DShower::getBeginningDlayer() const{
    int re_dlayer = -99; 
    std::vector<int> dlayers; dlayers.clear(); 
    for(int ish=0; ish<ShowerinLayer.size(); ish++) dlayers.push_back(ShowerinLayer[ish].getDlayer());
    re_dlayer = *std::min_element(dlayers.begin(), dlayers.end());
    
    return re_dlayer; 
  }

  int CRDCaloHit3DShower::getEndDlayer() const{
    int re_dlayer = -99;
    std::vector<int> dlayers; dlayers.clear();
    for(int ish=0; ish<ShowerinLayer.size(); ish++) dlayers.push_back(ShowerinLayer[ish].getDlayer());
    re_dlayer = *std::max_element(dlayers.begin(), dlayers.end());

    return re_dlayer;
  }


  void CRDCaloHit3DShower::FitProfile(){
    FitAxis();
    edm4hep::ConstCalorimeterHit initHit = getClusterInitialHit();
    TVector3 vec_init( initHit.getPosition().x, initHit.getPosition().y, initHit.getPosition().z );
    TVector3 vec_axis = axis;
    double X0=11.2; //unit: mm

    TH1D *h_hitz = new TH1D("h_hitz", "h_hitz", 28,0,28);
    for(int i=0;i<CaloHits.size();i++){
       TVector3 vec_ihit(CaloHits[i].getPosition().x, CaloHits[i].getPosition().y, CaloHits[i].getPosition().z);
       TVector3 vec_rel = vec_ihit-vec_init;
       h_hitz->Fill(vec_rel.Dot(vec_axis)/X0, CaloHits[i].getEnergy() );
    }
    TF1 *func = new TF1("fc1", "[2]*[1]*pow([1]*x, ([0]-1))*exp(-[1]*x)/TMath::Gamma([0])", 0, 25);
    func->SetParName(0, "alpha");
    func->SetParName(1, "beta");
    func->SetParName(2, "E0");
    func->SetParameter(0, 4.8);
    func->SetParameter(1, 0.5);
    func->SetParameter(2, 30);
    h_hitz->Fit("fc1","","",0,25);

    chi2 = func->GetChisquare()/func->GetNDF();
    alpha = func->GetParameter(0);
    beta = func->GetParameter(1);
    showerMax = (alpha-1)/beta;
  }


  void CRDCaloHit3DShower::FitAxis(){
    if(ShowerinLayer.size()==0) axis.SetXYZ(0,0,0);
    else if(ShowerinLayer.size()==1){ 
      axis.SetXYZ(ShowerinLayer[0].getPos().x(), ShowerinLayer[0].getPos().y(),ShowerinLayer[0].getPos().z());
      axis *= 1./axis.Mag();
    }
    else if( ShowerinLayer.size()==2 ){
      dd4hep::Position rpos = ShowerinLayer.back().getPos() - ShowerinLayer.front().getPos() ;
      double mag = sqrt(rpos.Mag2());
      axis.SetXYZ(rpos.x(), rpos.y(), rpos.z());
      axis *= 1./mag;
    }
    else if(ShowerinLayer.size()>2){
      track->clear();
      //track->setImpactParameter(0, 230.); //fix dr=0, dz=230.

      double barAngle = (ShowerinLayer[0].getModule()+2)*TMath::Pi()/4.;
      double posErr = 10./sqrt(12); 
      if(barAngle>=TMath::TwoPi()) barAngle = barAngle-TMath::TwoPi();
      track->setBarAngle(barAngle);
      for(int is=0;is<ShowerinLayer.size();is++){
        CRDEcalEDM::CRDCaloBarShower barsX = ShowerinLayer[is].getShowerX(); //U
        CRDEcalEDM::CRDCaloBarShower barsY = ShowerinLayer[is].getShowerY(); //Z
//printf("\t DEBUG: input pointX (%.3f, %.3f, %.3f) \n", barsX.getPos().x(), barsX.getPos().y(), barsX.getPos().z());
//printf("\t DEBUG: input pointY (%.3f, %.3f, %.3f) \n", barsY.getPos().x(), barsY.getPos().y(), barsY.getPos().z());
        track->setGlobalPoint(1, barsX.getPos().x(), posErr, barsX.getPos().y(), posErr, barsX.getPos().z(), posErr);
        track->setGlobalPoint(0, barsY.getPos().x(), posErr, barsY.getPos().y(), posErr, barsY.getPos().z(), posErr);
      }
    track->fitTrack();
    double fitPhi = track->getTrkPar(2);
    double fitTheta = track->getTrkPar(3);
//printf("\t DEBUG: fitted phi and theta: %.3f \t %.3f \n", fitPhi, fitTheta);

    axis.SetMag(1.);
    axis.SetPhi(fitPhi);
    axis.SetTheta(fitTheta);
  } 
  }

  void CRDCaloHit3DShower::AddShower(CRDEcalEDM::CRDCaloHit2DShower& shower){
    ShowerinLayer.push_back(shower);
    std::vector<edm4hep::ConstCalorimeterHit> m_hits = shower.getCaloHits();
    CaloHits.insert(CaloHits.end(), m_hits.begin(), m_hits.end());
    showerEnergy += shower.getShowerE();
    hitsEnergy += shower.getHitsE();
    FitAxis();
  }

  void CRDCaloHit3DShower::MergeCluster(CRDEcalEDM::CRDCaloHit3DShower& clus){
    for(int i=0;i<clus.ShowerinLayer.size();i++){
       AddShower(clus.ShowerinLayer[i]);
    }
  }


  void CRDCaloHit3DShower::IdentifyCluster(){
    const int Nlayer = ShowerinLayer.size(); 
    double En[Nlayer] = {0};
    double sumE = 0;
    for(int i=0; i<Nlayer; i++){ 
      En[i] = ShowerinLayer[i].getShowerE(); 
      sumE += En[i];
    }
    double aveE = sumE/Nlayer; 
    double sumE2 = 0;
    for(int i=0; i<Nlayer; i++)
      sumE2 += (En[i]-aveE)*(En[i]-aveE);
    double StdDev = sqrt( sumE2/Nlayer );


    if(ShowerinLayer.size()>=10 && fabs((aveE-0.02)/0.02)<0.2 && StdDev<0.4 ) type = 0; 
    else if( chi2<5 && ShowerinLayer.size() >5 ) type = 1; 
    else type = 2;  

  }

};
#endif
