#ifndef _CRD_ECAL_3DSHOWER_C
#define _CRD_ECAL_3DSHOWER_C

#include "Objects/CRDCaloHit3DCluster.h"
#include "Objects/TrackFitInEcal.h"
#include "TH1.h"
#include "TF1.h"
#include "TVector3.h"
#include "TCanvas.h"
#include <map>

namespace CRDEcalEDM{


 CRDCaloHit3DCluster:: CRDCaloHit3DCluster(){
    type=-1; 

/*    ID_BDT=-99;

    reader->AddVariable( "Lend", &Lend );
    reader->AddVariable( "LmaxWidth", &LmaxWidth );
    reader->AddVariable( "LmaxE", &LmaxE );
    reader->AddVariable( "maxE", &maxEnergy );
    reader->AddVariable( "maxWidth", &maxWidth );
    reader->AddVariable( "alpha", &alpha );
    reader->AddVariable( "beta", &beta );
    reader->AddVariable( "chi2", &chi2 );
    reader->BookMVA("BDT","/cefs/higgs/guofy/cepcsoft/CEPCSW_v9/Reconstruction/CRDEcalRec/include/Objects/TMVA_CID_BDTG.weights.xml");
*/
  }

  //Get the initial hit in this shower. defined as closest to IP. 
  edm4hep::ConstCalorimeterHit CRDCaloHit3DCluster::getClusterInitialHit() const{
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
  TVector3 CRDCaloHit3DCluster::getHitCenter() const{
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
  TVector3 CRDCaloHit3DCluster::getShowerCenter() const{
    dd4hep::Position spos(0,0,0);
    double totE = getShowerE();
    for(int i=0;i<ShowerinLayer.size(); i++) spos += ShowerinLayer[i].getPos()*ShowerinLayer[i].getShowerE()/totE;
    TVector3 vec(spos.x(), spos.y(), spos.z());
    return vec;
  }

  //Sum over the hit energy
  double CRDCaloHit3DCluster::getHitsE() const{
    double en=0;
    for(int i=0;i<CaloHits.size(); i++) en+=CaloHits[i].getEnergy();
    return en;
  }

  //Sum over the 2Dshower energy
  double CRDCaloHit3DCluster::getShowerE() const{
    double en=0;
    for(int i=0;i<ShowerinLayer.size(); i++) en+=ShowerinLayer[i].getShowerE();
    return en;
  }


  //Get the expected deposited energy in layer i. Now is the average of i+1 and i-1. Need to update with profile. 
  double CRDCaloHit3DCluster::getExpEnergy(int& dlayer) const{ 
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
  TVector3 CRDCaloHit3DCluster::getExpPos(int& dlayer) const{
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

  int CRDCaloHit3DCluster::getBeginningDlayer() const{
    int re_dlayer = -99; 
    std::vector<int> dlayers; dlayers.clear(); 
    for(int ish=0; ish<ShowerinLayer.size(); ish++) dlayers.push_back(ShowerinLayer[ish].getDlayer());
    re_dlayer = *std::min_element(dlayers.begin(), dlayers.end());
    
    return re_dlayer; 
  }

  int CRDCaloHit3DCluster::getEndDlayer() const{
    int re_dlayer = -99;
    std::vector<int> dlayers; dlayers.clear();
    for(int ish=0; ish<ShowerinLayer.size(); ish++) dlayers.push_back(ShowerinLayer[ish].getDlayer());
    re_dlayer = *std::max_element(dlayers.begin(), dlayers.end());

    return re_dlayer;
  }


  std::vector<double> CRDCaloHit3DCluster::getEnInLayer() const{
    map<int, double> orderedEn; orderedEn.clear();
    for(int il=0; il<ShowerinLayer.size(); il++){
      map<int, double>::iterator iter = orderedEn.find(ShowerinLayer[il].getDlayer());
      if(iter==orderedEn.end()) orderedEn[ShowerinLayer[il].getDlayer()] = ShowerinLayer[il].getShowerE();
      else orderedEn[ShowerinLayer[il].getDlayer()] += ShowerinLayer[il].getShowerE();
    }

    std::vector<double> m_EnInLayer; m_EnInLayer.clear();
    m_EnInLayer.resize(14); 
    for(int il=0; il<14; il++){
      if( orderedEn.find(il)==orderedEn.end() )  m_EnInLayer[il]=0;
      else m_EnInLayer[il] = orderedEn[il];
    }

    return m_EnInLayer; 
  }


  int CRDCaloHit3DCluster::getMaxELayer() const{

    std::vector<double> m_EnVec = getEnInLayer(); 

    int m_maxL=-1; 
    double maxE = -99; 
    for(int i=0; i<14; i++)
      if(m_EnVec[i]>maxE) { maxE=m_EnVec[i]; m_maxL=i; }
    
    return m_maxL; 
    
  }


  double CRDCaloHit3DCluster::getAveE() const{
    double sumE = getShowerE(); 
    std::vector<double> m_EnVec = getEnInLayer();
    int count=0; 
    for(int i=0; i<14; i++)
      if(m_EnVec[i]!=0) count++; 

    return sumE/(double)count; 
  }


  double CRDCaloHit3DCluster::getStdDevE() const{
    double aveE = getAveE();     
    std::vector<double> m_EnVec = getEnInLayer();
    double sumE2=0; 
    int count=0; 
    for(int i=0; i<14; i++){
      if( m_EnVec[i]==0 ) continue; 
      sumE2 += (m_EnVec[i]-aveE)*(m_EnVec[i]-aveE); 
      count++; 
    }

    double StdDev = sqrt( sumE2/(double)count );
    return StdDev; 
  }


  std::vector<double> CRDCaloHit3DCluster::getClusterWidth() const{
    std::vector<double> widthVec; widthVec.clear(); 
    widthVec.resize(14);

    std::map<int, std::vector<CRDEcalEDM::CRDCaloHit2DShower> > m_orderedShower; m_orderedShower.clear();
    for(int is=0;is<ShowerinLayer.size();is++){
      m_orderedShower[ShowerinLayer[is].getDlayer()].push_back(ShowerinLayer[is]);
    }
    for(int il=0; il<14; il++){
      std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_showers = m_orderedShower[il];
      if(m_showers.size()==0) { widthVec[il]=0; continue; }
      if(m_showers.size()==1) { widthVec[il]=m_showers[0].getHitsWidth(); continue; }

      std::vector<edm4hep::ConstCalorimeterHit> m_hits; m_hits.clear();
      double totE = 0;
      for(int is=0; is<m_showers.size(); is++){
        std::vector<edm4hep::ConstCalorimeterHit> m_calohits = m_showers[is].getCaloHits();
        m_hits.insert(m_hits.end(), m_calohits.begin(), m_calohits.end() );
        totE += m_showers[is].getHitsE();
      }

      TVector3 cent(0., 0., 0.);
      for(int i=0;i<m_hits.size();i++){
        TVector3 pos(m_hits[i].getPosition().x, m_hits[i].getPosition().y, m_hits[i].getPosition().z);
        cent += pos* (m_hits[i].getEnergy()/totE);
      }

      double width=0;
      for(int i=0;i<m_hits.size();i++){
        TVector3 pos(m_hits[i].getPosition().x, m_hits[i].getPosition().y, m_hits[i].getPosition().z);
        double r2 = (pos-cent).Mag2();
        width += r2*m_hits[i].getEnergy()/totE;
      }
      widthVec[il] = width;
    }    

    return widthVec; 
  }


  double CRDCaloHit3DCluster::getMaxWidth() const{
    double re_width = -99; 
    std::vector<double> widthVec = getClusterWidth();    
    re_width = *std::max_element(widthVec.begin(), widthVec.end());
    return re_width; 
  }


  void CRDCaloHit3DCluster::FitProfile(){
    FitAxis();
    edm4hep::ConstCalorimeterHit initHit = getClusterInitialHit();
    TVector3 vec_init( initHit.getPosition().x, initHit.getPosition().y, initHit.getPosition().z );
    TVector3 vec_axis = axis;
    double X0=11.2; //unit: mm

    //TCanvas *c1 = new TCanvas(); 
    //c1->cd(); 
    TH1D *h_hitz = new TH1D("h_hitz", "h_hitz", 14,0,28);
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
    func->SetParameter(2, 10);
//    func->SetParLimits(0, 2, 7);
//    func->SetParLimits(1, 0.3, 0.7);
    h_hitz->Fit("fc1","Q","",0,25);

    chi2 = func->GetChisquare()/func->GetNDF();
    alpha = func->GetParameter(0);
    beta = func->GetParameter(1);
    ExpEn = func->GetParameter(2);

    //h_hitz->Draw("same");
    //func->Draw("same"); 
    //c1->SaveAs("/cefs/higgs/guofy/cepcsoft/CEPCSW_v6/run/CheckClusterMerge/plots/Longi.C");
    //delete c1;
    delete h_hitz; 
  }


  void CRDCaloHit3DCluster::FitAxis(){
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

  void CRDCaloHit3DCluster::AddShower(CRDEcalEDM::CRDCaloHit2DShower& shower){
    ShowerinLayer.push_back(shower);
    std::vector<edm4hep::ConstCalorimeterHit> m_hits = shower.getCaloHits();
    CaloHits.insert(CaloHits.end(), m_hits.begin(), m_hits.end());
    showerEnergy += shower.getShowerE();
    hitsEnergy += shower.getHitsE();
    FitAxis();
  }

  void CRDCaloHit3DCluster::MergeCluster(CRDEcalEDM::CRDCaloHit3DCluster& clus){
    for(int i=0;i<clus.ShowerinLayer.size();i++){
       AddShower(clus.ShowerinLayer[i]);
    }
  }


  void CRDCaloHit3DCluster::IdentifyCluster(){
    float Lstart=0; 
    bool f_found = false;
    double stdDevE = getStdDevE();
    std::vector<double> m_EnVec = getEnInLayer();
    for(int i=0; i<m_EnVec.size(); i++){
      if(!f_found && m_EnVec[i]>0.1) { Lstart=i; f_found==true; }
      if(m_EnVec[i]>maxEnergy) { maxEnergy=m_EnVec[i]; LmaxE=(float)i; }
    }

    if(Lstart==0 && Lend>12) {type = 0; return; }  //MIP
    else if(stdDevE>0.2) {type = 1; return; }  //EM
    else {type = 2; return; }  //Had

/*    FitProfile(); 
    std::vector<double> m_widthVec = getClusterWidth();
    std::vector<double> m_EnVec = getEnInLayer();

    float Lstart=0; 
    bool f_found = false;
    LmaxWidth=0;
    LmaxE=0;
    maxEnergy=-99;
    maxWidth=-99;
    for(int i=0; i<m_EnVec.size(); i++){
      if(!f_found && m_EnVec[i]>0.1) { Lstart=i; f_found==true; }
      if(m_EnVec[i]>maxEnergy) { maxEnergy=m_EnVec[i]; LmaxE=(float)i; }
    }
    for(int i=0; i<m_widthVec.size(); i++)
      if(m_widthVec[i]>maxWidth) { maxWidth=m_widthVec[i]; LmaxWidth=(float)i; }

    Lend = getEndDlayer();
    alpha = getFitAlpha();
    beta = getFitBeta();
    chi2 = 1.;


    //TMVA::Reader *reader = new TMVA::Reader();
    //reader->AddVariable( "Lend", &Lend );
    //reader->AddVariable( "LmaxWidth", &LmaxWidth );
    //reader->AddVariable( "LmaxE", &LmaxE );
    //reader->AddVariable( "maxE", &maxEnergy );
    //reader->AddVariable( "maxWidth", &maxWidth );
    //reader->AddVariable( "alpha", &alpha );
    //reader->AddVariable( "beta", &beta );
    //reader->AddVariable( "chi2", &chi2 );

    //reader->BookMVA("BDT","/cefs/higgs/guofy/cepcsoft/CEPCSW_v9/Reconstruction/CRDEcalRec/include/Objects/TMVA_CID_BDTG.weights.xml");

    ID_BDT = reader->EvaluateMVA("BDT");

    if(ID_BDT>-0.04) type = 1; //Gam
    else type = 2;   //Had

    //delete reader; 
*/
  }


};
#endif
