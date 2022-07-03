#ifndef CALO_CLUSTER_C
#define CALO_CLUSTER_C

#include "Objects/CaloCluster.h"
#include "TH1.h"
#include "TF1.h"
#include "TVector3.h"
#include "TCanvas.h"
#include <map>

namespace PandoraPlus{

  void CaloCluster::Clear(){
    showers.clear();
    tracks.clear();
    daughter_clusters.clear();
    hits.clear();
    showerMax=0; chi2=0; alpha=0; beta=0; ExpEn=0; axis.SetXYZ(0,0,0); type=-1;
  }

  void CaloCluster::Clean(){
    for(int i=0; i<showers.size(); i++) { delete showers[i]; showers[i] = NULL; }
    for(int i=0; i<tracks.size(); i++)  { delete tracks[i];  tracks[i] = NULL; }
    for(int i=0; i<daughter_clusters.size(); i++) { delete daughter_clusters[i]; daughter_clusters[i] = NULL; }
    for(int i=0; i<hits.size(); i++) { delete hits[i]; hits[i] = NULL; }
    Clear();
  }

  void CaloCluster::Check(){
    for(int i=0; i<showers.size(); i++)
      if(!showers[i]) { showers.erase(showers.begin()+i); i--; }
    for(int i=0; i<tracks.size(); i++)
      if(!tracks[i]) { tracks.erase(tracks.begin()+i); i--; }
    for(int i=0; i<daughter_clusters.size(); i++)
      if(!daughter_clusters[i]) { daughter_clusters.erase(daughter_clusters.begin()+i); i--; }
    for(int i=0; i<hits.size(); i++)
      if(!hits[i]) { hits.erase(hits.begin()+i); i--; }
  }


  //Get the initial hit in this shower. defined as closest to IP. 
  const PandoraPlus::CaloHit* CaloCluster::getClusterInitialHit() const{
    double minr=10e9;
    PandoraPlus::CaloHit* initHit = nullptr; 
    int index = -1;
    for(int i=0;i<hits.size();i++){
       double rhit = hits[i]->getPosition().Mag();
       if(rhit<minr){ minr = rhit; index=i; }
    }
    if(index>0) return hits[index];
    else return nullptr;
  }
  
  //Get shower center by weighted hit position. 
  TVector3 CaloCluster::getHitCenter() const{
    TVector3 vec(0,0,0);
    double totE = getHitsE();
    for(int i=0;i<hits.size(); i++){
       TVector3 v_cent = hits[i]->getPosition();
       vec += v_cent * (hits[i]->getEnergy()/totE);
    }
    return vec;
  }

  //Get shower center by weighted 2Dshower(shower in layers) position. 
  TVector3 CaloCluster::getShowerCenter() const{
    TVector3 spos(0,0,0);
    double totE = getShowerE();
    for(int i=0;i<showers.size(); i++) spos += showers[i]->getPos()*(showers[i]->getShowerE()/totE);
    return spos;
  }

  //Sum over the hit energy
  double CaloCluster::getHitsE() const{
    double en=0;
    for(int i=0;i<hits.size(); i++) en+=hits[i]->getEnergy();
    return en;
  }

  //Sum over the 2Dshower energy
  double CaloCluster::getShowerE() const{
    double en=0;
    for(int i=0;i<showers.size(); i++) en+=showers[i]->getShowerE();
    return en;
  }


  int CaloCluster::getBeginningDlayer() const{
    int re_dlayer = -99; 
    std::vector<int> dlayers; dlayers.clear(); 
    for(int ish=0; ish<showers.size(); ish++) dlayers.push_back(showers[ish]->getDlayer());
    re_dlayer = *std::min_element(dlayers.begin(), dlayers.end());
    
    return re_dlayer; 
  }

  int CaloCluster::getEndDlayer() const{
    int re_dlayer = -99;
    std::vector<int> dlayers; dlayers.clear();
    for(int ish=0; ish<showers.size(); ish++) dlayers.push_back(showers[ish]->getDlayer());
    re_dlayer = *std::max_element(dlayers.begin(), dlayers.end());

    return re_dlayer;
  }


  std::vector<double> CaloCluster::getEnInLayer() const{
    map<int, double> orderedEn; orderedEn.clear();
    for(int il=0; il<showers.size(); il++){
      map<int, double>::iterator iter = orderedEn.find(showers[il]->getDlayer());
      if(iter==orderedEn.end()) orderedEn[showers[il]->getDlayer()] = showers[il]->getShowerE();
      else orderedEn[showers[il]->getDlayer()] += showers[il]->getShowerE();
    }

    std::vector<double> m_EnInLayer; m_EnInLayer.clear();
    m_EnInLayer.resize(14);  //TODO: hard-coded 14 layers! 
    for(int il=0; il<14; il++){
      if( orderedEn.find(il)==orderedEn.end() )  m_EnInLayer[il]=0;
      else m_EnInLayer[il] = orderedEn[il];
    }

    return m_EnInLayer; 
  }


  int CaloCluster::getMaxELayer() const{

    std::vector<double> m_EnVec = getEnInLayer(); 

    int m_maxL=-1; 
    double maxE = -99; 
    for(int i=0; i<14; i++)
      if(m_EnVec[i]>maxE) { maxE=m_EnVec[i]; m_maxL=i; }
    
    return m_maxL; 
    
  }


  double CaloCluster::getAveE() const{
    double sumE = getShowerE(); 
    std::vector<double> m_EnVec = getEnInLayer();
    int count=0; 
    for(int i=0; i<14; i++)
      if(m_EnVec[i]!=0) count++; 

    return sumE/(double)count; 
  }


  double CaloCluster::getStdDevE() const{
    double aveE = getAveE();     
    std::vector<double> m_EnVec = getEnInLayer();
    double sumE2=0; 
    int count=0; 
    for(int i=0; i<14; i++){   //TODO: Hard-coded 14 layers! 
      if( m_EnVec[i]==0 ) continue; 
      sumE2 += (m_EnVec[i]-aveE)*(m_EnVec[i]-aveE); 
      count++; 
    }

    double StdDev = sqrt( sumE2/(double)count );
    return StdDev; 
  }

  TVector3 CaloCluster::getAxis() {
    if(showers.size()==0) FitAxisHit();
    else FitAxis();
    return axis;
  }


/*
  std::vector<double> CaloCluster::getClusterWidth() const{
    std::vector<double> widthVec; widthVec.clear(); 
    widthVec.resize(14);

    std::map<int, std::vector<PandoraPlus::CRDCaloHitTransShower> > m_orderedShower; m_orderedShower.clear();
    for(int is=0;is<showers.size();is++){
      m_orderedShower[showers[is].getDlayer()].push_back(showers[is]);
    }
    for(int il=0; il<14; il++){
      std::vector<PandoraPlus::CRDCaloHitTransShower> m_showers = m_orderedShower[il];
      if(m_showers.size()==0) { widthVec[il]=0; continue; }
      if(m_showers.size()==1) { widthVec[il]=m_showers[0].getHitsWidth(); continue; }

      std::vector<edm4hep::ConstCalorimeterHit> m_hits; m_hits.clear();
      double totE = 0;
      for(int is=0; is<m_showers.size(); is++){
        std::vector<edm4hep::ConstCalorimeterHit> m_calohits = m_showers[is].gethits();
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
*/

  double CaloCluster::getMaxWidth() const{
    double re_width = -99; 
    //std::vector<double> widthVec = getClusterWidth();    
    //re_width = *std::max_element(widthVec.begin(), widthVec.end());
    return re_width; 
  }

/*
  void CaloCluster::FitProfile(){
    FitAxis();
    edm4hep::ConstCalorimeterHit initHit = getClusterInitialHit();
    TVector3 vec_init( initHit.getPosition().x, initHit.getPosition().y, initHit.getPosition().z );
    TVector3 vec_axis = axis;
    double X0=11.2; //unit: mm

    //TCanvas *c1 = new TCanvas(); 
    //c1->cd(); 
    TH1D *h_hitz = new TH1D("h_hitz", "h_hitz", 14,0,28);
    for(int i=0;i<hits.size();i++){
       TVector3 vec_ihit(hits[i].getPosition().x, hits[i].getPosition().y, hits[i].getPosition().z);
       TVector3 vec_rel = vec_ihit-vec_init;
       h_hitz->Fill(vec_rel.Dot(vec_axis)/X0, hits[i].getEnergy() );
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
    showerMax = (alpha-1)/beta;

    //h_hitz->Draw("same");
    //func->Draw("same"); 
    //c1->SaveAs("/cefs/higgs/guofy/cepcsoft/CEPCSW_v6/run/CheckClusterMerge/plots/Longi.C");
    //delete c1;
    delete h_hitz; 
  }
*/

  void CaloCluster::FitAxis(){

    if(showers.size()==0) axis.SetXYZ(0,0,0);

    else if(showers.size()==1){ 
      axis = showers[0]->getPos();
      axis *= 1./axis.Mag();
    }

    else if( showers.size()==2 ){

      axis = showers[1]->getPos() - showers[0]->getPos() ;
      int sign = (showers[1]->getDlayer()>showers[0]->getDlayer() ? 1 : -1);
      axis *= (double)sign/axis.Mag();
    }

    else if(showers.size()>2){
      trackFitter.clear();
      //track->setImpactParameter(0, 230.); //fix dr=0, dz=230.

      double barAngle = (showers[0]->getModule()+2)*TMath::Pi()/4.;
      double posErr = 10./sqrt(12); 
      if(barAngle>=TMath::TwoPi()) barAngle = barAngle-TMath::TwoPi();
      trackFitter.setBarAngle(barAngle);
      for(int is=0;is<showers.size();is++){
        TVector3 pos_barsX = showers[is]->getShowerX()->getPos(); //U
        TVector3 pos_barsY = showers[is]->getShowerY()->getPos(); //Z
//printf("\t DEBUG: input pointX (%.3f, %.3f, %.3f) \n", barsX.getPos().x(), barsX.getPos().y(), barsX.getPos().z());
//printf("\t DEBUG: input pointY (%.3f, %.3f, %.3f) \n", barsY.getPos().x(), barsY.getPos().y(), barsY.getPos().z());
        trackFitter.setGlobalPoint(1, pos_barsX.x(), posErr, pos_barsX.y(), posErr, pos_barsX.z(), posErr);
        trackFitter.setGlobalPoint(0, pos_barsY.x(), posErr, pos_barsY.y(), posErr, pos_barsY.z(), posErr);
      }
    trackFitter.fitTrack();
    double fitPhi =   trackFitter.getTrkPar(2);
    double fitTheta = trackFitter.getTrkPar(3);
//printf("\t DEBUG: fitted phi and theta: %.3f \t %.3f \n", fitPhi, fitTheta);

    axis.SetMag(1.);
    axis.SetPhi(fitPhi);
    axis.SetTheta(fitTheta);
  } 
  }

  void CaloCluster::FitAxisHit(){
    int Nhit = hits.size(); 
    if(Nhit==0) axis.SetXYZ(0,0,0);
    else if(Nhit==1){
      axis = hits[0]->getPosition();
      axis *= 1./axis.Mag();
    }
    else if( Nhit==2 ){
      axis = hits[1]->getPosition() - hits[0]->getPosition() ;
      int sign = (hits[1]->getLayer()>hits[0]->getLayer() ? 1 : -1);
      axis *= (double)sign/axis.Mag();
    }

    else{
/*      float sumX = 0;
      float sumY = 0;
      float sumZ = 0;
      float sumXZ = 0;
      float sumYZ = 0;
      float sumZZ = 0;
      for(int i=0; i<Nhit; i++){
        TVector3 pos = hits[i]->getPosition();
   
        sumX += pos.x();
        sumY += pos.y();
        sumZ += pos.z();
        sumXZ += pos.x()*pos.z();
        sumYZ += pos.y()*pos.z();
        sumZZ += pos.z()*pos.z();
      }
   
      float k1 = (Nhit*sumXZ - sumX*sumZ)/(Nhit*sumZZ-sumZ*sumZ);
      float k2 = (Nhit*sumYZ - sumY*sumZ)/(Nhit*sumZZ-sumZ*sumZ);
   
      axis.SetXYZ(k1, k2, 1);
      TVector3 pos0 = hits[0]->getPosition();
      int sign = (axis.Dot(pos0)>0  ? 1 : -1  );
      axis *= (double)sign/axis.Mag();
*/
      trackFitter.clear();
      //track->setImpactParameter(0, 230.); //fix dr=0, dz=230.

      //double barAngle = (showers[0]->getModule()+2)*TMath::Pi()/4.;
      //if(barAngle>=TMath::TwoPi()) barAngle = barAngle-TMath::TwoPi();
      //trackFitter.setBarAngle(barAngle);
      trackFitter.setBarAngle(0.);
      double posErr = 10./sqrt(12);
      for(int is=0;is<hits.size();is++){
        TVector3 pos = hits[is]->getPosition(); //U
//printf("\t DEBUG: input pointX (%.3f, %.3f, %.3f) \n", barsX.getPos().x(), barsX.getPos().y(), barsX.getPos().z());
//printf("\t DEBUG: input pointY (%.3f, %.3f, %.3f) \n", barsY.getPos().x(), barsY.getPos().y(), barsY.getPos().z());
        trackFitter.setGlobalPoint(1, pos.x(), posErr, pos.y(), posErr, pos.z(), posErr);
        trackFitter.setGlobalPoint(0, pos.x(), posErr, pos.y(), posErr, pos.z(), posErr);
      }
    trackFitter.fitTrack();
    double fitPhi =   trackFitter.getTrkPar(2);
    double fitTheta = trackFitter.getTrkPar(3);
//printf("\t DEBUG: fitted phi and theta: %.3f \t %.3f \n", fitPhi, fitTheta);

    axis.SetMag(1.);
    axis.SetPhi(fitPhi);
    axis.SetTheta(fitTheta);

    }
  }

  void CaloCluster::MergeCluster(const PandoraPlus::CaloCluster* clus){
    for(int i=0;i<clus->getShowers().size();i++){
       addShower(clus->getShowers()[i]);
    }
  }


  void CaloCluster::IdentifyCluster(){
    double stdDevE = getStdDevE(); 

    if( stdDevE<0.05 ) type = 0; 
    else if( stdDevE>0.2 ) type = 1; 
    else type = 2;  

  }


  void CaloCluster::Print() const{
    cout<<"---------Print out Cluster-----------"<<endl;

    cout<<"--Print showers: size = "<<showers.size()<<endl;
    for(int i=0; i<showers.size(); i++){
      cout<<"  Shower #"<<i<<": Address "<<showers[i]<<", Layer = "<<showers[i]->getDlayer(); 
      printf(", pos/E (%.2f, %.2f, %.2f, %.3f) \n", showers[i]->getPos().x(), showers[i]->getPos().y(), showers[i]->getPos().z(), showers[i]->getShowerE() );
      cout<<"  BarShowerX address: "<<showers[i]->getShowerX(); 
      printf(", pos/E (%.2f, %.2f, %.2f, %.3f) \n", showers[i]->getShowerX()->getPos().x(), 
                                                    showers[i]->getShowerX()->getPos().y(), 
                                                    showers[i]->getShowerX()->getPos().x(), 
                                                    showers[i]->getShowerX()->getE() );
      cout<<"  BarShowerY address: "<<showers[i]->getShowerY();
      printf(", pos/E (%.2f, %.2f, %.2f, %.3f) \n", showers[i]->getShowerY()->getPos().x(),
                                                    showers[i]->getShowerY()->getPos().y(),
                                                    showers[i]->getShowerY()->getPos().x(),
                                                    showers[i]->getShowerY()->getE() );
    }

    cout<<"---------End Print Cluster---------"<<endl;
    cout<<endl;

  }



};
#endif
