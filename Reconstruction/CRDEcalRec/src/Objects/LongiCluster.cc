#ifndef LONGICLUS_C
#define LONGICLUS_C

#include "Objects/LongiCluster.h"
#include <cmath>

namespace PandoraPlus{

  void LongiCluster::Clean(){
    for(int i=0; i<barShowerCol.size(); i++) {delete barShowerCol[i]; barShowerCol[i]=NULL; }
    Clear();
  }

  void LongiCluster::Check(){
    for(int i=0; i<barShowerCol.size(); i++) if(!barShowerCol[i]) { barShowerCol.erase(barShowerCol.begin()+i); i--; }
  }

  int LongiCluster::getSlayer() const{
    if(barShowerCol.size()>0) return barShowerCol[0]->getSlayer();
    else return -1;
  }

  TVector3 LongiCluster::getPos() const{
    TVector3 pos(0, 0, 0);
    double Etot = getEnergy();
    for(int i=0; i<barShowerCol.size(); i++){
      TVector3 m_pos(barShowerCol[i]->getPos().x(), barShowerCol[i]->getPos().y(), barShowerCol[i]->getPos().z());
      pos += m_pos * (barShowerCol[i]->getEnergy()/Etot);
    }
    return pos; 
  }

  double LongiCluster::getEnergy() const{
    double totE = 0;
    for(int i=0; i<barShowerCol.size(); i++) totE += barShowerCol[i]->getEnergy();
    return totE;
  }


  std::vector<const PandoraPlus::Calo1DCluster*> LongiCluster::getBarShowersInLayer(int _layer) const{
    std::vector<const PandoraPlus::Calo1DCluster*> outShowers; outShowers.clear(); 
    for(int i=0; i<barShowerCol.size(); i++)
      if(barShowerCol[i]->getDlayer()==_layer) outShowers.push_back(barShowerCol[i]);
    return outShowers; 
  }


  int LongiCluster::getBeginningDlayer() const{
    int Lstart = 99; 
    for(int i=0; i<barShowerCol.size(); i++)
      if(barShowerCol[i]->getDlayer()<Lstart) Lstart = barShowerCol[i]->getDlayer(); 
    if(Lstart==99) return -99;
    else return Lstart;  
  }


  int LongiCluster::getEndDlayer() const{
    int Lend = -99;
    for(int i=0; i<barShowerCol.size(); i++)
      if(barShowerCol[i]->getDlayer()>Lend) Lend = barShowerCol[i]->getDlayer(); 
    return Lend;
  }


  bool LongiCluster::isContinue() const{
    int ly_max = -1; 
    int ly_min = 99;
    std::vector<int> vec_layers; vec_layers.clear();
    for(int il=0; il<barShowerCol.size(); il++){
      vec_layers.push_back(barShowerCol[il]->getDlayer());
      if(barShowerCol[il]->getDlayer()>ly_max) ly_max = barShowerCol[il]->getDlayer();
      if(barShowerCol[il]->getDlayer()<ly_min) ly_min = barShowerCol[il]->getDlayer();
    }

    bool flag = true;
    std::sort(vec_layers.begin(), vec_layers.end());
    for(int il=0; il<vec_layers.size() && vec_layers[il]!=ly_max; il++)
      if( find(vec_layers.begin(), vec_layers.end(), vec_layers[il]+1) == vec_layers.end() ){ flag=false; break; }

    return flag;
  }

  bool LongiCluster::isContinueN(int n) const{
    if(n<=0) return true; 

    int ly_max = -1;
    int ly_min = 99;
    std::vector<int> vec_layers; vec_layers.clear();
    for(int il=0; il<barShowerCol.size(); il++){
      vec_layers.push_back(barShowerCol[il]->getDlayer());
      if(barShowerCol[il]->getDlayer()>ly_max) ly_max = barShowerCol[il]->getDlayer();
      if(barShowerCol[il]->getDlayer()<ly_min) ly_min = barShowerCol[il]->getDlayer();
    }
    if(n>(ly_max-ly_min+1)) return false; 

    bool flag = false;
    std::sort(vec_layers.begin(), vec_layers.end());
    for(int il=0; il<vec_layers.size(); il++){
      bool fl_continueN = true;
      for(int in=1; in<n; in++) 
        if( find(vec_layers.begin(), vec_layers.end(), vec_layers[il]+in)==vec_layers.end() ) { fl_continueN=false; break; }

      if(fl_continueN) {flag = true; break;}
    }
    return flag;
  }

  bool LongiCluster::isSubset(const LongiCluster* clus) const{

    for(int is=0; is<clus->getBarShowers().size(); is++)
      if( find(barShowerCol.begin(), barShowerCol.end(), clus->getBarShowers()[is])==barShowerCol.end() ) {return false; }
    return true;
  }

  double LongiCluster::OverlapRatioE( const LongiCluster* clus) const{
    double Eshare = 0.;
    for(int is=0; is<clus->getBarShowers().size(); is++)
      if( find(barShowerCol.begin(), barShowerCol.end(), clus->getBarShowers()[is])!=barShowerCol.end() ) Eshare += clus->getBarShowers()[is]->getEnergy(); 

    return Eshare / getEnergy() ;
  }

  
  void LongiCluster::FitAxis(){
    if(barShowerCol.size()==0){ axis.SetXYZ(0,0,0); return; }
    else if(barShowerCol.size()==1){ 
      axis.SetXYZ(barShowerCol[0]->getPos().x(), barShowerCol[0]->getPos().y(), barShowerCol[0]->getPos().z());
      axis = axis.Unit(); 
      return;
    }
    else if(barShowerCol.size()==2){
      TVector3 rpos = barShowerCol.back()->getPos() - barShowerCol.front()->getPos();
      axis.SetXYZ( rpos.x(), rpos.y(), rpos.z() );
      axis = axis.Unit();
      return; 
    }
    else{
      track->clear(); 
      double barAngle = (barShowerCol[0]->getTowerID()[0][0]+2)*TMath::Pi()/4.;
      double posErr = 10./sqrt(12);
      if(barAngle>=TMath::TwoPi()) barAngle = barAngle-TMath::TwoPi();
      track->setBarAngle(barAngle);
      for(int is=0; is<barShowerCol.size(); is++){
        TVector3 b_pos = barShowerCol[is]->getPos(); 
        track->setGlobalPoint(0, b_pos.x(), posErr, b_pos.y(), posErr, b_pos.z(), posErr);
        track->setGlobalPoint(1, b_pos.x(), posErr, b_pos.y(), posErr, b_pos.z(), posErr);
      }

      track->fitTrack();
      double fitPhi = track->getPhi();
      double fitTheta = track->getTheta();

      trk_dr = track->getDr();
      trk_dz = track->getDz();
      axis.SetMag(1.); 
      axis.SetPhi(fitPhi);
      axis.SetTheta(fitTheta);
    }
  }
  

  void LongiCluster::addBarShower( const PandoraPlus::Calo1DCluster* _shower, int option ){
    int index = -1;
    for(int i=0; i<barShowerCol.size(); i++)
      if(_shower->getDlayer()==barShowerCol[i]->getDlayer() ){ index=i; break; }
    
    if(index<0 || option==0){
      barShowerCol.push_back( _shower );
      //FitAxis();
    }
    else{
      PandoraPlus::Calo1DCluster* newshower = new PandoraPlus::Calo1DCluster();
      newshower->setBars( barShowerCol[index]->getBars() );
      newshower->setSeeds( barShowerCol[index]->getSeeds() );
      for(int i=0; i<_shower->getBars().size(); i++) newshower->addUnit( _shower->getBars()[i] );

      barShowerCol.erase( barShowerCol.begin()+index );  //WARNING: potential memory leakage! 
      barShowerCol.push_back( newshower );
    }
  }

  void LongiCluster::MergeCluster(const LongiCluster* clus ){
    for(int is=0; is<clus->getBarShowers().size(); is++){
      if( find(barShowerCol.begin(), barShowerCol.end(), clus->getBarShowers()[is])==barShowerCol.end() ) 
        barShowerCol.push_back( clus->getBarShowers()[is] );
    }
    //FitAxis();
  }

  void LongiCluster::RemoveShowers( std::vector<const PandoraPlus::Calo1DCluster*>& _showers ){
    for(int is=0; is<_showers.size(); is++){
      auto iter = std::remove(barShowerCol.begin(), barShowerCol.end(), _showers[is]);
      barShowerCol.erase( iter, barShowerCol.end() );
    }
  }


};
#endif
