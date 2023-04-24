#ifndef CALO_HALFCLUSTER_C
#define CALO_HALFCLUSTER_C

#include "Objects/CaloHalfCluster.h"

namespace PandoraPlus{

  void CaloHalfCluster::Clear() {
    type = -1; 
    slayer=-99;
    towerID.clear(); 
    m_1dclusters.clear(); 
    map_localMax.clear();
    map_halfClusCol.clear();
  }

  void CaloHalfCluster::Check(){
    for(int i=0; i<m_1dclusters.size(); i++)
      if(!m_1dclusters[i]) { m_1dclusters.erase(m_1dclusters.begin()+i); i--; }
  }

  void CaloHalfCluster::Clean(){
    //for(int i=0; i<m_1dclusters.size(); i++){ delete m_1dclusters[i]; m_1dclusters[i]=NULL; }
    for(auto it: map_localMax) { 
      //for(auto iter: it.second){ delete it.second[iter]; it.second[iter]=NULL; }
      it.second.clear();
    }
    for(auto it: map_halfClusCol) {
      //for(auto iter: it.second){ delete it.second[iter]; it.second[iter]=NULL; }
      it.second.clear();
    }
    Clear();
  }

  std::shared_ptr<PandoraPlus::CaloHalfCluster> CaloHalfCluster::Clone() const{
    std::shared_ptr<PandoraPlus::CaloHalfCluster> m_clus = std::make_shared<PandoraPlus::CaloHalfCluster>();
    for(int i1d=0; i1d<m_1dclusters.size(); i1d++)  m_clus->addUnit(m_1dclusters[i1d]);
    for(int itrk=0; itrk<m_TrackCol.size(); itrk++) m_clus->addAssociatedTrack(m_TrackCol[itrk]);
    for(auto iter:map_localMax)     m_clus->setLocalMax( iter.first, iter.second );
    for(auto iter:map_halfClusCol)  m_clus->setHalfClusters( iter.first, iter.second );
    m_clus->setHoughPars( Hough_alpha, Hough_rho );
    m_clus->setIntercept( Hough_intercept );
    m_clus->setType( type );

    return m_clus;
  }

  bool CaloHalfCluster::isNeighbor(const PandoraPlus::Calo1DCluster* m_1dcluster) const{
    assert(m_1dcluster->getBars().size() > 0 && getCluster().at(0)->getBars().size()>0 );
    if(m_1dcluster->getSlayer() != getSlayer()  ) return false; 

    for(int i1d=0; i1d<m_1dclusters.size(); i1d++)
    {
      for(int ibar=0; ibar<m_1dcluster->getBars().size(); ibar++)
      {
        for(int jbar=0; jbar<m_1dclusters.at(i1d)->getBars().size(); jbar++)
        {
          if( m_1dcluster->getBars().at(ibar)->isLongiNeighbor(m_1dclusters.at(i1d)->getBars().at(jbar)) ) return true;
          if( m_1dcluster->getBars().at(ibar)->isLongiModuleAdjacent(m_1dclusters.at(i1d)->getBars().at(jbar)) ) return true;
        }
      }
    }
    
    return false;
  }

 
  void CaloHalfCluster::addUnit(const Calo1DCluster* _1dcluster)
  {
    if(_1dcluster->getSlayer()==0) slayer=0;
    if(_1dcluster->getSlayer()==1) slayer=1;
    m_1dclusters.push_back(_1dcluster);

    std::vector< std::vector<int> > id = _1dcluster->getTowerID();
    for(int ii=0; ii<id.size(); ii++)
      if( find(towerID.begin(), towerID.end(), id[ii])==towerID.end() ) towerID.push_back(id[ii]);    

    fitAxis("");
  }
  
  std::vector<const CaloUnit*> CaloHalfCluster::getBars() const
  {
    std::vector<const CaloUnit*> results;
    results.clear();
    for(int i=0; i<m_1dclusters.size(); i++)
    {
      for(int j=0; j<m_1dclusters.at(i)->getBars().size(); j++)
      {
        results.push_back(m_1dclusters.at(i)->getBars().at(j));
      }
    }
    return results;
  }


  double CaloHalfCluster::getEnergy() const {
    double sumE = 0;
    for(int i=0; i<m_1dclusters.size(); i++)
    {
      sumE = sumE + m_1dclusters.at(i)->getEnergy();
    }
    return sumE;
  }

  TVector3 CaloHalfCluster::getPos() const{
    TVector3 pos(0, 0, 0);
    double Etot = getEnergy();
    for(int i=0; i<m_1dclusters.size(); i++){
      TVector3 m_pos(m_1dclusters[i]->getPos().x(), m_1dclusters[i]->getPos().y(), m_1dclusters[i]->getPos().z());
      pos += m_pos * (m_1dclusters[i]->getEnergy()/Etot);
    }
    return pos;
  }

  TVector3 CaloHalfCluster::getEnergyCenter() const{
    TVector3 pos = getPos();
    double maxEn = -99;
    for(int i=0; i<m_1dclusters.size(); i++){
      if(m_1dclusters[i]->getEnergy()>maxEn){
        maxEn = m_1dclusters[i]->getEnergy();
        pos = m_1dclusters[i]->getPos();
      }
    }
    return pos;
  }

  std::vector<const PandoraPlus::Calo1DCluster*> CaloHalfCluster::getLocalMaxCol(std::string name) const{
    std::vector<const PandoraPlus::Calo1DCluster*> emptyCol; emptyCol.clear(); 
    if(map_localMax.find(name)!=map_localMax.end()) emptyCol = map_localMax.at(name);
    return emptyCol;
  }

  std::vector<const Calo1DCluster*> CaloHalfCluster::getAllLocalMaxCol() const{
    std::vector<const PandoraPlus::Calo1DCluster*> emptyCol; emptyCol.clear();
    for(auto &iter: map_localMax) emptyCol.insert(emptyCol.end(), iter.second.begin(), iter.second.end());
    return emptyCol;
  }

  std::vector<const PandoraPlus::CaloHalfCluster*> CaloHalfCluster::getHalfClusterCol(std::string name) const{
    std::vector<const PandoraPlus::CaloHalfCluster*> emptyCol; emptyCol.clear(); 
    if(map_halfClusCol.find(name)!=map_halfClusCol.end()) emptyCol = map_halfClusCol.at(name);
    return emptyCol;
  }

  std::vector<const PandoraPlus::CaloHalfCluster*> CaloHalfCluster::getAllHalfClusterCol() const{
    std::vector<const PandoraPlus::CaloHalfCluster*> emptyCol; emptyCol.clear();
    for(auto &iter: map_halfClusCol) emptyCol.insert(emptyCol.end(), iter.second.begin(), iter.second.end());
    return emptyCol;
  }

  std::vector<const PandoraPlus::Calo1DCluster*> CaloHalfCluster::getClusterInLayer(int _layer) const{
    std::vector<const PandoraPlus::Calo1DCluster*> outShowers; outShowers.clear();
    for(int i=0; i<m_1dclusters.size(); i++)
      if(m_1dclusters[i]->getDlayer()==_layer) outShowers.push_back(m_1dclusters[i]);
    return outShowers;
  }


  int CaloHalfCluster::getBeginningDlayer() const{
    int Lstart = 99;
    for(int i=0; i<m_1dclusters.size(); i++)
      if(m_1dclusters[i]->getDlayer()<Lstart) Lstart = m_1dclusters[i]->getDlayer();
    if(Lstart==99) return -99;
    else return Lstart;
  }


  int CaloHalfCluster::getEndDlayer() const{
    int Lend = -99;
    for(int i=0; i<m_1dclusters.size(); i++)
      if(m_1dclusters[i]->getDlayer()>Lend) Lend = m_1dclusters[i]->getDlayer();
    return Lend;
  }  

  bool CaloHalfCluster::isContinue() const{
    int ly_max = -1;
    int ly_min = 99;
    std::vector<int> vec_layers; vec_layers.clear();
    for(int il=0; il<m_1dclusters.size(); il++){
      vec_layers.push_back(m_1dclusters[il]->getDlayer());
      if(m_1dclusters[il]->getDlayer()>ly_max) ly_max = m_1dclusters[il]->getDlayer();
      if(m_1dclusters[il]->getDlayer()<ly_min) ly_min = m_1dclusters[il]->getDlayer();
    }

    bool flag = true;
    std::sort(vec_layers.begin(), vec_layers.end());
    for(int il=0; il<vec_layers.size() && vec_layers[il]!=ly_max; il++)
      if( find(vec_layers.begin(), vec_layers.end(), vec_layers[il]+1) == vec_layers.end() ){ flag=false; break; }

    return flag;
  }

  bool CaloHalfCluster::isContinueN(int n) const{
    if(n<=0) return true;

    int ly_max = -1;
    int ly_min = 99;
    std::vector<int> vec_layers; vec_layers.clear();
    for(int il=0; il<m_1dclusters.size(); il++){
      vec_layers.push_back(m_1dclusters[il]->getDlayer());
      if(m_1dclusters[il]->getDlayer()>ly_max) ly_max = m_1dclusters[il]->getDlayer();
      if(m_1dclusters[il]->getDlayer()<ly_min) ly_min = m_1dclusters[il]->getDlayer();
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


  bool CaloHalfCluster::isSubset(const CaloHalfCluster* clus) const{

    for(int is=0; is<clus->getCluster().size(); is++)
      if( find(m_1dclusters.begin(), m_1dclusters.end(), clus->getCluster()[is])==m_1dclusters.end() ) {return false; }

    if(m_1dclusters.size() > clus->getCluster().size())
      return true;
    else{
      if( TMath::Abs(Hough_rho) <= TMath::Abs(clus->getHoughRho()) ) { return true; }
      else { return false; }
    }
  }

  double CaloHalfCluster::OverlapRatioE( const CaloHalfCluster* clus) const{
    double Eshare = 0.;
    for(int is=0; is<clus->getCluster().size(); is++)
      if( find(m_1dclusters.begin(), m_1dclusters.end(), clus->getCluster()[is])!=m_1dclusters.end() ) Eshare += clus->getCluster()[is]->getEnergy();

    return Eshare / getEnergy() ;
  }


  void CaloHalfCluster::fitAxis( std::string name ){
    std::vector<const PandoraPlus::Calo1DCluster*> barShowerCol; barShowerCol.clear();
    if(!name.empty() && map_localMax.find(name)!=map_localMax.end() ) barShowerCol = map_localMax.at(name);
    else barShowerCol = m_1dclusters; 

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
      axis.SetPhi(fitPhi);
      axis.SetTheta(fitTheta);
      axis.SetMag(1.);
    }
  }

  void CaloHalfCluster::mergeHalfCluster(const CaloHalfCluster* clus ){
    for(int is=0; is<clus->getCluster().size(); is++){
      if( find(m_1dclusters.begin(), m_1dclusters.end(), clus->getCluster()[is])==m_1dclusters.end() )
        addUnit( clus->getCluster()[is] );
    }

    for(int itrk=0; itrk<clus->getAssociatedTracks().size(); itrk++){
      if( find(m_TrackCol.begin(), m_TrackCol.end(), clus->getAssociatedTracks()[itrk])==m_TrackCol.end() )
        m_TrackCol.push_back( clus->getAssociatedTracks()[itrk] );
    }

    for(auto iter:clus->getLocalMaxMap() ){
      if(map_localMax.find(iter.first)==map_localMax.end()) map_localMax[iter.first] = iter.second;
      else{
        for(int il=0; il<iter.second.size(); il++)
          if( find(map_localMax[iter.first].begin(), map_localMax[iter.first].end(), iter.second[il])==map_localMax[iter.first].end() ) 
            map_localMax[iter.first].push_back( iter.second[il] );
      }
    }

    for(auto iter:clus->getHalfClusterMap() ){
      if(map_halfClusCol.find(iter.first)==map_halfClusCol.end()) map_halfClusCol[iter.first] = iter.second;
      else{
        for(int il=0; il<iter.second.size(); il++)
          if( find(map_halfClusCol[iter.first].begin(), map_halfClusCol[iter.first].end(), iter.second[il])==map_halfClusCol[iter.first].end() )
            map_halfClusCol[iter.first].push_back( iter.second[il] );
      }
    }

    fitAxis("");
  }


};
#endif
