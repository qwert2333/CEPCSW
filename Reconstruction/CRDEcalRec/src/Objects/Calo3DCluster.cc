#ifndef CALO_3DCLUSTER_C
#define CALO_3DCLUSTER_C

#include "Objects/Calo3DCluster.h"
#include <cmath>
using namespace std;
namespace PandoraPlus{

  void Calo3DCluster::Clear() 
  {
	  m_2dclusters.clear();
	  m_towers.clear();
    towerID.clear(); 
	  //m_modules.clear();
	  //m_parts.clear();
	  //m_staves.clear();
  }

  void Calo3DCluster::Clean(){
    //for(int i=0; i<m_2dclusters.size(); i++) { delete m_2dclusters[i]; m_2dclusters[i]=NULL; }
    //for(int i=0; i<m_towers.size(); i++) { delete m_towers[i]; m_towers[i]=NULL; }
    //std::vector<int>().swap(m_modules);
    //std::vector<int>().swap(m_parts);
    //std::vector<int>().swap(m_staves);
    Clear();
  }

  void Calo3DCluster::Check()
  {
    for(int i=0; i<m_2dclusters.size(); i++)
    if(!m_2dclusters[i]) { m_2dclusters.erase(m_2dclusters.begin()+i); i--; }
    //for(int i=0; i<m_towers.size(); i++)
    //if(!m_towers[i]) { m_towers.erase(m_towers.begin()+i); i--; }
  }

  //TODO: This function is sooooooo time consumeing now!!
  bool Calo3DCluster::isNeighbor(const PandoraPlus::Calo2DCluster* m_2dcluster) const
  {

    //Inner module
    for(int i=0; i<m_2dcluster->getTowerID().size(); i++){
    for(int j=0; j<towerID.size(); j++){
      if( m_2dcluster->getTowerID()[i]==towerID[j] ) return true;
    }}

    //for(int i=0; i<m_2dcluster->getModules().size(); i++){
    //  for(int j=0; j<m_modules.size(); j++){
    //    if(m_2dcluster->getModules().at(i)==m_modules.at(j) && m_2dcluster->getParts().at(i)==m_parts.at(j) && m_2dcluster->getStaves().at(i)==m_staves.at(j)){
    //      return true;
    //  }}
    //}

    //Inner module but in adjacent Dlayer
    for(int i=0; i<m_2dclusters.size(); i++){
      //If not in the same module:
      bool inSameModule = false; 
      std::vector<int> ID_module; ID_module.clear();
      for(int it=0; it<towerID.size(); it++) ID_module.push_back(towerID[it][0]); 
      for(int it=0; it<m_2dcluster->getTowerID().size(); it++){
        if( find(ID_module.begin(), ID_module.end(), m_2dcluster->getTowerID()[it][0])!=ID_module.end() ) { inSameModule=true; break; }
      }
      if(!inSameModule) continue;

      //If not the adjacent Dlayer: 
      if(fabs(m_2dcluster->getDlayer()-m_2dclusters[i]->getDlayer())>1) continue; 

      std::vector<const PandoraPlus::CaloUnit*> m_EdgeBarU_clus1; m_EdgeBarU_clus1.clear(); //Income 2DCluster
      std::vector<const PandoraPlus::CaloUnit*> m_EdgeBarV_clus1; m_EdgeBarV_clus1.clear();
      std::vector<const PandoraPlus::CaloUnit*> m_EdgeBarU_clus2; m_EdgeBarU_clus2.clear(); //2DCluster in present 3DCluster. 
      std::vector<const PandoraPlus::CaloUnit*> m_EdgeBarV_clus2; m_EdgeBarV_clus2.clear();
      
      for(int ib=0; ib<m_2dcluster->getBarUCol().size(); ib++){
        if(m_2dcluster->getBarUCol()[ib]->isAtLowerEdgeZ() || m_2dcluster->getBarUCol()[ib]->isAtUpperEdgeZ()) m_EdgeBarU_clus1.push_back(m_2dcluster->getBarUCol()[ib]);
      }
      for(int ib=0; ib<m_2dcluster->getBarVCol().size(); ib++){
        if(m_2dcluster->getBarVCol()[ib]->isAtLowerEdgePhi() || m_2dcluster->getBarVCol()[ib]->isAtUpperEdgePhi()) m_EdgeBarV_clus1.push_back(m_2dcluster->getBarVCol()[ib]);
      }
      for(int ib=0; ib<m_2dclusters[i]->getBarUCol().size(); ib++){
        if(m_2dclusters[i]->getBarUCol()[ib]->isAtLowerEdgeZ() || m_2dclusters[i]->getBarUCol()[ib]->isAtUpperEdgeZ()) m_EdgeBarU_clus2.push_back(m_2dclusters[i]->getBarUCol()[ib]);
      }
      for(int ib=0; ib<m_2dclusters[i]->getBarVCol().size(); ib++){
        if(m_2dclusters[i]->getBarVCol()[ib]->isAtLowerEdgePhi() || m_2dclusters[i]->getBarVCol()[ib]->isAtUpperEdgePhi()) m_EdgeBarV_clus2.push_back(m_2dclusters[i]->getBarVCol()[ib]);
      }


      for(int ib1=0; ib1<m_EdgeBarU_clus1.size(); ib1++){
        for(int ib2=0; ib2<m_2dclusters[i]->getBarVCol().size(); ib2++){
          if( m_EdgeBarU_clus1[ib1]->getPart()==m_2dclusters[i]->getBarVCol()[ib2]->getPart() && abs(m_EdgeBarU_clus1[ib1]->getStave()-m_2dclusters[i]->getBarVCol()[ib2]->getStave())==1 ) 
            return true; 
      }}
      for(int ib1=0; ib1<m_EdgeBarV_clus1.size(); ib1++){
        for(int ib2=0; ib2<m_2dclusters[i]->getBarUCol().size(); ib2++){
          if( m_EdgeBarV_clus1[ib1]->getStave()==m_2dclusters[i]->getBarUCol()[ib2]->getStave() && abs(m_EdgeBarV_clus1[ib1]->getPart()-m_2dclusters[i]->getBarUCol()[ib2]->getPart())==1 )
            return true;
      }}
      for(int ib1=0; ib1<m_EdgeBarU_clus2.size(); ib1++){
        for(int ib2=0; ib2<m_2dcluster->getBarVCol().size(); ib2++){
          if( m_EdgeBarU_clus2[ib1]->getPart()==m_2dcluster->getBarVCol()[ib2]->getPart() && abs(m_EdgeBarU_clus2[ib1]->getStave()-m_2dcluster->getBarVCol()[ib2]->getStave())==1 )
            return true;
      }}
      for(int ib1=0; ib1<m_EdgeBarV_clus2.size(); ib1++){
        for(int ib2=0; ib2<m_2dcluster->getBarUCol().size(); ib2++){
          if( m_EdgeBarV_clus2[ib1]->getStave()==m_2dcluster->getBarUCol()[ib2]->getStave() && abs(m_EdgeBarV_clus2[ib1]->getPart()-m_2dcluster->getBarUCol()[ib2]->getPart())==1 )
            return true;
      }}

	  }

    //Over modules
    std::vector<const PandoraPlus::CaloUnit*> bars_2d = m_2dcluster->getBars();
    for(int ib2d=0; ib2d<bars_2d.size(); ib2d++){
      for(int ic=0; ic<m_2dclusters.size(); ic++){
        for(int ib3d=0; ib3d<m_2dclusters[ic]->getBars().size(); ib3d++){
          if(bars_2d[ib2d]->isModuleAdjacent(m_2dclusters[ic]->getBars()[ib3d])) return true;
        }
      }
    }
	  return false;
  }

  void Calo3DCluster::addUnit(const Calo2DCluster* _2dcluster){

    m_2dclusters.push_back(_2dcluster);
    std::vector< std::vector<int> > id = _2dcluster->getTowerID();
    for(int ii=0; ii<id.size(); ii++)
      if( find(towerID.begin(), towerID.end(), id[ii])==towerID.end() ) towerID.push_back(id[ii]);
  }


  void Calo3DCluster::mergeCluster( const PandoraPlus::Calo3DCluster* _clus ){
    for(int i=0; i<_clus->getCluster().size(); i++)
      addUnit( _clus->getCluster()[i] );

    for(int i=0; i<_clus->getTowers().size(); i++)
      addTower( _clus->getTowers()[i] );

    for(int itrk=0; itrk<_clus->getAssociatedTracks().size(); itrk++){
      if( find(m_TrackCol.begin(), m_TrackCol.end(), _clus->getAssociatedTracks()[itrk])==m_TrackCol.end() ) 
        m_TrackCol.push_back( _clus->getAssociatedTracks()[itrk] );
    }
    
    for(auto iter:_clus->getLocalMaxUMap() ){
      if(map_localMaxU.find(iter.first)==map_localMaxU.end()) map_localMaxU[iter.first] = iter.second;
      else{
        for(int il=0; il<iter.second.size(); il++)
          if( find(map_localMaxU[iter.first].begin(), map_localMaxU[iter.first].end(), iter.second[il])==map_localMaxU[iter.first].end() )
            map_localMaxU[iter.first].push_back( iter.second[il] );
      }
    }
    for(auto iter:_clus->getLocalMaxVMap() ){
      if(map_localMaxV.find(iter.first)==map_localMaxV.end()) map_localMaxV[iter.first] = iter.second;
      else{
        for(int il=0; il<iter.second.size(); il++)
          if( find(map_localMaxV[iter.first].begin(), map_localMaxV[iter.first].end(), iter.second[il])==map_localMaxV[iter.first].end() )
            map_localMaxV[iter.first].push_back( iter.second[il] );
      }
    }

    for(auto iter:_clus->getHalfClusterUMap() ){
      if(map_halfClusUCol.find(iter.first)==map_halfClusUCol.end()) map_halfClusUCol[iter.first] = iter.second;
      else{
        for(int il=0; il<iter.second.size(); il++)
          if( find(map_halfClusUCol[iter.first].begin(), map_halfClusUCol[iter.first].end(), iter.second[il])==map_halfClusUCol[iter.first].end() )
            map_halfClusUCol[iter.first].push_back( iter.second[il] );
      }
    }
    for(auto iter:_clus->getHalfClusterVMap() ){
      if(map_halfClusVCol.find(iter.first)==map_halfClusVCol.end()) map_halfClusVCol[iter.first] = iter.second;
      else{
        for(int il=0; il<iter.second.size(); il++)
          if( find(map_halfClusVCol[iter.first].begin(), map_halfClusVCol[iter.first].end(), iter.second[il])==map_halfClusVCol[iter.first].end() )
            map_halfClusVCol[iter.first].push_back( iter.second[il] );
      }
    }

  }


  std::vector<const PandoraPlus::CaloUnit*> Calo3DCluster::getBars() const{
    std::vector<const PandoraPlus::CaloUnit*> results; results.clear();
    for(int i=0; i<m_2dclusters.size(); i++){
      for(int j=0; j<m_2dclusters.at(i)->getBars().size(); j++){
        results.push_back(m_2dclusters.at(i)->getBars().at(j));
      }
    }
    return results;
  }

  double Calo3DCluster::getEnergy() const{
    double result = 0;
    for(int m=0; m<m_2dclusters.size(); m++)
      result += m_2dclusters[m]->getEnergy();
    return result;
  }

  double Calo3DCluster::getHitsE() const{
    double en=0;
    for(int i=0;i<hits.size(); i++) en+=hits[i]->getEnergy();
    return en;
  }

  double Calo3DCluster::getLongiE() const{
    double en=0;
    for(auto iter: map_halfClusUCol){
      if(iter.first!="LinkedLongiCluster") continue;
      for(auto iclus: iter.second)
        en += iclus->getEnergy();
    }
    for(auto iter: map_halfClusVCol){
      if(iter.first!="LinkedLongiCluster") continue;
      for(auto iclus: iter.second)
        en += iclus->getEnergy();
    }
    return en;
  }

  TVector3 Calo3DCluster::getHitCenter() const{
    TVector3 vec(0,0,0);
    double totE = getHitsE();
    for(int i=0;i<hits.size(); i++){
       TVector3 v_cent = hits[i]->getPosition();
       vec += v_cent * (hits[i]->getEnergy()/totE);
    }
    return vec;
  }

  TVector3 Calo3DCluster::getShowerCenter() const{
    TVector3 spos(0,0,0);
    double totE = 0.;
    for(int i=0;i<m_2dclusters.size(); i++){ spos += m_2dclusters[i]->getPos()*m_2dclusters[i]->getEnergy(); totE += m_2dclusters[i]->getEnergy(); }
    spos = spos*(1./totE);
    return spos;
  }

  int Calo3DCluster::getBeginningDlayer() const{
    int re_dlayer = -99;
    std::vector<int> dlayers; dlayers.clear();
    if(hits.size()!=0)  for(int ih=0; ih<hits.size(); ih++) dlayers.push_back(hits[ih]->getLayer());
    else                for(int ish=0; ish<m_2dclusters.size(); ish++) dlayers.push_back(m_2dclusters[ish]->getDlayer());
    re_dlayer = *std::min_element(dlayers.begin(), dlayers.end());

    return re_dlayer;
  }

  int Calo3DCluster::getEndDlayer() const{
    int re_dlayer = -99;
    std::vector<int> dlayers; dlayers.clear();
    if(hits.size()!=0)  for(int ih=0; ih<hits.size(); ih++) dlayers.push_back(hits[ih]->getLayer());
    else                for(int ish=0; ish<m_2dclusters.size(); ish++) dlayers.push_back(m_2dclusters[ish]->getDlayer());
    re_dlayer = *std::max_element(dlayers.begin(), dlayers.end());

    return re_dlayer;
  }


  std::vector<const PandoraPlus::CaloHalfCluster*> Calo3DCluster::getHalfClusterUCol(std::string name) const {
    std::vector<const CaloHalfCluster*> emptyCol; emptyCol.clear(); 
    if(map_halfClusUCol.find(name)!=map_halfClusUCol.end()) emptyCol = map_halfClusUCol.at(name);
    return emptyCol;
  }

  std::vector<const PandoraPlus::CaloHalfCluster*> Calo3DCluster::getHalfClusterVCol(std::string name) const {
    std::vector<const CaloHalfCluster*> emptyCol; emptyCol.clear(); 
    if(map_halfClusVCol.find(name)!=map_halfClusVCol.end()) emptyCol = map_halfClusVCol.at(name);
    return emptyCol;
  }

  std::vector<const PandoraPlus::Calo1DCluster*> Calo3DCluster::getLocalMaxUCol(std::string name) const{
    std::vector<const Calo1DCluster*> emptyCol; emptyCol.clear(); 
    if(map_localMaxU.find(name)!=map_localMaxU.end()) emptyCol = map_localMaxU.at(name);
    return emptyCol;
  }

  std::vector<const PandoraPlus::Calo1DCluster*> Calo3DCluster::getLocalMaxVCol(std::string name) const{
    std::vector<const Calo1DCluster*> emptyCol; emptyCol.clear(); 
    if(map_localMaxV.find(name)!=map_localMaxV.end()) emptyCol = map_localMaxV.at(name);
    return emptyCol; 
  }

  std::vector< std::pair<edm4hep::MCParticle, float> > Calo3DCluster::getLinkedMCPfromHFCluster(std::string name){
    MCParticleWeight.clear();

    std::map<edm4hep::MCParticle, float> map_truthP_totE; map_truthP_totE.clear();
    for(auto icl: map_halfClusUCol[name]){
      for(auto ipair: icl->getLinkedMCP()) map_truthP_totE[ipair.first] += icl->getEnergy()*ipair.second;
    }
    for(auto icl: map_halfClusVCol[name]){
      for(auto ipair: icl->getLinkedMCP()) map_truthP_totE[ipair.first] += icl->getEnergy()*ipair.second;
    }

    for(auto imcp: map_truthP_totE){
      MCParticleWeight.push_back( std::make_pair(imcp.first, imcp.second/getLongiE()) );
    }

    return MCParticleWeight;
  }


  void Calo3DCluster::FitAxis(){
    if(m_2dclusters.size()==0) axis.SetXYZ(0,0,0);

    else if(m_2dclusters.size()==1){
      axis = m_2dclusters[0]->getPos();
      axis *= 1./axis.Mag();
    }

    else if( m_2dclusters.size()==2 ){
      TVector3 pos1 = m_2dclusters[0]->getPos();
      TVector3 pos2 = m_2dclusters[1]->getPos();

      axis = ( pos1.Mag()>pos2.Mag() ? pos1-pos2 : pos2-pos1 );
      axis *= 1./axis.Mag();
    }

    else{
      trackFitter.clear();
      //track->setImpactParameter(0., 0.); //fix dr=0, dz=0.

      double barAngle = (towerID[0][0]+2)*TMath::Pi()/4.;
      double posErr = 10./sqrt(12);
      if(barAngle>=TMath::TwoPi()) barAngle = barAngle-TMath::TwoPi();
      trackFitter.setBarAngle(barAngle);
      for(int is=0;is<m_2dclusters.size();is++){
        TVector3 pos_barsX = m_2dclusters[is]->getShowerUCol()[0]->getPos(); //U
        TVector3 pos_barsY = m_2dclusters[is]->getShowerVCol()[0]->getPos(); //Z
//printf("\t DEBUG: input pointX (%.3f, %.3f, %.3f) \n", barsX.getPos().x(), barsX.getPos().y(), barsX.getPos().z());
//printf("\t DEBUG: input pointY (%.3f, %.3f, %.3f) \n", barsY.getPos().x(), barsY.getPos().y(), barsY.getPos().z());
        trackFitter.setGlobalPoint(1, pos_barsX.x(), posErr, pos_barsX.y(), posErr, pos_barsX.z(), posErr);
        trackFitter.setGlobalPoint(0, pos_barsY.x(), posErr, pos_barsY.y(), posErr, pos_barsY.z(), posErr);
      }
      trackFitter.fitTrack();
      double fitPhi =   trackFitter.getTrkPar(2);
      double fitTheta = trackFitter.getTrkPar(3);
//printf("\t DEBUG: fitted phi and theta: %.3f \t %.3f \n", fitPhi, fitTheta);

      axis.SetPhi(fitPhi);
      axis.SetTheta(fitTheta);
      axis.SetMag(1.);
    }
  }

  //void Calo3DCluster::FitAxisHit(){
  //}

  //void Calo3DCluster::FitProfile(){
  //}


};
#endif
