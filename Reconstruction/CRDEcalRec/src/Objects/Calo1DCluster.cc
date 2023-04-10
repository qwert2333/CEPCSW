#ifndef CALO_1DCLUSTER_C
#define CALO_1DCLUSTER_C

#include "Objects/CaloUnit.h"
#include "Objects/Calo1DCluster.h"
using namespace std; 
namespace PandoraPlus{

  void Calo1DCluster::Clear(){
    Bars.clear();
    Seeds.clear();
    Energy=0.;
    pos.SetXYZ(0.,0.,0.);
    towerID.clear();
    //m_modules.clear();
    //m_parts.clear();
    //m_staves.clear();
  }

  void Calo1DCluster::Clean() {
    for(int i=0; i<Bars.size(); i++) { delete Bars[i]; Bars[i]=NULL; }
    for(int i=0; i<Seeds.size(); i++) { delete Seeds[i]; Seeds[i]=NULL; }
    //std::vector<int>().swap(m_modules);
    //std::vector<int>().swap(m_parts);
    //std::vector<int>().swap(m_staves);
    Clear();
  }

  void Calo1DCluster::Check() {
    for(int i=0; i<Bars.size(); i++)
      if(!Bars[i]) { Bars.erase(Bars.begin()+i); i--; }
    for(int i=0; i<Seeds.size(); i++)
      if(!Seeds[i]) { Seeds.erase(Seeds.begin()+i); i--; }
  }

  PandoraPlus::Calo1DCluster* Calo1DCluster::Clone() const{
    PandoraPlus::Calo1DCluster* p_cluster = new PandoraPlus::Calo1DCluster();
    p_cluster->setBars( Bars );
    p_cluster->setSeeds( Seeds );
    p_cluster->setIDInfo();
    for(int i=0; i<CousinClusters.size(); i++) p_cluster->addCousinCluster( CousinClusters[i] );
    for(int i=0; i<ChildClusters.size(); i++) p_cluster->addChildCluster( ChildClusters[i] );
   
    return p_cluster;
  }


  bool Calo1DCluster::isNeighbor(const PandoraPlus::CaloUnit* m_bar) const
  {
    for(int i1d = 0; i1d<Bars.size(); i1d++){
      if(Bars[i1d]->isNeighbor(m_bar)){ /*cout<<" isNeighbor! "<<endl;*/ return true;}
    }
    return false;
  }

  bool Calo1DCluster::inCluster(const PandoraPlus::CaloUnit* iBar) const{
    return (find(Bars.begin(), Bars.end(), iBar)!=Bars.end() );
  }

  double Calo1DCluster::getEnergy() const{
    double E=0;
    for(int i=0;i<Bars.size();i++) E+=Bars[i]->getEnergy();
    return E;
  }

  TVector3 Calo1DCluster::getPos() const{
    TVector3 pos(0,0,0);
    double Etot=getEnergy();
    for(int i=0;i<Bars.size();i++) pos += Bars[i]->getPosition() * (Bars[i]->getEnergy()/Etot);
    return pos;
  }

  double Calo1DCluster::getT1() const{
    double T1=0;
    double Etot = getEnergy();
    for(int i=0;i<Bars.size();i++) T1 += (Bars[i]->getT1() * Bars[i]->getEnergy())/Etot;
    return T1;
  }

  double Calo1DCluster::getT2() const{
    double T2=0;
    double Etot = getEnergy();
    for(int i=0;i<Bars.size();i++) T2 += (Bars[i]->getT2() * Bars[i]->getEnergy())/Etot;
    return T2;
  }

  double Calo1DCluster::getWidth() const{
    TVector3 centPos = getPos();
    double Etot = getEnergy();
    double sigmax=0;
    double sigmay=0;
    double sigmaz=0;
    for(int i=0; i<Bars.size(); i++){
      double wi = Bars[i]->getEnergy()/Etot;
      sigmax += wi*(Bars[i]->getPosition().x()-centPos.x())*(Bars[i]->getPosition().x()-centPos.x());
      sigmay += wi*(Bars[i]->getPosition().y()-centPos.y())*(Bars[i]->getPosition().y()-centPos.y());
      sigmaz += wi*(Bars[i]->getPosition().z()-centPos.z())*(Bars[i]->getPosition().z()-centPos.z());
    }

    if(sigmaz!=0) return sigmaz;  //sLayer=1, bars along z-axis.
    else if(sigmax==0 && sigmaz==0) return sigmay; //Module 2, 6
    else if(sigmay==0 && sigmaz==0) return sigmax; //Module 0, 4
    else if(sigmax!=0 && sigmay!=0 && sigmaz==0) return sqrt(sigmax*sigmax+sigmay*sigmay); //Module 1, 3, 5, 7;
    else return 0.;
  }

  double Calo1DCluster::getScndMoment() const{
    if(Bars.size()<=1) return 0.;
    TVector3 pos = getPos();
    double Etot = getEnergy();
    double scndM = 0;
    for(int i=0;i<Bars.size();i++) scndM += (Bars[i]->getEnergy() * (pos-Bars[i]->getPosition()).Mag2()) / Etot;
    return scndM;
  }

  bool Calo1DCluster::getGlobalRange( double& xmin,  double& ymin, double& zmin, double& xmax, double& ymax, double& zmax ) const{
    if(Bars.size()==0) return false; 

    std::vector<double> posx; posx.clear();
    std::vector<double> posy; posx.clear();
    std::vector<double> posz; posx.clear();
    for(int i=0; i<Bars.size(); i++){  
      posx.push_back(Bars[i]->getPosition().x());  
      posy.push_back(Bars[i]->getPosition().y());
      posz.push_back(Bars[i]->getPosition().z());
    }
    auto xminptr = std::min_element( std::begin(posx), std::end(posx) );
    auto xmaxptr = std::max_element( std::begin(posx), std::end(posx) );
    auto yminptr = std::min_element( std::begin(posy), std::end(posy) );
    auto ymaxptr = std::max_element( std::begin(posy), std::end(posy) );
    auto zminptr = std::min_element( std::begin(posz), std::end(posz) );
    auto zmaxptr = std::max_element( std::begin(posz), std::end(posz) );

    int slayer = Bars[0]->getSlayer(); 
    if(slayer==0){
      xmin = -9999.;    xmax = 9999.; 
      ymin = -9999.;    ymax = 9999.; 
      zmin = *zminptr-6.;  zmax = *zmaxptr+6.;
    }
    else if(slayer==1){
      xmin = *xminptr-6;  xmax = *xmaxptr+6;
      ymin = *yminptr-6;  ymax = *ymaxptr+6;
      zmin = -9999;     zmax = 9999; 
    }
    else return false; 

    return true;
  }

  int Calo1DCluster::getLeftEdge(){
    std::sort(Bars.begin(), Bars.end());
    if(Bars.size()==0) return -99;
    int edge = -99;
    if( Bars[0]->getSlayer()==0 ) edge = Bars[0]->getBar() + Bars[0]->getStave()*PandoraPlus::CaloUnit::NbarZ; 
    if( Bars[0]->getSlayer()==1 ) edge = Bars[0]->getBar() + Bars[0]->getPart()*(CaloUnit::NbarPhi-Bars[0]->getDlayer()+1);
    return edge;
  }

  int Calo1DCluster::getRightEdge(){
    std::sort(Bars.begin(), Bars.end());
    if(Bars.size()==0) return -99;
    int edge = -99;
    if( Bars[Bars.size()-1]->getSlayer()==0 ) edge = Bars[Bars.size()-1]->getBar() + Bars[Bars.size()-1]->getStave()*CaloUnit::NbarZ;
    if( Bars[Bars.size()-1]->getSlayer()==1 ) edge = Bars[Bars.size()-1]->getBar() + Bars[Bars.size()-1]->getPart()*(CaloUnit::NbarPhi-Bars[Bars.size()-1]->getDlayer()+1);
    return edge;
  }

  void Calo1DCluster::addUnit(const PandoraPlus::CaloUnit* _bar ) 
  {
    Bars.push_back(_bar);
    std::vector<int> id(3);
    id[0] = _bar->getModule();
    id[1] = _bar->getPart();
    id[2] = _bar->getStave();
    if(find(towerID.begin(), towerID.end(), id)==towerID.end()) towerID.push_back(id);

    //m_modules.push_back(_bar->getModule());
    //m_parts.push_back(_bar->getPart());
    //m_staves.push_back(_bar->getStave());
    //if( find( m_modules.begin(), m_modules.end(), _bar->getModule())==m_modules.end() ) m_modules.push_back(_bar->getModule());
    //if( find( m_parts.begin(), m_parts.end(), _bar->getPart())==m_parts.end() )         m_parts.push_back(_bar->getPart());
    //if( find( m_staves.begin(), m_staves.end(), _bar->getStave())==m_staves.end() )     m_staves.push_back(_bar->getStave());
  }


  void Calo1DCluster::setIDInfo() {
    for(int i=0; i<Bars.size(); i++){
      std::vector<int> id(3);
      id[0] = Bars[i]->getModule();
      id[1] = Bars[i]->getPart();
      id[2] = Bars[i]->getStave();
      if(find(towerID.begin(), towerID.end(), id)==towerID.end()) towerID.push_back(id);  
    }
  }

  void Calo1DCluster::setSeed(){
    double maxE = -1; 
    int index = -1; 
    for(int i=0; i<Bars.size(); ++i){
      if(Bars[i]->getEnergy()>maxE) { maxE = Bars[i]->getEnergy(); index = i; }
    }
    //if(Seeds.size()==0 && index>=0 && maxE>0.005 ) Seeds.push_back( Bars[index] );
    if(Seeds.size()==0 && index>=0 ) Seeds.push_back( Bars[index] );
  }

};
#endif
