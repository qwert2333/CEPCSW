#ifndef _CRD_ECAL_TRANSSHOWER_C
#define _CRD_ECAL_TRANSSHOWER_C
#define PI 3.141592653

#include "Objects/CRDCaloHitTransShower.h"

using namespace std;
namespace CRDEcalEDM{

  CRDCaloHitTransShower::CRDCaloHitTransShower(CRDEcalEDM::CRDCaloBarShower _bars1, CRDEcalEDM::CRDCaloBarShower _bars2 /*, std::vector<edm4hep::ConstCalorimeterHit> _hits*/ ){
    barShowerX = _bars1;
    barShowerY = _bars2; 
    //CaloHits = _hits;
  }

  TVector3 CRDCaloHitTransShower::getPos() const{
    float rotAngle = module*PI/4.;
    TVector3 m_vecX(0,0,0);
    TVector3 m_vecY(0,0,0);
    m_vecX.SetXYZ(barShowerX.getPos().x(), barShowerX.getPos().y(),barShowerX.getPos().z());
    m_vecY.SetXYZ(barShowerY.getPos().x(), barShowerY.getPos().y(),barShowerY.getPos().z());
    m_vecX.RotateZ(rotAngle);
    m_vecY.RotateZ(rotAngle);
    TVector3 m_pos( m_vecY.x(), (m_vecX.y()+m_vecY.y())/2 , m_vecX.z() );
    m_pos.RotateZ(-rotAngle);
    TVector3 pos(m_pos.x(), m_pos.y(), m_pos.z());
    return pos;
  }

/*
  double CRDCaloHitTransShower::getHitsE() const{
    double en=0;
    for(int i=0;i<CaloHits.size();i++) en+=CaloHits[i].getEnergy();
    return en;
  }

  double CRDCaloHitTransShower::getHitsWidth() const{
    TVector3 cent(0., 0., 0.); 
    double totE = getHitsE(); 
    for(int i=0;i<CaloHits.size();i++){
      TVector3 pos(CaloHits[i].getPosition().x, CaloHits[i].getPosition().y, CaloHits[i].getPosition().z);
      cent += pos* (CaloHits[i].getEnergy()/totE); 
    }

    double width=0; 
    for(int i=0;i<CaloHits.size();i++){
      TVector3 pos(CaloHits[i].getPosition().x, CaloHits[i].getPosition().y, CaloHits[i].getPosition().z);
      double r2 = (pos-cent).Mag2(); 
      width += r2*CaloHits[i].getEnergy()/totE;
    }

    return width; 
  }
*/
  void CRDCaloHitTransShower::setShadowClusType( int _flag ){
    if(_flag==1)      { m_isTrkShower=true; m_isNeuShower=false; }
    else if(_flag==0) { m_isTrkShower=false; m_isNeuShower=true; }
    else              { m_isTrkShower=false; m_isNeuShower=false;}
  }

};
#endif
