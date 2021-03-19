#ifndef _CRD_ECAL_2DSHOWER_C
#define _CRD_ECAL_2DSHOWER_C
#define PI 3.141592653

#include "CRDEcalEDMSvc/CRDCaloHit2DShower.h"

using namespace std;
namespace CRDEcalEDM{

  CRDCaloHit2DShower::CRDCaloHit2DShower(CRDEcalEDM::CRDCaloBarShower _bars1, CRDEcalEDM::CRDCaloBarShower _bars2, std::vector<edm4hep::ConstCalorimeterHit> _hits ){
    barShowerX = _bars1;
    barShowerY = _bars2; 
    CaloHits = _hits;
  }

  dd4hep::Position CRDCaloHit2DShower::getPos() const{
    float rotAngle = module*PI/4.;
    TVector3 m_vecX(0,0,0);
    TVector3 m_vecY(0,0,0);
    m_vecX.SetXYZ(barShowerX.getPos().x(), barShowerX.getPos().y(),barShowerX.getPos().z());
    m_vecY.SetXYZ(barShowerY.getPos().x(), barShowerY.getPos().y(),barShowerY.getPos().z());
    m_vecX.RotateZ(rotAngle);
    m_vecY.RotateZ(rotAngle);
    TVector3 m_pos( m_vecY.x(), (m_vecX.y()+m_vecY.y())/2 , m_vecX.z() );
    m_pos.RotateZ(-rotAngle);
    dd4hep::Position pos(m_pos.x(), m_pos.y(), m_pos.z());
    return pos;
  }

  double CRDCaloHit2DShower::getHitsE() const{
    double en=0;
    for(int i=0;i<CaloHits.size();i++) en+=CaloHits[i].getEnergy();
    return en;
  }

};
#endif
