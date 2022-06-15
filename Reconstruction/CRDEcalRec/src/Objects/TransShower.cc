#ifndef TRANSSHOWER_C
#define TRANSSHOWER_C

#include "Objects/TransShower.h"
#include "TMath.h"

using namespace std;
namespace PandoraPlus{

  TVector3 TransShower::getPos() const{
    float rotAngle = module*TMath::Pi()/4.;
    TVector3 m_vecX = barShowerX->getPos();
    TVector3 m_vecY = barShowerY->getPos();
    m_vecX.RotateZ(rotAngle);
    m_vecY.RotateZ(rotAngle);
    TVector3 m_pos( m_vecY.x(), (m_vecX.y()+m_vecY.y())/2 , m_vecX.z() );
    m_pos.RotateZ(-rotAngle);
    return m_pos;
  }

/*
  double TransShower::getHitsE() const{
    double en=0;
    for(int i=0;i<CaloHits.size();i++) en+=CaloHits[i].getEnergy();
    return en;
  }

  double TransShower::getHitsWidth() const{
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

};
#endif
