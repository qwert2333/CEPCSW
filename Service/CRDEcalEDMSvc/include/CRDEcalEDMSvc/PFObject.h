#ifndef _PFOBJECT_H
#define _PFOBJECT_H

#include "TLorentzVector.h"
#include "TMath.h"
namespace CRDEcalEDM {

  class PFObject{
  public:
    PFObject(){};
    ~PFObject(){};

    int    getPdgID()  const { return m_particleID; }
    int    getCharge() const { return m_charge; }
    double getMass()   const { return m_4vec.M(); }
    double getEnergy() const { return m_4vec.E(); }
    CRDEcalEDM::CRDCaloHit3DShower getECalShower() const { return m_EcalShower; }
    TLorentzVector getP4() const { return m_4vec; }

    void Clear() { m_EcalShower.Clear(); m_particleID=-99; m_charge=-99; m_4vec.SetXYZT(0.,0.,0.,0.); }

    void setEcalShower( CRDEcalEDM::CRDCaloHit3DShower& _sh) { m_EcalShower = _sh; }
    void setPdgID( int _id ) { m_particleID = _id; }
    void setCharge( int _ch) { m_charge = _ch; }
    void setP4( TLorentzVector& _vec ) { m_4vec = _vec; }

  private:
    CRDEcalEDM::CRDCaloHit3DShower   m_EcalShower; 

    int    m_particleID; 
    int    m_charge;
    TLorentzVector m_4vec;

  };
};
#endif 
