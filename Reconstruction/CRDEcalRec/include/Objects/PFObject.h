#ifndef _PFOBJECT_H
#define _PFOBJECT_H

#include "Objects/Track.h"
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
    CRDEcalEDM::CRDCaloHit3DCluster getECalShower() const { return m_EcalShower; }
    CRDEcalEDM::Track getTrack() const { return m_Track; }
    TLorentzVector getP4() const { return m_4vec; }
    bool ContainTrack() const { return f_hasTrk; }

    void Clear() { m_EcalShower.Clear(); m_Track.Clear(); m_particleID=-99; m_charge=-99; m_4vec.SetXYZT(0.,0.,0.,0.); f_hasTrk=false; }

    void setEcalShower( CRDEcalEDM::CRDCaloHit3DCluster& _sh) { m_EcalShower = _sh; }
    void setTrack( CRDEcalEDM::Track _trk) { m_Track = _trk; f_hasTrk = true;}
    void setPdgID( int _id ) { m_particleID = _id; }
    void setCharge( int _ch) { m_charge = _ch; }
    void setP4( TLorentzVector& _vec ) { m_4vec = _vec; }

  private:
    CRDEcalEDM::CRDCaloHit3DCluster   m_EcalShower; 
    CRDEcalEDM::Track                m_Track;

    int    m_particleID; 
    int    m_charge;
    TLorentzVector m_4vec;

    bool f_hasTrk; 
  };
};
#endif 
