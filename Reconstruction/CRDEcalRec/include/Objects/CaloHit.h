#ifndef CALO_HIT_H
#define CALO_HIT_H

#include "TVector3.h"
#include "Objects/TransShower.h"

namespace PandoraPlus{
  class TransShower; 

  class CaloHit{
  public:
    CaloHit() {};
    ~CaloHit() { Clear(); };

    void Clear() { cellID=0; position.SetXYZ(0.,0.,0.); energy=-1; layer=-1; ParentShower=nullptr; }

    TVector3 getPosition() const { return position; }
    double   getEnergy() const { return energy; } 
    int getLayer() const {return layer;}

    void setcellID(unsigned long long _id) { cellID=_id; }
    void setEnergy(double _en) { energy=_en; }
    void setPosition( TVector3 _vec ) { position=_vec; }
    void setLayer(int _l) { layer = _l; }
    void setParentShower( PandoraPlus::TransShower* _p ) { ParentShower=_p; }

  private: 
    int layer; 
    unsigned long long cellID; 
    TVector3 position;
    double   energy; 
    PandoraPlus::TransShower* ParentShower; 
  };

};
#endif
