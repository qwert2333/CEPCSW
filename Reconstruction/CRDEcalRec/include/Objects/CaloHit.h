#ifndef CALO_HIT_H
#define CALO_HIT_H

#include "TVector3.h"

namespace PandoraPlus{

  class CaloHit{
  public:
    CaloHit() {};
    ~CaloHit() {};

    void Clear() { cellID=0; position.SetXYZ(0.,0.,0.); energy=-1; }

    TVector3 getPosition() const { return position; }
    double   getEnergy() const { return energy; } 

    void setcellID(unsigned long long& _id) { cellID=_id; }
    void setEnergy(double& _en) { energy=_en; }
    void setPosition( TVector3& _vec ) { position=_vec; }

  private: 
    unsigned long long cellID; 
    TVector3 position;
    double   energy; 
  };

};
#endif
