#ifndef CALO_HIT_H
#define CALO_HIT_H

#include "TVector3.h"
#include "edm4hep/MCParticle.h"
#include "Objects/Calo2DCluster.h"

namespace PandoraPlus{
  class Calo2DCluster; 

  class CaloHit{
  public:
    CaloHit() {};
    ~CaloHit() { Clear(); };

    void Clear() { cellID=0; position.SetXYZ(0.,0.,0.); energy=-1; layer=-1; ParentShower=nullptr; }

    TVector3 getPosition() const { return position; }
    double   getEnergy() const { return energy; } 
    int getLayer() const {return layer;}
    std::vector< std::pair<edm4hep::MCParticle, float> > getLinkedMCP() const { return MCParticleWeight; }
    edm4hep::MCParticle getLeadingMCP() const;
    float getLeadingMCPweight() const;

    void setcellID(unsigned long long _id) { cellID=_id; }
    void setEnergy(double _en) { energy=_en; }
    void setPosition( TVector3 _vec ) { position=_vec; }
    void setLayer(int _l) { layer = _l; }
    void setParentShower( PandoraPlus::Calo2DCluster* _p ) { ParentShower=_p; }
    void addLinkedMCP( std::pair<edm4hep::MCParticle, float> _pair ) { MCParticleWeight.push_back(_pair); }
    void setLinkedMCP( std::vector<std::pair<edm4hep::MCParticle, float>> _pairVec ) { MCParticleWeight.clear(); MCParticleWeight = _pairVec; }

  private: 
    int layer; 
    unsigned long long cellID; 
    TVector3 position;
    double   energy; 
    PandoraPlus::Calo2DCluster* ParentShower; 
    std::vector< std::pair<edm4hep::MCParticle, float> > MCParticleWeight;
  };

};
#endif
