#ifndef CALO_HIT_H
#define CALO_HIT_H

#include "TVector3.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/CalorimeterHit.h"
#include "Objects/Calo2DCluster.h"

namespace Cyber{
  class Calo2DCluster; 

  class CaloHit{
  public:
    CaloHit() {};
    ~CaloHit() { Clear(); };

    void Clear() { cellID=0; position.SetXYZ(0.,0.,0.); energy=-1; system = -1; module=-1; layer=-1; ParentShower=nullptr; }
    std::shared_ptr<CaloHit> Clone() const;

    void setOriginHit( edm4hep::CalorimeterHit& _hit ) { m_hit = _hit; }
    TVector3 getPosition() const { return position; }
    double   getEnergy() const { return energy; } 
    int getLayer() const {return layer;}
    int getModule() const {return module;}
    int getSystem() const {return system;}
    std::vector< std::pair<edm4hep::MCParticle, float> > getLinkedMCP() const { return MCParticleWeight; }
    edm4hep::MCParticle getLeadingMCP() const;
    float getLeadingMCPweight() const;
    edm4hep::CalorimeterHit getOriginHit() const { return m_hit; }

    void setcellID(unsigned long long _id) { cellID=_id; }
    void setcellID(int _s, int _m, int _l) { system=_s; module=_m; layer=_l; }
    void setEnergy(double _en) { energy=_en; }
    void setPosition( TVector3 _vec ) { position=_vec; }
    void setSystem(int _s) { system = _s; }
    void setModule(int _m ) { module = _m; }
    void setLayer(int _l) { layer = _l; }
    void setParentShower( Cyber::Calo2DCluster* _p ) { ParentShower=_p; }
    void addLinkedMCP( std::pair<edm4hep::MCParticle, float> _pair ) {MCParticleWeight.push_back(_pair); } 
    void setLinkedMCP( std::vector<std::pair<edm4hep::MCParticle, float>> _pairVec ) { MCParticleWeight.clear(); MCParticleWeight = _pairVec; }

  private: 
    int system; 
    int module; 
    int layer; 
    unsigned long long cellID; 
    TVector3 position;
    double   energy; 
    edm4hep::CalorimeterHit m_hit;
    Cyber::Calo2DCluster* ParentShower; 
    std::vector< std::pair<edm4hep::MCParticle, float> > MCParticleWeight;
  };

};
#endif
