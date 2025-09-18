#ifndef CALO_HIT_C
#define CALO_HIT_C

#include "Objects/CaloHit.h"

namespace Cyber{


  std::shared_ptr<CaloHit> CaloHit::Clone() const{
    std::shared_ptr<CaloHit> _newhit = std::make_shared<CaloHit>();
    _newhit->setcellID(cellID);
    _newhit->setLayer(layer);
    _newhit->setPosition(position);
    _newhit->setEnergy(energy);
    _newhit->setParentShower(ParentShower);
    _newhit->setLinkedMCP(MCParticleWeight);
    edm4hep::CalorimeterHit originhit = getOriginHit();
    _newhit->setOriginHit( originhit );
    return _newhit;
  }

  edm4hep::MCParticle CaloHit::getLeadingMCP() const{
    float maxWeight = -1.;
    edm4hep::MCParticle mcp;
    for(auto& iter: MCParticleWeight){
      if(iter.second>maxWeight){
        mcp = iter.first;
        maxWeight = iter.second;
      }
    }

    return mcp;
  }

  float CaloHit::getLeadingMCPweight() const{
    float maxWeight = -1.;
    for(auto& iter: MCParticleWeight){
      if(iter.second>maxWeight){
        maxWeight = iter.second;
      }
    }
    return maxWeight;
  }

};
#endif
