#ifndef CALO_HIT_C
#define CALO_HIT_C

#include "Objects/CaloHit.h"

namespace PandoraPlus{


  std::shared_ptr<CaloHit> CaloHit::Clone() const{
    std::shared_ptr<CaloHit> m_hit = std::make_shared<CaloHit>();
    m_hit->setcellID(cellID);
    m_hit->setLayer(layer);
    m_hit->setPosition(position);
    m_hit->setEnergy(energy);
    m_hit->setParentShower(ParentShower);
    m_hit->setLinkedMCP(MCParticleWeight);
    return m_hit;
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
