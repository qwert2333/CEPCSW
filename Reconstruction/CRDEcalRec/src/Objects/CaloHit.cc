#ifndef CALO_HIT_H
#define CALO_HIT_H

#include "Objects/CaloHit.h"

namespace PandoraPlus{

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
