#ifndef CALOBAR_C
#define CALOBAR_C

#include "Objects/CaloUnit.h"
#include <cmath>

namespace Cyber{

  bool CaloUnit::isNeighbor(const CaloUnit* x) const {
    if( cellID==x->getcellID() ) return false;
    if( system!=x->getSystem() || dlayer!=x->getDlayer() || slayer!=x->getSlayer() ) return false; 
    if( system==System_Endcap){ // endcap
      if( module!=x->getModule() ) return false; 
      if( slayer==0){
        if(isAtLowerEdgeZ()){
          if( ( ( getPosition().x()-x->getPosition().x() <= ecal_endcap_barsize + 2*ecal_endcap_deadarea && getPosition().x()-x->getPosition().x()>1e-6)
                || ( x->getPosition().x()-getPosition().x() <= ecal_endcap_barsize + 1e-6 && x->getPosition().x()-getPosition().x()>1e-6) )
              && fabs(getPosition().y()-x->getPosition().y()) <= (getBarLength() + x->getBarLength())/2. + 2*ecal_endcap_deadarea + ecal_endcap_barsize) return true;
        }
        else if(isAtUpperEdgeZ()){
          if( ( (getPosition().x()-x->getPosition().x() <= ecal_endcap_barsize + 1e-6 && getPosition().x()-x->getPosition().x()>1e-6) 
                || (x->getPosition().x()-getPosition().x() <= ecal_endcap_barsize + 2*ecal_endcap_deadarea && x->getPosition().x()-getPosition().x()>1e-6) )
              && fabs(getPosition().y()-x->getPosition().y()) <= (getBarLength() + x->getBarLength())/2. + 2*ecal_endcap_deadarea + ecal_endcap_barsize) return true;
        }
        else{
          if(fabs(getPosition().x()-x->getPosition().x()) <= ecal_endcap_barsize + 1e-6 
          && fabs(getPosition().y()-x->getPosition().y()) <= (getBarLength() + x->getBarLength())/2. + 2*ecal_endcap_deadarea + 2*ecal_endcap_barsize) return true;
        }
      } 
      if( slayer==1) {
        if(isAtLowerEdgePhi()){
          if( ( (getPosition().y()-x->getPosition().y() <= ecal_endcap_barsize + 2*ecal_endcap_deadarea && getPosition().y()-x->getPosition().y()>1e-6) 
                || ( x->getPosition().y()-getPosition().y() <= ecal_endcap_barsize + 1e-6 && x->getPosition().y()-getPosition().y()>1e-6) )
              && fabs(getPosition().x()-x->getPosition().x()) <= (getBarLength() + x->getBarLength())/2. + 2*ecal_endcap_deadarea + ecal_endcap_barsize) return true;
        }
        else if(isAtUpperEdgePhi()){
          if( ( (getPosition().y()-x->getPosition().y() <= ecal_endcap_barsize + 1e-6 && getPosition().y()-x->getPosition().y()>1e-6) 
                || (x->getPosition().y()-getPosition().y() <= ecal_endcap_barsize + 2*ecal_endcap_deadarea && x->getPosition().y()-getPosition().y()>1e-6) )
              && fabs(getPosition().x()-x->getPosition().x()) <= (getBarLength() + x->getBarLength())/2. + 2*ecal_endcap_deadarea + ecal_endcap_barsize) return true;
        }
        else{
          if(fabs(getPosition().y()-x->getPosition().y()) <= ecal_endcap_barsize + 1e-6 
          && fabs(getPosition().x()-x->getPosition().x()) <= (getBarLength() + x->getBarLength())/2. + 2*ecal_endcap_deadarea + 2*ecal_endcap_barsize) return true;
        }
      } 
    }
    else{ //barrel
      if( module==x->getModule() && stave==x->getStave() && fabs(bar - x->getBar())==1 ) return true;
      if( slayer==0 && stave==x->getStave() && (fabs(module-x->getModule())<=1 || fabs(module-x->getModule())==Nmodule-1 ) && fabs(bar - x->getBar())<=1 ) return true;
      if( slayer==1 && module==x->getModule() && fabs(stave-x->getStave())<=1 && fabs(bar - x->getBar())<=1 ) return true;

      if( isAtLowerEdgeZ()   && x->isAtUpperEdgeZ()   && x->getStave()==stave-1 && (fabs(x->getModule()-module)<=1 || fabs(x->getModule()-module)==Nmodule-1 ) ) return true;
      if( isAtUpperEdgeZ()   && x->isAtLowerEdgeZ()   && x->getStave()==stave+1 && (fabs(x->getModule()-module)<=1 || fabs(x->getModule()-module)==Nmodule-1 ) ) return true;
      if( isAtLowerEdgePhi() && x->isAtUpperEdgePhi() && (x->getModule()==module-1 || (x->getModule() - module)==Nmodule-1 ) && fabs(x->getStave()-stave)<=1 ) return true;
      if( isAtUpperEdgePhi() && x->isAtLowerEdgePhi() && (x->getModule()==module+1 || (module - x->getModule())==Nmodule-1 ) && fabs(x->getStave()-stave)<=1 ) return true;
    }
    return false;
  }

  bool CaloUnit::isLongiNeighbor(const CaloUnit* x) const {
    if( cellID==x->getcellID() ) return false;
    if( system!=x->getSystem() || slayer!=x->getSlayer() ) return false; 
    if( abs(dlayer-x->getDlayer())!=1) return false;
    if( system==System_Endcap ){ // endcap
      if( slayer==0 ){
        if(fabs(getPosition().x()-x->getPosition().x()) <= ecal_endcap_barsize + 1e-6 
          && fabs(getPosition().y()-x->getPosition().y()) <= (getBarLength() + x->getBarLength())/2. + 2*ecal_endcap_deadarea + 2*ecal_endcap_barsize) return true;
      }
      if( slayer==1 ){
        if(fabs(getPosition().y()-x->getPosition().y()) <= ecal_endcap_barsize + 1e-6 
          && fabs(getPosition().x()-x->getPosition().x()) <= (getBarLength() + x->getBarLength())/2. + 2*ecal_endcap_deadarea + 2*ecal_endcap_barsize) return true;
      }
    }
    else{
      if( slayer==0 )
      {
        if( stave==x->getStave() && (fabs(module-x->getModule())<=1 || fabs(module-x->getModule())==Nmodule-1 ) && fabs(bar - x->getBar())<=1 ) return true;
        if( isAtLowerEdgeZ()   && x->isAtUpperEdgeZ()   && x->getStave()==stave-1 && (fabs(x->getModule()-module)<=1 || fabs(x->getModule()-module)==Nmodule-1 ) ) return true;
        if( isAtUpperEdgeZ()   && x->isAtLowerEdgeZ()   && x->getStave()==stave+1 && (fabs(x->getModule()-module)<=1 || fabs(x->getModule()-module)==Nmodule-1 ) ) return true;
      }
      if( slayer==1 )
      {
        if(module==x->getModule() && fabs(stave-x->getStave())<=1 && fabs(bar - x->getBar())<=1 ) return true;
        if( isAtLowerEdgePhi() && x->isAtUpperEdgePhi() && (x->getModule()==module-1 || (x->getModule() - module)==Nmodule-1 ) && fabs(x->getStave()-stave)<=1 ) return true;
        if( isAtUpperEdgePhi() && x->isAtLowerEdgePhi() && (x->getModule()==module+1 || (module - x->getModule())==Nmodule-1 ) && fabs(x->getStave()-stave)<=1 ) return true;
      }
    }

    
    return false;
  }

  bool CaloUnit::isAtLowerEdgePhi() const{
    return ( slayer==1 && bar==0 ); 
  }


  bool CaloUnit::isAtUpperEdgePhi() const{
    bool isEdge = false; 
    if( system==System_Endcap ){
      if(slayer==1 && bar==nBarInLayer-1) isEdge = true;
    }
    else{
      if(module%2==0){
        if(slayer==1 && bar==NbarPhi_even[dlayer-1]-1) isEdge = true;
      }
      else{
        if(slayer==1 && bar==NbarPhi_odd[dlayer-1]-1) isEdge = true;
      }
    }
    return isEdge; 
  }


  bool CaloUnit::isAtLowerEdgeZ() const{
    return ( slayer==0 && bar==0 );
  }


  bool CaloUnit::isAtUpperEdgeZ() const{
    if( system==System_Endcap ){
      return ( slayer==0 && bar==nBarInLayer-1 );
    }
    else{
      return ( slayer==0 && bar==NbarZ-1 );
    } 
  }

/*
  bool CaloUnit::isModuleAdjacent( const CaloUnit* x ) const{

    if(module==x->getModule()) return false;
    if( fabs(module-x->getModule())>1 && abs(module-x->getModule())!=Nmodule-1 ) return false; 

    int dlayer_lo, slayer_lo, part_lo, stave_lo, bar_lo;
    int dlayer_hi, slayer_hi, part_hi, stave_hi, bar_hi;
    if( module - x->getModule()==1 || x->getModule()-module==Nmodule-1 ) {
      dlayer_lo=x->getDlayer(); slayer_lo=x->getSlayer(); part_lo=x->getPart(); stave_lo=x->getStave(); bar_lo=x->getBar(); 
      dlayer_hi=dlayer; slayer_hi=slayer; part_hi=part; stave_hi=stave; bar_hi=bar; 
    }
    else{
      dlayer_lo=dlayer; slayer_lo=slayer; part_lo=part; stave_lo=stave; bar_lo=bar;
      dlayer_hi=x->getDlayer(); slayer_hi=x->getSlayer(); part_hi=x->getPart(); stave_hi=x->getStave(); bar_hi=x->getBar();
    }


    if( dlayer_lo!=1 || slayer_lo!=0 || part_lo!=Npart || part_hi!=1 ) return false;
    if( slayer_hi==0 ){
      if( (stave_lo==stave_hi && abs(bar_lo-bar_hi)<=1) || 
          (abs(stave_lo-stave_hi)==1 && isAtLowerEdgeZ() && x->isAtUpperEdgeZ() ) ||
          (abs(stave_lo-stave_hi)==1 && isAtUpperEdgeZ() && x->isAtLowerEdgeZ() )  ) return true; 
    }
    else if( slayer_hi==1 && stave_lo==stave_hi && bar_hi==1 ) return true;

    return false; 
  }
*/

  std::shared_ptr<CaloUnit> CaloUnit::Clone() const{
    std::shared_ptr<CaloUnit> m_bar = std::make_shared<CaloUnit>();
    m_bar->setcellID(cellID);
    m_bar->setcellID( system, module, stave, part, dlayer, slayer, bar );
    m_bar->setPosition(position);
    m_bar->setQ(Q1, Q2);
    m_bar->setT(T1, T2);
    for(auto ilink : MCParticleWeight) m_bar->addLinkedMCP(ilink);
    return m_bar;
  }

  edm4hep::MCParticle CaloUnit::getLeadingMCP() const{
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

  float CaloUnit::getLeadingMCPweight() const{
    float maxWeight = -1.;
    for(auto& iter: MCParticleWeight){
      if(iter.second>maxWeight){
        maxWeight = iter.second;
      }
    }
    return maxWeight;
  }

/*  bool CaloUnit::isLongiModuleAdjacent( const CaloUnit* x ) const{

    if(module==x->getModule()) return false;
    if( fabs(module-x->getModule())>1 && abs(module-x->getModule())!=Nmodule-1 ) return false; 

    int dlayer_lo, slayer_lo, part_lo, stave_lo, bar_lo;
    int dlayer_hi, slayer_hi, part_hi, stave_hi, bar_hi;

    if( module - x->getModule()==1 || x->getModule()-module==Nmodule-1 ) {
      dlayer_lo=x->getDlayer(); slayer_lo=x->getSlayer(); part_lo=x->getPart(); stave_lo=x->getStave(); bar_lo=x->getBar(); 
      dlayer_hi=dlayer; slayer_hi=slayer; part_hi=part; stave_hi=stave; bar_hi=bar; 
    }
    else{
      dlayer_lo=dlayer; slayer_lo=slayer; part_lo=part; stave_lo=stave; bar_lo=bar;
      dlayer_hi=x->getDlayer(); slayer_hi=x->getSlayer(); part_hi=x->getPart(); stave_hi=x->getStave(); bar_hi=x->getBar();
    }

    // if( dlayer_lo!=1 || slayer_lo!=0 || part_lo!=Npart || part_hi!=1 ) return false;
    // if( slayer_hi==0 ){
    //   if( (stave_lo==stave_hi && abs(bar_lo-bar_hi)<=1) || 
    //       (abs(stave_lo-stave_hi)==1 && isAtLowerEdgeZ() && x->isAtUpperEdgeZ() ) ||
    //       (abs(stave_lo-stave_hi)==1 && isAtUpperEdgeZ() && x->isAtLowerEdgeZ() )  ) return true; 
    // }
    // else if( slayer_hi==1 && stave_lo==stave_hi && bar_hi==1 ) return true;
    
    if(slayer_lo==0)
    {
      if(part_lo==Npart && part_hi==1 && dlayer_lo==1)
      {
        if( (stave_lo==stave_hi && abs(bar_lo-bar_hi)<=1) || 
          (abs(stave_lo-stave_hi)==1 && isAtLowerEdgeZ() && x->isAtUpperEdgeZ() ) ||
          (abs(stave_lo-stave_hi)==1 && isAtUpperEdgeZ() && x->isAtLowerEdgeZ() )  ) return true; 
      } 
    }
    else
    {
      if(part_lo==Npart && part_hi==1 && dlayer_lo==1)
      {
        if(bar_hi==1 && abs(stave_lo-stave_hi)<=1)
        {
          if(abs(bar_lo-over_module[(dlayer_hi-1)*2])<=over_module_set || abs(bar_lo-over_module[(dlayer_hi-1)*2+1])<=over_module_set)
          {
            return true;
          }
        } 
        //method1: abs(dlayer_hi*2*sqrt(2)-(bar_lo-14+1))<=3
        //method2: abs(bar_lo - (14+(dlayer_hi-1)*3))<=3
      }
    }
    return false; 
  }
*/  
};
#endif
