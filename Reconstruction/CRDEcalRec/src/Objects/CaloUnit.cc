#ifndef CALOBAR_C
#define CALOBAR_C

#include "Objects/CaloUnit.h"
#include <cmath>

namespace PandoraPlus{

  bool CaloUnit::isNeighbor(const CaloUnit* x) const {
    if( cellID==x->getcellID() ) return false;
    if( system!=x->getSystem() || module!=x->getModule() || dlayer!=x->getDlayer() || slayer!=x->getSlayer() ) return false; 

    if( part==x->getPart() && stave==x->getStave() && fabs(bar - x->getBar())==1 ) return true;
    if( slayer==0 && stave==x->getStave() && fabs(part-x->getPart())<=1   && fabs(bar - x->getBar())<=1 ) return true;
    if( slayer==1 && part==x->getPart()   && fabs(stave-x->getStave())<=1 && fabs(bar - x->getBar())<=1 ) return true;

    if( isAtLowerEdgeZ()   && x->isAtUpperEdgeZ()   && x->getStave()==stave-1 && (x->getPart()-part)<=1 ) return true;
    if( isAtUpperEdgeZ()   && x->isAtLowerEdgeZ()   && x->getStave()==stave+1 && (x->getPart()-part)<=1 ) return true;
    if( isAtLowerEdgePhi() && x->isAtUpperEdgePhi() && x->getPart()==part-1   && (x->getStave()-stave)<=1 ) return true;
    if( isAtUpperEdgePhi() && x->isAtLowerEdgePhi() && x->getPart()==part+1   && (x->getStave()-stave)<=1 ) return true;

    return false;
  }


  bool CaloUnit::isAtLowerEdgePhi() const{
    return ( slayer==1 && bar==1 ); 
  }


  bool CaloUnit::isAtUpperEdgePhi() const{
    return ( slayer==1 && bar==NbarPhi ); 
  }


  bool CaloUnit::isAtLowerEdgeZ() const{
    return ( slayer==0 && bar==1 );
  }


  bool CaloUnit::isAtUpperEdgeZ() const{
    return ( slayer==0 && bar==NbarZ );
  }

  bool CaloUnit::isModuleAdjacent( const CaloUnit* x ) const{
    if(module==x->getModule() || (fabs(module-x->getModule())>1 && fabs(module-x->getModule())!=Nmodule-1 ) ) return false; 

    int dlayer_lo, slayer_lo, part_lo, stave_lo, bar_lo;
    int dlayer_hi, slayer_hi, part_hi, stave_hi, bar_hi;
    if( module > x->getModule() ) {
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


  CaloUnit* CaloUnit::Clone() const{
    CaloUnit* m_bar = new CaloUnit();
    m_bar->setcellID(cellID);
    m_bar->setcellID( system, module, stave, dlayer, part, slayer, bar );
    m_bar->setPosition(position);
    m_bar->setQ(Q1, Q2);
    m_bar->setT(T1, T2);
    return m_bar;
  }

};
#endif
