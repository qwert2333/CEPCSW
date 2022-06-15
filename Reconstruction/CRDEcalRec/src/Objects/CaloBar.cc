#ifndef CALOBAR_C
#define CALOBAR_C

#include "Objects/CaloBar.h"
#include <cmath>

namespace PandoraPlus{

  bool CaloBar::isNeighbor(const CaloBar* x) const {
    if(x->getcellID() != cellID &&
       x->getSystem() == system &&
       x->getModule() == module &&
       x->getDlayer() == dlayer &&
       x->getPart()   == part   &&
       x->getStave()  == stave  &&
       x->getSlayer() == slayer &&
       (x->getBar()==bar+1 || x->getBar()==bar-1)
      ) return true;
      return false;
  }

  CaloBar* CaloBar::Clone() const{
    CaloBar* m_bar = new CaloBar();
    m_bar->setcellID(cellID);
    m_bar->setcellID( system, module, stave, dlayer, part, slayer, bar );
    m_bar->setPosition(position);
    m_bar->setQ(Q1, Q2);
    m_bar->setT(T1, T2);
    return m_bar;
  }

};
#endif
