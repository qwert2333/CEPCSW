#ifndef HOUGHSPACE_C
#define HOUGHSPACE_C

#include "Objects/HoughSpace.h"
#include <cmath>

namespace PandoraPlus{

  bool HoughSpace::isFromIP(PandoraPlus::HoughSpace::HoughHill& _hill) const{
    for(int i=0; i<_hill.getIndexAlpha().size(); ++i)
    {
      double alpha_low = m_sapceMap.GetXaxis()->GetBinLowEdge( _hill.getIndexAlpha()[i] );
      double alpha_high = alpha_low + m_sapceMap.GetXaxis()->GetBinWidth(_hill.getIndexAlpha()[i]);
      double rho_low = m_sapceMap.GetYaxis()->GetBinLowEdge( _hill.getIndexRho()[i] );
      double rho_high = rho_low +  m_sapceMap.GetYaxis()->GetBinWidth(_hill.getIndexRho()[i]);

      if( (IPHoughLine_outer.Eval(alpha_low)<rho_low  && IPHoughLine_outer.Eval(alpha_high)>rho_low) ||  
          (IPHoughLine_outer.Eval(alpha_low)<rho_high && IPHoughLine_outer.Eval(alpha_low)>rho_high) ||
          (IPHoughLine_inner.Eval(alpha_low)<rho_low  && IPHoughLine_inner.Eval(alpha_high)>rho_low) ||
          (IPHoughLine_inner.Eval(alpha_low)<rho_high && IPHoughLine_inner.Eval(alpha_low)>rho_high) )
        return true; 

      double alpha = m_sapceMap.GetXaxis()->GetBinCenter( _hill.getIndexAlpha()[i] );
      double rho = m_sapceMap.GetYaxis()->GetBinCenter( _hill.getIndexRho()[i] );
      if( (IPHoughLine_outer.Eval(alpha)>rho && IPHoughLine_inner.Eval(alpha)<rho) ||
          (IPHoughLine_inner.Eval(alpha)>rho && IPHoughLine_outer.Eval(alpha)<rho)   )
        return true;
    }
    return false;

  }


};
#endif
