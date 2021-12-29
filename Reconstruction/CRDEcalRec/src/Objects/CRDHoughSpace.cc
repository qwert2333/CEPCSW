#ifndef _CRD_CALOHOUGHSPACE_C
#define _CRD_CALOHOUGHSPACE_C

#include "Objects/CRDHoughSpace.h"

namespace CRDEcalEDM{

  bool CRDHoughSpace::HoughHill::isNeighbor( CRDEcalEDM::CRDHoughSpace::HoughCell _cell ){
    for(int ic=0; ic<m_cells.size(); ic++){
      if( fabs(m_cells[ic].getIndexAlpha()-_cell.getIndexAlpha())<=1 && 
          fabs(m_cells[ic].getIndexRho() - _cell.getIndexRho())<=1 )
      return true; 
    }
    return false; 
  }

  CRDEcalEDM::CRDCaloHitLongiCluster CRDHoughSpace::HoughHill::TransformToCluster(){
    CRDEcalEDM::CRDCaloHitLongiCluster m_clus; m_clus.Clear();
    for(int ic=0; ic<m_cells.size(); ic++){
    for(int iobj=0; iobj<m_cells[ic].getObjects().size(); iobj++){
      CRDEcalEDM::CRDHoughObject m_obj = m_cells[ic].getObjects()[iobj];
      m_clus.AddBarShower( *m_obj.getLocalMax() );
    }}
    m_clus.FitAxis();
    return m_clus; 
  }

};
#endif
