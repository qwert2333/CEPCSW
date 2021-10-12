#ifndef _BASICCLUSTERID_ALG_C
#define _BASICCLUSTERID_ALG_C

#include "Algorithm/BasicClusterIDAlg.h"

void BasicClusterIDAlg::Settings::SetInitialValue(){
  fl_departShower = false; 
}

StatusCode BasicClusterIDAlg::Initialize(){

  return StatusCode::SUCCESS; 
}

StatusCode BasicClusterIDAlg::RunAlgorithm(BasicClusterIDAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){
  settings = m_settings; 

  //Very simple cluster ID: only use stdDevE of the cluster
  for(int ic=0; ic<m_datacol.GoodClus3DCol.size(); ic++){
    m_datacol.GoodClus3DCol[ic].IdentifyCluster(); 
//    if(m_datacol.GoodClus3DCol[ic].getStdDevE()<0.05)      m_datacol.GoodClus3DCol[ic].setType(0); 
//    else if(m_datacol.GoodClus3DCol[ic].getStdDevE()>0.2)  m_datacol.GoodClus3DCol[ic].setType(1);
//    else                                                   m_datacol.GoodClus3DCol[ic].setType(2);
  }


  if(settings.fl_departShower){
  //separate 2Dshowers into different type
  m_datacol.ClearShower(); 
  std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_showerCol; m_showerCol.clear(); 
  std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_MIPshowerCol; m_MIPshowerCol.clear();
  std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_EMshowerCol;  m_EMshowerCol.clear(); 

  for(int ic=0; ic<m_datacol.GoodClus3DCol.size(); ic++){
    if(m_datacol.GoodClus3DCol[ic].getType()==0){
      std::vector<CRDEcalEDM::CRDCaloHitTransShower> tmp_col = m_datacol.GoodClus3DCol[ic].get2DShowers(); 
      m_MIPshowerCol.insert( m_MIPshowerCol.end(), tmp_col.begin(), tmp_col.end() );
    }
  
    else if(m_datacol.GoodClus3DCol[ic].getType()==1){
      std::vector<CRDEcalEDM::CRDCaloHitTransShower> tmp_col = m_datacol.GoodClus3DCol[ic].get2DShowers(); 
      m_EMshowerCol.insert( m_EMshowerCol.end(), tmp_col.begin(), tmp_col.end() );

    }
    else{
      std::vector<CRDEcalEDM::CRDCaloHitTransShower> tmp_col = m_datacol.GoodClus3DCol[ic].get2DShowers(); 
      m_showerCol.insert( m_showerCol.end(), tmp_col.begin(), tmp_col.end() );
    }

  }

  for(int ic=0; ic<m_datacol.BadClus3DCol.size(); ic++){
    std::vector<CRDEcalEDM::CRDCaloHitTransShower> tmp_col = m_datacol.BadClus3DCol[ic].get2DShowers(); 
    m_showerCol.insert( m_showerCol.end(), tmp_col.begin(), tmp_col.end() );
  }

  m_datacol.Shower2DCol = m_showerCol; 
  m_datacol.MIPShower2DCol = m_MIPshowerCol;
  m_datacol.EMShower2DCol = m_EMshowerCol; 

  }

  return StatusCode::SUCCESS;
}

StatusCode BasicClusterIDAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS; 
}


#endif
