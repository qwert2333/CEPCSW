#ifndef _BASICCLUSTERID_ALG_C
#define _BASICCLUSTERID_ALG_C

#include "Algorithm/"

void BasicClusterIDAlg::Settings::SetInitialValue(){

}

StatusCode BasicClusterIDAlg::Initialize(){

  return StatusCode::SUCCESS; 
}

StatusCode BasicClusterIDAlg::RunAlgorithm(BasicClusterIDAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){

  settings = m_settings; 

  


  return StatusCode::SUCCESS;
}

StatusCode BasicClusterIDAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS; 
}


#endif
