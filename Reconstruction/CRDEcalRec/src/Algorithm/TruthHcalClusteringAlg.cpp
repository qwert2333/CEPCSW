#ifndef _TRUTHHCALCLUS_ALG_C
#define _TRUTHHCALCLUS_ALG_C

#include "Algorithm/TruthHcalClusteringAlg.h"

StatusCode TruthHcalClusteringAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters

  return StatusCode::SUCCESS;
};

StatusCode TruthHcalClusteringAlg::Initialize( PandoraPlusDataCol& m_datacol ){

  return StatusCode::SUCCESS;
};

StatusCode TruthHcalClusteringAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  return StatusCode::SUCCESS;
};

StatusCode TruthHcalClusteringAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
};



#endif
