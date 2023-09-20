#ifndef _CLUSTERMERGING_ALG_C
#define _CLUSTERMERGING_ALG_C

#include "Algorithm/ClusterMergingAlg.h"

StatusCode ClusterMergingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters

  return StatusCode::SUCCESS;
};

StatusCode ClusterMergingAlg::Initialize( PandoraPlusDataCol& m_datacol ){

  return StatusCode::SUCCESS;
};

StatusCode ClusterMergingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){


  return StatusCode::SUCCESS;
};

StatusCode ClusterMergingAlg::ClearAlgorithm(){


  return StatusCode::SUCCESS;
};


#endif
