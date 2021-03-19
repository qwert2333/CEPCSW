#ifndef _CLUSTERMERGING_ALG_C
#define _CLUSTERMERGING_ALG_C

#include "Algorithm/ClusterMergingAlg.h"

void ClusterMergingAlg::Settings::SetInitialValue(){

}

StatusCode ClusterMergingAlg::Initialize(){

  return StatusCode::SUCCESS;
}

StatusCode ClusterMergingAlg::RunAlgorithm(ClusterMergingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){

  std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_goodClus = m_datacol.GoodClus3DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_badClus  = m_datacol.BadClus3DCol;


  HaveGhostHits( m_goodClus, m_badClus  );



}


bool ClusterMergingAlg::HaveGhostHits(  )

#endif
