#include "EcalClusterReconstruction.h"

void EcalClusterReconstruction::Initialize(){

  //Initialize settings
  m_ESAlgSettings = new EnergySplittingAlg::Settings();
  m_ETAlgSettings = new EnergyTimeMatchingAlg::Settings();
  m_CCAlgSettings = new ConeClusteringAlg::Settings();
  m_CMAlgSettings = new CandidateMakingAlg::Settings(); 
  m_CLAlgSettings = new ClusterMergingAlg::Settings(); 
  m_ACAlgSettings = new ArborClusteringAlg::Settings();

  //Initialize Algorthm
  m_energysplittingAlg = new EnergySplittingAlg();
  m_etmatchingAlg      = new EnergyTimeMatchingAlg();
  m_coneclusterAlg     = new ConeClusteringAlg();
  m_candidatemakingAlg = new CandidateMakingAlg();
  m_clustermergingAlg  = new ClusterMergingAlg();  
  m_arborclusteringAlg = new ArborClusteringAlg();
}

StatusCode EcalClusterReconstruction::RunAlgorithm( EcalClusterReconstruction::Settings& settings,  PandoraPlusDataCol& dataCol ){

  //Iteration 0: pre-treatment, get initial level candidates
  m_ESAlgSettings->SetInitialValue();
  m_ACAlgSettings->SetInitialValue();  
 
  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol ); 
  m_arborclusteringAlg->RunAlgorithm( *m_ACAlgSettings, dataCol );

  dataCol.PrintArborTree();

  return StatusCode::SUCCESS;

}


StatusCode EcalClusterReconstruction::ClearAlgorithm(){

  delete m_energysplittingAlg;
  delete m_etmatchingAlg;
  delete m_coneclusterAlg;
  delete m_candidatemakingAlg;
  delete m_clustermergingAlg; 
  delete m_arborclusteringAlg;

  delete m_ESAlgSettings;
  delete m_ETAlgSettings;
  delete m_CCAlgSettings;
  delete m_CMAlgSettings;
  delete m_CLAlgSettings; 
  delete m_ACAlgSettings;

  return StatusCode::SUCCESS;

}


