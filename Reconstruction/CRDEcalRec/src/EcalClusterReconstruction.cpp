#include "EcalClusterReconstruction.h"

void EcalClusterReconstruction::Initialize(){

  //Initialize settings
  m_ESAlgSettings = new EnergySplittingAlg::Settings();
  m_ETAlgSettings = new EnergyTimeMatchingAlg::Settings();
  m_CCAlgSettings = new ConeClusteringAlg::Settings();

  //Initialize Algorthm
  m_energysplittingAlg = new EnergySplittingAlg();
  m_etmatchingAlg = new EnergyTimeMatchingAlg();
  m_coneclusterAlg = new ConeClusteringAlg();

 
}

StatusCode EcalClusterReconstruction::RunAlgorithm( EcalClusterReconstruction::Settings& settings,  PandoraPlusDataCol& dataCol ){

  //Iteration 1: Get initial information. 
  m_ESAlgSettings->SetInitialValue();
  m_ETAlgSettings->SetInitialValue();
  m_CCAlgSettings->SetInitialValue();

  std::cout<<"ESAlg::Running"<<std::endl;
  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol ); 
  std::cout<<"ETAlg::Running"<<std::endl;
  m_etmatchingAlg->     RunAlgorithm( *m_ETAlgSettings, dataCol );
  std::cout<<"CCAlg::Running"<<std::endl;
  m_coneclusterAlg->    RunAlgorithm( *m_CCAlgSettings, dataCol);
  std::cout<<"Iteration 0 end!"<<std::endl;


  //Iteration: 

  return StatusCode::SUCCESS;

}


StatusCode EcalClusterReconstruction::ClearAlgorithm(){

  delete m_energysplittingAlg;
  delete m_etmatchingAlg;
  delete m_coneclusterAlg;

  delete m_ESAlgSettings;
  delete m_ETAlgSettings;
  delete m_CCAlgSettings;

  return StatusCode::SUCCESS;

}


