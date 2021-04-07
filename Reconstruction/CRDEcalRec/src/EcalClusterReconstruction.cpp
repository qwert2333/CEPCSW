#include "EcalClusterReconstruction.h"

void EcalClusterReconstruction::Initialize(){

  //Initialize settings
  m_ESAlgSettings = new EnergySplittingAlg::Settings();
  m_ETAlgSettings = new EnergyTimeMatchingAlg::Settings();
  m_CCAlgSettings = new ConeClusteringAlg::Settings();
  m_CMAlgSettings = new ClusterMergingAlg::Settings(); 

  //Initialize Algorthm
  m_energysplittingAlg = new EnergySplittingAlg();
  m_etmatchingAlg      = new EnergyTimeMatchingAlg();
  m_coneclusterAlg     = new ConeClusteringAlg();
  m_clustermergingAlg  = new ClusterMergingAlg();
 
  m_clustermergingAlg->Initialize();

}

StatusCode EcalClusterReconstruction::RunAlgorithm( EcalClusterReconstruction::Settings& settings,  PandoraPlusDataCol& dataCol ){

  //Iteration 0: Get initial information. 
  m_ESAlgSettings->SetInitialValue();
  m_ETAlgSettings->SetInitialValue();
  m_CCAlgSettings->SetInitialValue();

  //std::cout<<"ESAlg::Running"<<std::endl;
  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol ); 
  //std::cout<<"ETAlg::Running"<<std::endl;
  m_etmatchingAlg->     RunAlgorithm( *m_ETAlgSettings, dataCol );
  //std::cout<<"CCAlg::Running"<<std::endl;
  m_coneclusterAlg->    RunAlgorithm( *m_CCAlgSettings, dataCol);
  //std::cout<<"Iteration 0 is done!"<<std::endl;


  //dataCol.PrintLayer();
  //dataCol.PrintShower();
  //dataCol.Print3DClus();

  //Iteration 1: 
  m_CMAlgSettings->SetInitialValue(); 
  m_clustermergingAlg->RunAlgorithm( *m_CMAlgSettings, dataCol );

  //std::cout<<"DEBUG: come back! iter = "<<dataCol.Flag_Iter<<std::endl;

  if(dataCol.Flag_Iter){
    dataCol.ClearLayer();
    dataCol.ClearShower();
    dataCol.ClearCluster();
    dataCol.ClearPFO();

    m_ESAlgSettings->SetValues();
    m_CCAlgSettings->SetMergeBadClus(true);

    m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol );
    m_etmatchingAlg->     RunAlgorithm( *m_ETAlgSettings, dataCol );
    m_coneclusterAlg->    RunAlgorithm( *m_CCAlgSettings, dataCol); 
  }
  //std::cout<<"Iteration 1 is done!"<<std::endl;


  //Iteration 2: 
/*  m_CMAlgSettings->SetInitialValue();
  m_clustermergingAlg->RunAlgorithm( *m_CMAlgSettings, dataCol );
  if(dataCol.Flag_Iter){
    dataCol.ClearLayer();
    dataCol.ClearShower();
    dataCol.ClearCluster();
    dataCol.ClearPFO();

    m_ESAlgSettings->SetValues();
    m_CCAlgSettings->SetMergeBadClus(true);

    m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol );
    m_etmatchingAlg->     RunAlgorithm( *m_ETAlgSettings, dataCol );
    m_coneclusterAlg->    RunAlgorithm( *m_CCAlgSettings, dataCol);
  }
  //std::cout<<"Iteration 2 is done!"<<std::endl;
*/

  return StatusCode::SUCCESS;

}


StatusCode EcalClusterReconstruction::ClearAlgorithm(){

  delete m_energysplittingAlg;
  delete m_etmatchingAlg;
  delete m_coneclusterAlg;
  delete m_clustermergingAlg;

  delete m_ESAlgSettings;
  delete m_ETAlgSettings;
  delete m_CCAlgSettings;
  delete m_CMAlgSettings;

  return StatusCode::SUCCESS;

}


