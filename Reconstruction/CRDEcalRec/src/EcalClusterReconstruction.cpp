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
  m_CMAlgSettings->SetInitialValue(); 
  //m_CMAlgSettings->SetDebug(1);

  dataCol.BlockVec_raw = dataCol.BlockVec;
  std::cout<<"Start Iteration 0"<<std::endl;
  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol ); 
  m_etmatchingAlg->     RunAlgorithm( *m_ETAlgSettings, dataCol );
  m_coneclusterAlg->    RunAlgorithm( *m_CCAlgSettings, dataCol);
  //std::cout<<"Iteration 0 is done!"<<std::endl;
  m_clustermergingAlg->RunAlgorithm( *m_CMAlgSettings, dataCol );
  dataCol.Flag_Iter++;


    dataCol.BlockVec_iter0       = dataCol.BlockVec;
    dataCol.LayerCol_iter0       = dataCol.LayerCol;
    dataCol.Shower2DCol_iter0    = dataCol.Shower2DCol;
    dataCol.GoodClus3DCol_iter0  = dataCol.GoodClus3DCol;
    dataCol.BadClus3DCol_iter0   = dataCol.BadClus3DCol;
    dataCol.Clus3DCol_iter0      = dataCol.Clus3DCol;

/*  //Iteration 1: 
  std::cout<<"Start Iteration 1"<<std::endl;
  dataCol.ClearLayer();
  dataCol.ClearShower();
  dataCol.ClearCluster();
  dataCol.ClearPFO();

  m_ESAlgSettings->SetValues();
  m_CCAlgSettings->SetGoodClusLevel(1);
  m_CCAlgSettings->SetMergeBadClus(false);

  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol );
  m_etmatchingAlg->     RunAlgorithm( *m_ETAlgSettings, dataCol );
  m_coneclusterAlg->    RunAlgorithm( *m_CCAlgSettings, dataCol ); 
  m_clustermergingAlg-> RunAlgorithm( *m_CMAlgSettings, dataCol );
  dataCol.Flag_Iter++;

    dataCol.BlockVec_iter1       = dataCol.BlockVec;
    dataCol.LayerCol_iter1       = dataCol.LayerCol;
    dataCol.Shower2DCol_iter1    = dataCol.Shower2DCol;
    dataCol.GoodClus3DCol_iter1  = dataCol.GoodClus3DCol;
    dataCol.BadClus3DCol_iter1   = dataCol.BadClus3DCol;
    dataCol.Clus3DCol_iter1      = dataCol.Clus3DCol;  


  //Itertion 2: 
  std::cout<<"Start Iteration 2"<<std::endl;
  dataCol.ClearLayer();
  dataCol.ClearShower();
  dataCol.ClearCluster();
  dataCol.ClearPFO();

  m_ESAlgSettings->SetValues();
  m_CCAlgSettings->SetGoodClusLevel(2);
  m_CCAlgSettings->SetMergeBadClus(true);

  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol );
  m_etmatchingAlg->     RunAlgorithm( *m_ETAlgSettings, dataCol );
  m_coneclusterAlg->    RunAlgorithm( *m_CCAlgSettings, dataCol );
  dataCol.Flag_Iter++;

  std::cout<<"Finish Iteration! "<<std::endl;
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


