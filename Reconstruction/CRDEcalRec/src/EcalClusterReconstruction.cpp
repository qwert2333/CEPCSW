#include "EcalClusterReconstruction.h"

void EcalClusterReconstruction::Initialize(){

  //Initialize settings
  m_ESAlgSettings = new EnergySplittingAlg::Settings();
  m_ETAlgSettings = new EnergyTimeMatchingAlg::Settings();
  m_CCAlgSettings = new ConeClusteringAlg::Settings();
  m_CMAlgSettings = new CandidateMakingAlg::Settings(); 
  m_CLAlgSettings = new ClusterMergingAlg::Settings(); 

  //Initialize Algorthm
  m_energysplittingAlg = new EnergySplittingAlg();
  m_etmatchingAlg      = new EnergyTimeMatchingAlg();
  m_coneclusterAlg     = new ConeClusteringAlg();
  m_candidatemakingAlg = new CandidateMakingAlg();
  m_clustermergingAlg  = new ClusterMergingAlg();  

}

StatusCode EcalClusterReconstruction::RunAlgorithm( EcalClusterReconstruction::Settings& settings,  PandoraPlusDataCol& dataCol ){

  //Iteration 0: pre-treatment, get initial level candidates
  m_ESAlgSettings->SetInitialValue();
  m_ETAlgSettings->SetInitialValue();
  m_CCAlgSettings->SetInitialValue();
  m_CMAlgSettings->SetInitialValue(); 
  m_CLAlgSettings->SetInitialValue(); 
  m_CMAlgSettings->SetDebug(0);

  std::cout<<"Start Iteration 0"<<std::endl;
  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol ); 
  m_etmatchingAlg->     RunAlgorithm( *m_ETAlgSettings, dataCol );
  m_coneclusterAlg->    RunAlgorithm( *m_CCAlgSettings, dataCol);
  m_candidatemakingAlg->RunAlgorithm( *m_CMAlgSettings, dataCol );

    dataCol.BlockVec_iter0       = dataCol.BlockVec;
    dataCol.LayerCol_iter0       = dataCol.LayerCol;
    dataCol.Shower2DCol_iter0    = dataCol.Shower2DCol;
    dataCol.GoodClus3DCol_iter0  = dataCol.GoodClus3DCol;
    dataCol.BadClus3DCol_iter0   = dataCol.BadClus3DCol;
    dataCol.Clus3DCol_iter0      = dataCol.Clus3DCol;

//for(int ib=0; ib<dataCol.BlockVec.size(); ib++){ cout<<"  Block #"<<ib<<'\n';  dataCol.BlockVec[ib].PrintCandidates(); }

  //Iteration 1: Add tracks to get candidates. 
  dataCol.ClearLayer();
  dataCol.ClearShower();
  dataCol.ClearCluster();
  dataCol.ClearPFO();

  std::cout<<"Start Iteration 1"<<std::endl;
  m_ESAlgSettings->SetValues();
  //m_ESAlgSettings->SetEseedRel(0.45);
  m_ESAlgSettings->SetDebug(0);
  m_ESAlgSettings->SetUseCandidate(true);
  m_ETAlgSettings->SetUseCandidate(1);
  m_ETAlgSettings->SetUseChi2(false);
  m_CCAlgSettings->SetGoodClusLevel(2);
  m_CLAlgSettings->SetMergeGoodCluster(true);
  m_CLAlgSettings->SetMergeBadCluster(false);
  m_CMAlgSettings->SetUseTrk(true);

  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol );
  m_etmatchingAlg->     RunAlgorithm( *m_ETAlgSettings, dataCol );
  m_coneclusterAlg->    RunAlgorithm( *m_CCAlgSettings, dataCol ); 
  m_clustermergingAlg-> RunAlgorithm( *m_CLAlgSettings, dataCol );
  m_candidatemakingAlg->RunAlgorithm( *m_CMAlgSettings, dataCol );

    dataCol.BlockVec_iter1       = dataCol.BlockVec;
    dataCol.LayerCol_iter1       = dataCol.LayerCol;
    dataCol.Shower2DCol_iter1    = dataCol.Shower2DCol;
    dataCol.GoodClus3DCol_iter1  = dataCol.GoodClus3DCol;
    dataCol.BadClus3DCol_iter1   = dataCol.BadClus3DCol;
    dataCol.Clus3DCol_iter1      = dataCol.Clus3DCol;

dataCol.PrintLayer(); 
dataCol.PrintShower(); 
dataCol.Print3DClus(); 
for(int ib=0; ib<dataCol.BlockVec.size(); ib++){ cout<<"  Block #"<<ib<<'\n';  dataCol.BlockVec[ib].PrintCandidates(); }

  //Iteration 2: Use all candidates for reconstruction. 
  std::cout<<"Start Iteration 2"<<std::endl;
  dataCol.ClearLayer();
  dataCol.ClearShower();
  dataCol.ClearCluster();
  dataCol.ClearPFO();

  m_ESAlgSettings->SetValues();
  //m_ESAlgSettings->SetEseedRel(0.45);
  m_ESAlgSettings->SetUseCandidate(true);
  m_ETAlgSettings->SetUseCandidate(0); 
  m_ETAlgSettings->SetUseChi2(true);    //Not use candidate, only use chi2 for matching. (Old way)
  m_CCAlgSettings->SetConeValue( PI/6., 30., PI/6., 30.);
  m_CCAlgSettings->SetUseCandidate(1);
  m_CCAlgSettings->SetGoodClusLevel(3);
  m_CLAlgSettings->SetMergeGoodCluster(true);
  m_CLAlgSettings->SetMergeBadCluster(true);

cout<<"EnergySplittingAlg"<<endl;
  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol );
cout<<"ETMatchingAlg"<<endl;
  m_etmatchingAlg->     RunAlgorithm( *m_ETAlgSettings, dataCol );
cout<<"ConeClusteringAlg"<<endl;
  m_coneclusterAlg->    RunAlgorithm( *m_CCAlgSettings, dataCol );
cout<<"ClsuterMergingAlg"<<endl;
  m_clustermergingAlg-> RunAlgorithm( *m_CLAlgSettings, dataCol );
cout<<"End Reconstruction"<<endl;

    dataCol.BlockVec_iter2       = dataCol.BlockVec;
    dataCol.LayerCol_iter2       = dataCol.LayerCol;
    dataCol.Shower2DCol_iter2    = dataCol.Shower2DCol;
    dataCol.GoodClus3DCol_iter2  = dataCol.GoodClus3DCol;
    dataCol.BadClus3DCol_iter2   = dataCol.BadClus3DCol;
    dataCol.Clus3DCol_iter2      = dataCol.Clus3DCol;  


dataCol.PrintLayer(); 
dataCol.PrintShower(); 
dataCol.Print3DClus(); 
//for(int ib=0; ib<dataCol.BlockVec.size(); ib++){ cout<<"  Block #"<<ib<<'\n';  dataCol.BlockVec[ib].PrintCandidates(); }

  std::cout<<"Finish Iteration! "<<std::endl;

  return StatusCode::SUCCESS;

}


StatusCode EcalClusterReconstruction::ClearAlgorithm(){

  delete m_energysplittingAlg;
  delete m_etmatchingAlg;
  delete m_coneclusterAlg;
  delete m_candidatemakingAlg;
  delete m_clustermergingAlg; 

  delete m_ESAlgSettings;
  delete m_ETAlgSettings;
  delete m_CCAlgSettings;
  delete m_CMAlgSettings;
  delete m_CLAlgSettings; 

  return StatusCode::SUCCESS;

}


