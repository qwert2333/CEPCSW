#include "EcalClusterReconstruction.h"

void EcalClusterReconstruction::Initialize(){

  //Initialize settings
  m_ESAlgSettings = new EnergySplittingAlg::Settings();
  m_ETAlgSettings = new EnergyTimeMatchingAlg::Settings();
  m_CCAlgSettings = new ConeClusteringAlg::Settings();
  m_CMAlgSettings = new CandidateMakingAlg::Settings(); 
  m_CLAlgSettings = new ClusterMergingAlg::Settings(); 
  m_ACAlgSettings = new ArborClusteringAlg::Settings();
  m_AMAlgSettings = new ArborTreeMergingAlg::Settings(); 

  //Initialize Algorthm
  m_energysplittingAlg  = new EnergySplittingAlg();
  m_etmatchingAlg       = new EnergyTimeMatchingAlg();
  m_coneclusterAlg      = new ConeClusteringAlg();
  m_candidatemakingAlg  = new CandidateMakingAlg();
  m_clustermergingAlg   = new ClusterMergingAlg();  
  m_arborclusteringAlg  = new ArborClusteringAlg();
  m_arbortreemergingAlg = new ArborTreeMergingAlg(); 
}

StatusCode EcalClusterReconstruction::RunAlgorithm( EcalClusterReconstruction::Settings& settings,  PandoraPlusDataCol& dataCol ){

/*  //Iteration 0: pre-treatment, get initial level candidates
  m_ESAlgSettings->SetInitialValue();
  m_ETAlgSettings->SetInitialValue();
  m_CCAlgSettings->SetInitialValue();
  m_CMAlgSettings->SetInitialValue();
  m_CLAlgSettings->SetInitialValue();
  m_CMAlgSettings->SetDebug(0);

  m_CCAlgSettings->SetConeValue( PI/4., 50., PI/6., 30.);

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
  m_CCAlgSettings->SetConeValue( PI/4., 40., PI/6., 30.);
  m_CCAlgSettings->SetUseCandidate(1);
  m_CCAlgSettings->SetGoodClusLevel(3);
  m_CLAlgSettings->SetMergeGoodCluster(true);
  m_CLAlgSettings->SetMergeBadCluster(true);


  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol );
  m_etmatchingAlg->     RunAlgorithm( *m_ETAlgSettings, dataCol );
  m_coneclusterAlg->    RunAlgorithm( *m_CCAlgSettings, dataCol );
  //m_clustermergingAlg-> RunAlgorithm( *m_CLAlgSettings, dataCol );

    dataCol.BlockVec_iter2       = dataCol.BlockVec;
    dataCol.LayerCol_iter2       = dataCol.LayerCol;
    dataCol.Shower2DCol_iter2    = dataCol.Shower2DCol;
    dataCol.GoodClus3DCol_iter2  = dataCol.GoodClus3DCol;
    dataCol.BadClus3DCol_iter2   = dataCol.BadClus3DCol;
    dataCol.Clus3DCol_iter2      = dataCol.Clus3DCol;
*/


  //Iteration 0: pre-treatment, get initial level candidates
  m_ESAlgSettings->SetInitialValue();
  m_ETAlgSettings->SetInitialValue(); 
  m_CCAlgSettings->SetInitialValue();
  m_ACAlgSettings->SetInitialValue();  
  m_CLAlgSettings->SetInitialValue(); 
  m_CMAlgSettings->SetInitialValue();
  m_AMAlgSettings->SetInitialValue(); 
  m_CCAlgSettings->SetConeValue( PI/4., 50., PI/6., 30.);
  m_CMAlgSettings->SetUseTrk(true);
  //m_ESAlgSettings->SetDebug(2);
  //m_CMAlgSettings->SetDebug(2); 

  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol ); 
  m_etmatchingAlg     ->RunAlgorithm( *m_ETAlgSettings, dataCol );
  m_coneclusterAlg->    RunAlgorithm( *m_CCAlgSettings, dataCol);
  m_candidatemakingAlg->RunAlgorithm( *m_CMAlgSettings, dataCol );
  //m_clustermergingAlg ->RunAlgorithm( *m_CLAlgSettings, dataCol );

    dataCol.BlockVec_iter0       = dataCol.BlockVec;
    dataCol.LayerCol_iter0       = dataCol.LayerCol;
    dataCol.Shower2DCol_iter0    = dataCol.Shower2DCol;
    dataCol.GoodClus3DCol_iter0  = dataCol.GoodClus3DCol;
    dataCol.BadClus3DCol_iter0   = dataCol.BadClus3DCol;
    dataCol.Clus3DCol_iter0      = dataCol.Clus3DCol;

//dataCol.PrintLayer();
//dataCol.PrintShower();
//dataCol.Print3DClus(); 
//for(int ib=0; ib<dataCol.BlockVec.size(); ib++){ cout<<"  Block #"<<ib<<'\n';  dataCol.BlockVec[ib].PrintCandidates(); }

  //Iteration 
  dataCol.ClearLayer();
  dataCol.ClearShower();
  dataCol.ClearCluster();
  dataCol.ClearPFO();

  m_ESAlgSettings->SetValues();
  m_ESAlgSettings->SetUseCandidate(true); 
  m_ETAlgSettings->SetUseChi2(true);
  m_ACAlgSettings->SetRnodeThres(70, 2.);
  m_ACAlgSettings->SetKappaOrderWeight(2, 1, 0);
  m_AMAlgSettings->SetGoodTreeLevel(3, 3);
  m_AMAlgSettings->SetMergeTrees(true);
  m_CLAlgSettings->SetMergeGoodCluster(true);
  m_CLAlgSettings->SetMergeBadCluster(true);

  //m_ACAlgSettings->SetDebug(2);

  m_energysplittingAlg ->RunAlgorithm( *m_ESAlgSettings, dataCol );
  m_etmatchingAlg      ->RunAlgorithm( *m_ETAlgSettings, dataCol );
  m_arborclusteringAlg ->RunAlgorithm( *m_ACAlgSettings, dataCol );
  m_arbortreemergingAlg->RunAlgorithm( *m_AMAlgSettings, dataCol );
  m_clustermergingAlg  ->RunAlgorithm( *m_CLAlgSettings, dataCol );

  return StatusCode::SUCCESS;

}


StatusCode EcalClusterReconstruction::ClearAlgorithm(){

  delete m_energysplittingAlg;
  delete m_etmatchingAlg;
  delete m_coneclusterAlg;
  delete m_candidatemakingAlg;
  delete m_clustermergingAlg; 
  delete m_arborclusteringAlg;
  delete m_arbortreemergingAlg; 

  delete m_ESAlgSettings;
  delete m_ETAlgSettings;
  delete m_CCAlgSettings;
  delete m_CMAlgSettings;
  delete m_CLAlgSettings; 
  delete m_ACAlgSettings;
  delete m_AMAlgSettings;

  return StatusCode::SUCCESS;

}


