#include "EcalClusterReconstruction.h"

StatusCode EcalClusterReconstruction::Initialize(){

  //Initialize settings
  m_ESAlgSettings  = new EnergySplittingAlg::Settings();
  m_ETAlgSettings  = new EnergyTimeMatchingAlg::Settings();
  m_CCAlgSettings  = new ConeClusteringAlg::Settings();
  m_CMAlgSettings  = new ShadowMakingAlg::Settings(); 
  m_CLAlgSettings  = new ClusterMergingAlg::Settings(); 
  m_ACAlgSettings  = new ArborClusteringAlg::Settings();
  m_AMAlgSettings  = new ArborTreeMergingAlg::Settings(); 
  m_CIDAlgSettings = new BasicClusterIDAlg::Settings();

  //Initialize Algorthm
  m_energysplittingAlg  = new EnergySplittingAlg();
  m_etmatchingAlg       = new EnergyTimeMatchingAlg();
  m_coneclusterAlg      = new ConeClusteringAlg();
  m_shadowmakingAlg     = new ShadowMakingAlg();
  m_clustermergingAlg   = new ClusterMergingAlg();  
  m_arborclusteringAlg  = new ArborClusteringAlg();
  m_arbortreemergingAlg = new ArborTreeMergingAlg(); 
  m_clusteridAlg        = new BasicClusterIDAlg(); 

  return StatusCode::SUCCESS;
}

StatusCode EcalClusterReconstruction::RunAlgorithm( EcalClusterReconstruction::Settings& settings,  PandoraPlusDataCol& dataCol ){


std::cout<<"Start preparation"<<std::endl;
  //Iteration 0: pre-treatment, get initial level candidates
  m_ESAlgSettings->SetInitialValue();
  m_ETAlgSettings->SetInitialValue(); 
  m_CCAlgSettings->SetInitialValue();
  m_ACAlgSettings->SetInitialValue();  
  m_CLAlgSettings->SetInitialValue(); 
  m_CMAlgSettings->SetInitialValue();
  m_AMAlgSettings->SetInitialValue(); 
  m_CCAlgSettings->SetConeValue( PI/4., 50., PI/6., 30.);
  m_CCAlgSettings->SetGoodClusLevel(2); 
  m_CMAlgSettings->SetUseTrk(true);
  m_CIDAlgSettings->SetInitialValue(); 
  //m_CMAlgSettings->SetDebug(2); 

  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol ); 
  m_etmatchingAlg     ->RunAlgorithm( *m_ETAlgSettings, dataCol );
  m_coneclusterAlg    ->RunAlgorithm( *m_CCAlgSettings, dataCol);
  m_shadowmakingAlg   ->RunAlgorithm( *m_CMAlgSettings, dataCol );
  //m_clustermergingAlg ->RunAlgorithm( *m_CLAlgSettings, dataCol );

  //Iteration 
  std::cout<<"Start iteration"<<std::endl;

  dataCol.ClearLayer();
  dataCol.ClearShower();
  dataCol.ClearCluster();
  dataCol.ClearPFO();

  m_ESAlgSettings->SetValues();
  m_ESAlgSettings->SetUseCandidate(true);
  m_ETAlgSettings->SetUseChi2(true);
  m_ETAlgSettings->SetUseCandidate(0);
  m_CCAlgSettings->SetConeValue( PI/5., 40., PI/6., 30.);
  m_CCAlgSettings->SetGoodClusLevel(4);
  //m_ESAlgSettings->SetDebug(3);

  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol );
  m_etmatchingAlg->     RunAlgorithm( *m_ETAlgSettings, dataCol );
  m_coneclusterAlg->    RunAlgorithm( *m_CCAlgSettings, dataCol);
//  m_arborclusteringAlg ->RunAlgorithm( *m_ACAlgSettings, dataCol );
//  m_arbortreemergingAlg->RunAlgorithm( *m_AMAlgSettings, dataCol );
//dataCol.PrintShower();
//dataCol.Print3DClus(); 
//dataCol.PrintLayer(); 

std::cout<<"   Before cluster ID and depart: shower col size = "<<dataCol.Shower2DCol.size()<<std::endl;
std::cout<<"                                 good cluster size = "<<dataCol.GoodClus3DCol.size();
std::cout<<", bad cluster size = "<<dataCol.BadClus3DCol.size();
std::cout<<", total cluster size = "<<dataCol.Clus3DCol.size()<<std::endl;


std::cout<<"ECalRec: Merge EM tails"<<std::endl;
//  m_CLAlgSettings->SetMergeGoodCluster(false);
//  m_CLAlgSettings->SetMergeBadCluster(false);
//  m_CLAlgSettings->SetMergeEMTailShower(true);
//  m_clustermergingAlg  ->RunAlgorithm( *m_CLAlgSettings, dataCol );
//dataCol.PrintShower();
//dataCol.Print3DClus();


std::cout<<"ECalRec: clusterID 1st time"<<std::endl;

  m_CIDAlgSettings->SetDepartShower(true);
  m_clusteridAlg->RunAlgorithm( *m_CIDAlgSettings, dataCol );


std::cout<<"    MIP shower col size: "<<dataCol.MIPShower2DCol.size()<<std::endl;
std::cout<<"    EM shower col size: "<<dataCol.EMShower2DCol.size()<<std::endl;
std::cout<<"    other shower col size: "<<dataCol.Shower2DCol.size()<<std::endl;
//dataCol.PrintShower();
//dataCol.Print3DClus(); 


std::cout<<"ECalRec: MIP clustering"<<std::endl;
  m_CCAlgSettings->SetClusterType("MIP");
  m_CCAlgSettings->SetConeValue(PI/4., 50., PI/6., 30.);
  m_CCAlgSettings->SetGoodClusLevel(5);
  m_CCAlgSettings->SetUseCandidate(1);
  m_CCAlgSettings->SetOverwrite(true);
  m_coneclusterAlg->RunAlgorithm( *m_CCAlgSettings, dataCol);

std::cout<<"ECalRec: EM clustering"<<std::endl;
  m_CCAlgSettings->SetClusterType("EM");
  m_CCAlgSettings->SetUseCandidate(1);
  m_CCAlgSettings->SetGoodClusLevel(4);
  m_CCAlgSettings->SetOverwrite(false);
  m_coneclusterAlg->RunAlgorithm( *m_CCAlgSettings, dataCol);

std::cout<<"ECalRec: Others clustering"<<std::endl;
  m_ACAlgSettings->SetClusterType("");
  m_ACAlgSettings->SetRthStart(2, 30.);
  m_ACAlgSettings->SetRthPars(60, 0.);
  m_ACAlgSettings->SetKappaOrderWeight(2, 1, 0);
  m_AMAlgSettings->SetClusterType("");
  m_AMAlgSettings->SetGoodTreeLevel(5, 4, 12);
  m_AMAlgSettings->SetMergeTrees(true);
  m_AMAlgSettings->SetMergeTreePars(40, PI/10.);
  m_AMAlgSettings->SetOverwrite(false);
  m_arborclusteringAlg ->RunAlgorithm( *m_ACAlgSettings,  dataCol );
  m_arbortreemergingAlg->RunAlgorithm( *m_AMAlgSettings,  dataCol );
std::cout<<"  After clustering other: good cluster size = "<<dataCol.GoodClus3DCol.size();
std::cout<<", bad cluster size = "<<dataCol.BadClus3DCol.size();
std::cout<<", total cluster size = "<<dataCol.Clus3DCol.size()<<std::endl;

//dataCol.Print3DClus(); 

std::cout<<"EcalRec: Last Cluster merging"<<std::endl;
  m_CLAlgSettings->SetMergeGoodCluster(true);
  m_CLAlgSettings->SetMergeBadCluster(true);
  m_CLAlgSettings->SetMergeEMTailCluster(false);
  m_clustermergingAlg  ->RunAlgorithm( *m_CLAlgSettings, dataCol );


  return StatusCode::SUCCESS;

}


StatusCode EcalClusterReconstruction::ClearAlgorithm(){

  delete m_energysplittingAlg;
  delete m_etmatchingAlg;
  delete m_coneclusterAlg;
  delete m_shadowmakingAlg;
  delete m_clustermergingAlg; 
  delete m_arborclusteringAlg;
  delete m_arbortreemergingAlg; 
  delete m_clusteridAlg; 

  delete m_ESAlgSettings;
  delete m_ETAlgSettings;
  delete m_CCAlgSettings;
  delete m_CMAlgSettings;
  delete m_CLAlgSettings; 
  delete m_ACAlgSettings;
  delete m_AMAlgSettings;
  delete m_CIDAlgSettings; 

  return StatusCode::SUCCESS;

}


