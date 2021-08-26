#include "EcalClusterReconstruction.h"

void EcalClusterReconstruction::Initialize(){

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
}

StatusCode EcalClusterReconstruction::RunAlgorithm( EcalClusterReconstruction::Settings& settings,  PandoraPlusDataCol& dataCol ){


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
  //m_ESAlgSettings->SetDebug(2);
  //m_CMAlgSettings->SetDebug(2); 

  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol ); 
  m_etmatchingAlg     ->RunAlgorithm( *m_ETAlgSettings, dataCol );
  m_coneclusterAlg    ->RunAlgorithm( *m_CCAlgSettings, dataCol);
  m_shadowmakingAlg   ->RunAlgorithm( *m_CMAlgSettings, dataCol );
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
//std::cout<<"Start iteration"<<std::endl;

  dataCol.ClearLayer();
  dataCol.ClearShower();
  dataCol.ClearCluster();
  dataCol.ClearPFO();

  m_ESAlgSettings->SetValues();
  m_ESAlgSettings->SetUseCandidate(true);
  m_ETAlgSettings->SetUseChi2(true);
  m_CCAlgSettings->SetGoodClusLevel(5);

  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol );
  m_etmatchingAlg->     RunAlgorithm( *m_ETAlgSettings, dataCol );
  m_coneclusterAlg->    RunAlgorithm( *m_CCAlgSettings, dataCol);
//  m_arborclusteringAlg ->RunAlgorithm( *m_ACAlgSettings, dataCol );
//  m_arbortreemergingAlg->RunAlgorithm( *m_AMAlgSettings, dataCol );
//dataCol.PrintShower();


std::cout<<"   Before cluster ID and depart: shower col size = "<<dataCol.Shower2DCol.size()<<std::endl;
std::cout<<"                                 good cluster size = "<<dataCol.GoodClus3DCol.size();
std::cout<<", bad cluster size = "<<dataCol.BadClus3DCol.size();
std::cout<<", total cluster size = "<<dataCol.Clus3DCol.size()<<std::endl;
std::cout<<"  Good cluster stdDevE: "; 
for(int ig=0 ;ig<dataCol.GoodClus3DCol.size(); ig++) std::cout<<dataCol.GoodClus3DCol[ig].getStdDevE()<<'\t';
std::cout<<std::endl;


  m_CIDAlgSettings->SetDepartShower(true);
  m_clusteridAlg->RunAlgorithm( *m_CIDAlgSettings, dataCol );

std::cout<<"    MIP shower col size: "<<dataCol.MIPShower2DCol.size()<<std::endl;
std::cout<<"    EM shower col size: "<<dataCol.EMShower2DCol.size()<<std::endl;
std::cout<<"    other shower col size: "<<dataCol.Shower2DCol.size()<<std::endl;
//dataCol.Print3DClus(); 


  m_CCAlgSettings->SetClusterType("MIP");
  m_CCAlgSettings->SetConeValue(PI/4., 50., PI/6., 30.);
  m_CCAlgSettings->SetGoodClusLevel(5);
  m_CCAlgSettings->SetUseCandidate(1);
  m_CCAlgSettings->SetOverwrite(true);
  m_coneclusterAlg->RunAlgorithm( *m_CCAlgSettings, dataCol);
std::cout<<"  After clustering MIP: good cluster size = "<<dataCol.GoodClus3DCol.size();
std::cout<<", bad cluster size = "<<dataCol.BadClus3DCol.size();
std::cout<<", total cluster size = "<<dataCol.Clus3DCol.size()<<std::endl;


  m_CCAlgSettings->SetClusterType("EM");
  m_CCAlgSettings->SetUseCandidate(1);
  m_CCAlgSettings->SetGoodClusLevel(4);
  m_CCAlgSettings->SetOverwrite(false);
  m_coneclusterAlg->RunAlgorithm( *m_CCAlgSettings, dataCol);  
std::cout<<"  After clustering EM: good cluster size = "<<dataCol.GoodClus3DCol.size();
std::cout<<", bad cluster size = "<<dataCol.BadClus3DCol.size();
std::cout<<", total cluster size = "<<dataCol.Clus3DCol.size()<<std::endl;


  m_ACAlgSettings->SetClusterType("");
  m_ACAlgSettings->SetRthStart(2, 30.);
  m_ACAlgSettings->SetRthPars(60, 0.);
  m_ACAlgSettings->SetKappaOrderWeight(2, 1, 0);
  m_AMAlgSettings->SetClusterType("");
  m_AMAlgSettings->SetGoodTreeLevel(5, 4, 10);
  m_AMAlgSettings->SetOverwrite(false);
  m_arborclusteringAlg ->RunAlgorithm( *m_ACAlgSettings,  dataCol );
  m_arbortreemergingAlg->RunAlgorithm( *m_AMAlgSettings,  dataCol );
std::cout<<"  After clustering other: good cluster size = "<<dataCol.GoodClus3DCol.size();
std::cout<<", bad cluster size = "<<dataCol.BadClus3DCol.size();
std::cout<<", total cluster size = "<<dataCol.Clus3DCol.size()<<std::endl;

//dataCol.Print3DClus(); 

  m_CLAlgSettings->SetMergeGoodCluster(true);
  m_CLAlgSettings->SetMergeBadCluster(true);
  m_CLAlgSettings->SetMergeEMTail(false);
  m_clustermergingAlg  ->RunAlgorithm( *m_CLAlgSettings, dataCol );

  //Cluster ID
  //Shadow cluster making
  
  //Next iteration: 
  //Cluster splitting
  //E-T matching (chi2 and shadow cluster)
  //clustering (with shadow cluster)
  //Identify MIP & EM
  //mask MIP and EM
  //clustering others again
  //cluster merging (hadrons only)
  //End. 


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


