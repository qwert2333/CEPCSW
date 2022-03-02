#include "EcalClusterReconstruction.h"

void EcalClusterReconstruction::Initialize(){

  //Initialize settings
  m_LFAlgSettings  = new LocalMaxFindingAlg::Settings();
  m_ESAlgSettings  = new EnergySplittingAlg::Settings();
  m_CC2AlgSettings = new ConeClustering2DAlg::Settings();
  m_HCAlgSettings  = new HoughClusteringAlg::Settings();
  m_CMAlgSettings  = new ShadowMakingAlg::Settings(); 
  m_ETAlgSettings  = new EnergyTimeMatchingAlg::Settings();
  m_CCAlgSettings  = new ConeClusteringAlg::Settings();
  m_CLAlgSettings  = new ClusterMergingAlg::Settings();
  m_CIDAlgSettings = new BasicClusterIDAlg::Settings();


  //Initialize Algorthm
  m_localmaxfindingAlg  = new LocalMaxFindingAlg();
  m_energysplittingAlg  = new EnergySplittingAlg();
  m_coneclus2DAlg       = new ConeClustering2DAlg();
  m_houghclusAlg        = new HoughClusteringAlg();
  m_shadowmakingAlg     = new ShadowMakingAlg(); 
  m_etmatchingAlg       = new EnergyTimeMatchingAlg();
  m_coneclusterAlg      = new ConeClusteringAlg();
  m_clustermergingAlg   = new ClusterMergingAlg();
  m_clusteridAlg        = new BasicClusterIDAlg();  

}

StatusCode EcalClusterReconstruction::RunAlgorithm( EcalClusterReconstruction::Settings& settings,  PandoraPlusDataCol& dataCol ){


  //Iteration 0: pre-treatment, get initial level candidates
  m_LFAlgSettings->SetInitialValue();
  m_ESAlgSettings->SetInitialValue();
  m_CC2AlgSettings->SetInitialValue();
  m_HCAlgSettings->SetInitialValue();
  m_CMAlgSettings->SetInitialValue(); 
  m_ETAlgSettings->SetInitialValue();
  m_CCAlgSettings->SetInitialValue();
  m_CLAlgSettings->SetInitialValue();
  m_CIDAlgSettings->SetInitialValue();


  //Stage 1
  cout<<"Finding Local Maximum"<<endl;
  m_localmaxfindingAlg->RunAlgorithm( *m_LFAlgSettings, dataCol );
  cout<<"Hough Clustering"<<endl;
  m_houghclusAlg      ->RunAlgorithm( *m_HCAlgSettings, dataCol );

//  cout<<"Energy splitting"<<endl;
//  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, dataCol );
//  cout<<"Chi2 matching"<<endl;
//  m_etmatchingAlg     ->RunAlgorithm( *m_ETAlgSettings, dataCol );


/*
cout<<"Block Info: "<<endl;
for(int ib=0; ib<dataCol.BlockVec.size(); ib++){
  printf("  #%d Block: ID=(%d, %d, %d, %d). Bar shower size: X=%d, Y=%d \n", ib, dataCol.BlockVec[ib].getModule(), dataCol.BlockVec[ib].getStave(), dataCol.BlockVec[ib].getPart(), dataCol.BlockVec[ib].getDlayer(), dataCol.BlockVec[ib].getShowerXCol().size(), dataCol.BlockVec[ib].getShowerYCol().size() );

  cout<<"  Print BarShowersX "<<endl;
  for(int is=0 ;is<dataCol.BlockVec[ib].getShowerXCol().size(); is++) 
  printf("    #%d showerX: cellID=(%d, %d, %d, %d) \n", is, dataCol.BlockVec[ib].getShowerXCol()[is].getModule(), dataCol.BlockVec[ib].getShowerXCol()[is].getStave(), dataCol.BlockVec[ib].getShowerXCol()[is].getPart(), dataCol.BlockVec[ib].getShowerXCol()[is].getDlayer() );

  cout<<"  Print BarShowersY "<<endl;
  for(int is=0 ;is<dataCol.BlockVec[ib].getShowerYCol().size(); is++)
  printf("    #%d showerX: cellID=(%d, %d, %d, %d) \n", is, dataCol.BlockVec[ib].getShowerYCol()[is].getModule(), dataCol.BlockVec[ib].getShowerYCol()[is].getStave(), dataCol.BlockVec[ib].getShowerYCol()[is].getPart(), dataCol.BlockVec[ib].getShowerYCol()[is].getDlayer() );
}

cout<<endl;
cout<<endl;

cout<<"Tower Info: "<<endl;
for(int it=0; it<dataCol.TowerCol.size(); it++){
  printf("  #%d Tower: ID=(%d, %d, %d). Block size: %d \n", dataCol.TowerCol[it].getModule(), dataCol.TowerCol[it].getStave(), dataCol.TowerCol[it].getPart(), dataCol.TowerCol[it].getBlocks().size() );

for(int ib=0; ib<dataCol.TowerCol[it].getBlocks().size(); ib++){
  printf("  #%d Block: ID=(%d, %d, %d, %d). Bar shower size: X=%d, Y=%d \n", ib, dataCol.TowerCol[it].getBlocks()[ib].getModule(), dataCol.TowerCol[it].getBlocks()[ib].getStave(), dataCol.TowerCol[it].getBlocks()[ib].getPart(), dataCol.TowerCol[it].getBlocks()[ib].getDlayer(), dataCol.TowerCol[it].getBlocks()[ib].getShowerXCol().size(), dataCol.TowerCol[it].getBlocks()[ib].getShowerYCol().size() );

  cout<<"  Print BarShowersX "<<endl;
  for(int is=0 ;is<dataCol.TowerCol[it].getBlocks()[ib].getShowerXCol().size(); is++)
  printf("    #%d showerX: cellID=(%d, %d, %d, %d) \n", is, dataCol.TowerCol[it].getBlocks()[ib].getShowerXCol()[is].getModule(), dataCol.TowerCol[it].getBlocks()[ib].getShowerXCol()[is].getStave(), dataCol.TowerCol[it].getBlocks()[ib].getShowerXCol()[is].getPart(), dataCol.TowerCol[it].getBlocks()[ib].getShowerXCol()[is].getDlayer() );

  cout<<"  Print BarShowersY "<<endl;
  for(int is=0 ;is<dataCol.TowerCol[it].getBlocks()[ib].getShowerYCol().size(); is++)
  printf("    #%d showerX: cellID=(%d, %d, %d, %d) \n", is, dataCol.TowerCol[it].getBlocks()[ib].getShowerYCol()[is].getModule(), dataCol.TowerCol[it].getBlocks()[ib].getShowerYCol()[is].getStave(), dataCol.TowerCol[it].getBlocks()[ib].getShowerYCol()[is].getPart(), dataCol.TowerCol[it].getBlocks()[ib].getShowerYCol()[is].getDlayer() );
}

}
*/

  dataCol.PrintShower();
  dataCol.Print3DClus();

  return StatusCode::SUCCESS;

}


StatusCode EcalClusterReconstruction::ClearAlgorithm(){

  delete m_energysplittingAlg;
  delete m_coneclus2DAlg; 
  delete m_houghclusAlg;
  delete m_shadowmakingAlg; 
  delete m_etmatchingAlg;
  delete m_coneclusterAlg;
  delete m_clustermergingAlg;
  delete m_clusteridAlg;

  delete m_ESAlgSettings;
  delete m_CC2AlgSettings; 
  delete m_HCAlgSettings;
  delete m_CMAlgSettings; 
  delete m_ETAlgSettings;
  delete m_CCAlgSettings;
  delete m_CLAlgSettings;
  delete m_CIDAlgSettings;


  return StatusCode::SUCCESS;

}


