#ifndef _LOCALMAXFINDING_ALG_C
#define _LOCALMAXFINDING_ALG_C

#include "Algorithm/LocalMaxFindingAlg.h"
using namespace PandoraPlus; 

StatusCode LocalMaxFindingAlg::ReadSettings(PandoraPlus::Settings& m_settings){
  settings = m_settings;

  if(settings.map_floatPars.find("Eth_localMax")==settings.map_floatPars.end()) settings.map_floatPars["Eth_localMax"] = 0.005;
  if(settings.map_floatPars.find("Eth_MaxWithNeigh")==settings.map_floatPars.end()) settings.map_floatPars["Eth_MaxWithNeigh"] = 0.;
  if(settings.map_stringPars.find("OutputLocalMaxName")==settings.map_stringPars.end()) settings.map_stringPars["OutputLocalMaxName"] = "AllLocalMax";
  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::Initialize(){

  return StatusCode::SUCCESS;
}


StatusCode LocalMaxFindingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol){
/*
  std::vector<PandoraPlus::Calo2DCluster*> m_2dClusCol = m_datacol.Cluster2DCol;

cout<<"LocalMaxFinding: input 2DCluster size = "<<m_2dClusCol.size()<<endl;  
  for(int ic=0; ic<m_2dClusCol.size(); ic++){
    GetLocalMax( m_2dClusCol[ic] );

    for(int is=0; is<m_2dClusCol[ic]->getShowerUCol().size(); is++)
      m_datacol.bk_BarShowerCol.push_back( const_cast<PandoraPlus::Calo1DCluster *>(m_2dClusCol[ic]->getShowerUCol()[is]) );
    for(int is=0; is<m_2dClusCol[ic]->getShowerVCol().size(); is++)
      m_datacol.bk_BarShowerCol.push_back( const_cast<PandoraPlus::Calo1DCluster *>(m_2dClusCol[ic]->getShowerVCol()[is]) );
  }
*/
  std::vector<PandoraPlus::Calo3DCluster*>* p_3DClusters = &(m_datacol.Cluster3DCol);
  if(!p_3DClusters) {std::cout<<"ERROR: No 3DCluster in present data collection! "<<std::endl; return StatusCode::FAILURE; }

  for(int i3d=0; i3d<p_3DClusters->size(); i3d++){
    std::vector<const PandoraPlus::Calo1DCluster*> m_localMaxUCol; m_localMaxUCol.clear();
    std::vector<const PandoraPlus::Calo1DCluster*> m_localMaxVCol; m_localMaxVCol.clear();

    std::vector<const Calo2DCluster*> m_2dClusCol = p_3DClusters->at(i3d)->getCluster();
    for(int i2d=0; i2d<m_2dClusCol.size(); i2d++) GetLocalMax(m_2dClusCol[i2d], m_localMaxUCol, m_localMaxVCol);

    p_3DClusters->at(i3d)->setLocalMax(settings.map_stringPars["OutputLocalMaxName"], m_localMaxUCol, settings.map_stringPars["OutputLocalMaxName"], m_localMaxVCol);

//cout<<"LocalMaxFinding: LocalMax in 3DCluster #"<<i3d<<endl;
//cout<<"LocalMaxU: "<<endl;
//for(int i=0; i<m_localMaxUCol.size(); i++)
//  printf("  #%d: (%.3f, %.3f, %.3f) \n ", i, m_localMaxUCol[i]->getPos().x(), m_localMaxUCol[i]->getPos().y(), m_localMaxUCol[i]->getPos().z());
//cout<<"LocalMaxV: "<<endl;
//for(int i=0; i<m_localMaxVCol.size(); i++)
//  printf("  #%d: (%.3f, %.3f, %.3f) \n ", i, m_localMaxVCol[i]->getPos().x(), m_localMaxVCol[i]->getPos().y(), m_localMaxVCol[i]->getPos().z());
    m_datacol.bk_Cluster1DCol.insert(m_datacol.bk_Cluster1DCol.end(), m_localMaxUCol.begin(), m_localMaxUCol.end());
    m_datacol.bk_Cluster1DCol.insert(m_datacol.bk_Cluster1DCol.end(), m_localMaxVCol.begin(), m_localMaxVCol.end());
  }
  p_3DClusters = nullptr;

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::GetLocalMax( const PandoraPlus::Calo2DCluster* m_2dClus, 
                                            std::vector<const PandoraPlus::Calo1DCluster*>& m_outputU, 
                                            std::vector<const PandoraPlus::Calo1DCluster*>& m_outputV){

  if(m_2dClus->getBarUCol().size()==0 && m_2dClus->getBarVCol().size()==0) return StatusCode::SUCCESS;

  std::vector<const PandoraPlus::CaloUnit*> m_barUCol = m_2dClus->getBarUCol();
  std::vector<const PandoraPlus::CaloUnit*> m_barVCol = m_2dClus->getBarVCol();

//cout<<"  LocalMaxFindingAlg::GetLocalMax: Input bar collection size: "<<m_barUCol.size()<<"  "<<m_barVCol.size()<<endl;

  std::vector<const PandoraPlus::CaloUnit*> localMaxUCol; localMaxUCol.clear();
  std::vector<const PandoraPlus::CaloUnit*> localMaxVCol; localMaxVCol.clear();
  GetLocalMaxBar( m_barUCol, localMaxUCol );
  GetLocalMaxBar( m_barVCol, localMaxVCol );

//cout<<"  LocalMaxFindingAlg::GetLocalMax: Found local max bar size: "<<localMaxUCol.size()<<"  "<<localMaxVCol.size()<<endl;
//cout<<"  Transfer bar to barShower"<<endl;

  for(int j=0; j<localMaxUCol.size(); j++){
    PandoraPlus::Calo1DCluster* m_shower = new PandoraPlus::Calo1DCluster();
    m_shower->addUnit( localMaxUCol[j] );
    m_shower->addSeed( localMaxUCol[j] );
    m_outputU.push_back(m_shower);
  }
  for(int j=0; j<localMaxVCol.size(); j++){
    PandoraPlus::Calo1DCluster* m_shower = new PandoraPlus::Calo1DCluster();
    m_shower->addUnit( localMaxVCol[j] );
    m_shower->addSeed( localMaxVCol[j] );
    m_outputV.push_back(m_shower);
  }
//cout<<"  Output bar shower size: "<<m_outputU.size()<<"  "<<m_outputV.size()<<endl;
//cout<<endl;

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::GetLocalMaxBar( std::vector<const PandoraPlus::CaloUnit*>& barCol, std::vector<const PandoraPlus::CaloUnit*>& localMaxCol ){
  //std::sort( barCol.begin(), barCol.end(), compBar );

  for(int ib=0; ib<barCol.size(); ib++){
    std::vector<const PandoraPlus::CaloUnit*> m_neighbors = getNeighbors( barCol[ib], barCol );
    if( m_neighbors.size()==0 && barCol[ib]->getEnergy()>settings.map_floatPars["Eth_localMax"] ) { 
      localMaxCol.push_back( barCol[ib] ); continue; 
    }

    bool isLocalMax=true;
    double Eneigh=0;
    if(barCol[ib]->getEnergy()<settings.map_floatPars["Eth_localMax"]) isLocalMax = false;
    for(int j=0;j<m_neighbors.size();j++){
      if(m_neighbors[j]->getEnergy()>barCol[ib]->getEnergy()) isLocalMax=false;
      Eneigh += m_neighbors[j]->getEnergy();
    }
    if( (barCol[ib]->getEnergy()/(barCol[ib]->getEnergy()+Eneigh))<settings.map_floatPars["Eth_MaxWithNeigh"] ) isLocalMax = false;

    if(isLocalMax) localMaxCol.push_back( barCol[ib] );
  }

  return StatusCode::SUCCESS;
}


std::vector<const PandoraPlus::CaloUnit*> LocalMaxFindingAlg::getNeighbors( const PandoraPlus::CaloUnit* seed, std::vector<const PandoraPlus::CaloUnit*>& barCol){
  std::vector<const PandoraPlus::CaloUnit*> m_neighbor; m_neighbor.clear();
  for(int i=0;i<barCol.size();i++){
    bool fl_neighbor = false; 
    if( seed->getModule()==barCol[i]->getModule() && 
        seed->getPart()==barCol[i]->getPart() && 
        seed->getStave()==barCol[i]->getStave() && 
        seed->getDlayer()==barCol[i]->getDlayer() &&
        seed->getSlayer()==barCol[i]->getSlayer() &&
        abs( seed->getBar()-barCol[i]->getBar() )==1 ) fl_neighbor=true;
    else if( seed->getModule()==barCol[i]->getModule() && 
             seed->getStave()==barCol[i]->getStave() && 
             ( ( seed->getPart()-barCol[i]->getPart()==1 && seed->isAtLowerEdgePhi() && barCol[i]->isAtUpperEdgePhi() ) ||
               ( barCol[i]->getPart()-seed->getPart()==1 && seed->isAtUpperEdgePhi() && barCol[i]->isAtLowerEdgePhi() ) ) ) fl_neighbor=true;
    else if( seed->getModule()==barCol[i]->getModule() && 
             seed->getPart()==barCol[i]->getPart() &&
             ( ( seed->getStave()-barCol[i]->getStave()==1 && seed->isAtLowerEdgeZ() && barCol[i]->isAtUpperEdgeZ() ) || 
               ( barCol[i]->getStave()-seed->getStave()==1 && seed->isAtUpperEdgeZ() && barCol[i]->isAtLowerEdgeZ() ) ) ) fl_neighbor=true;


    if(fl_neighbor) m_neighbor.push_back( barCol[i] );

  }
  if(m_neighbor.size()>2) std::cout<<"WARNING: more than 2 hits in neighborCol!!"<<std::endl;

  return m_neighbor;
}

#endif
