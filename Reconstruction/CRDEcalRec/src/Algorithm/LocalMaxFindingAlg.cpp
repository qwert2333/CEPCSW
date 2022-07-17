#ifndef _LOCALMAXFINDING_ALG_C
#define _LOCALMAXFINDING_ALG_C

#include "Algorithm/LocalMaxFindingAlg.h"
using namespace PandoraPlus; 

StatusCode LocalMaxFindingAlg::ReadSettings(PandoraPlus::Settings& m_settings){
  settings = m_settings;

  if(settings.map_floatPars.find("Eth_localMax")==settings.map_floatPars.end()) settings.map_floatPars["Eth_localMax"] = 0.005;
  if(settings.map_floatPars.find("Eth_MaxWithNeigh")==settings.map_floatPars.end()) settings.map_floatPars["Eth_MaxWithNeigh"] = 0.;
  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::Initialize(){

  return StatusCode::SUCCESS;
}


StatusCode LocalMaxFindingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol){

  std::vector<PandoraPlus::Calo2DCluster*> m_2dClusCol = m_datacol.Cluster2DCol;

cout<<"LocalMaxFinding: input 2DCluster size = "<<m_2dClusCol.size()<<endl;  
  for(int ic=0; ic<m_2dClusCol.size(); ic++){
    GetLocalMax( m_2dClusCol[ic] );

    for(int is=0; is<m_2dClusCol[ic]->getShowerUCol().size(); is++)
      m_datacol.bk_BarShowerCol.push_back( const_cast<PandoraPlus::CaloBarShower *>(m_2dClusCol[ic]->getShowerUCol()[is]) );
    for(int is=0; is<m_2dClusCol[ic]->getShowerVCol().size(); is++)
      m_datacol.bk_BarShowerCol.push_back( const_cast<PandoraPlus::CaloBarShower *>(m_2dClusCol[ic]->getShowerVCol()[is]) );
  }

//cout<<"  LocalMaxFinding: tower size = "<<m_TowerCol.size();
//cout<<"  #Blocks in tower: ";
//for(int i=0; i<m_TowerCol.size(); i++) cout<<m_TowerCol[i]->getBlocks().size()<<"  ";
//cout<<endl;

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::GetLocalMax(PandoraPlus::Calo2DCluster* m_2dClus){
  if(m_2dClus->getBarUCol().size()==0 && m_2dClus->getBarVCol().size()==0) return StatusCode::SUCCESS;

  std::vector<const PandoraPlus::CaloUnit*> m_barUCol = m_2dClus->getBarUCol();
  std::vector<const PandoraPlus::CaloUnit*> m_barVCol = m_2dClus->getBarVCol();

  std::vector<const PandoraPlus::CaloUnit*> localMaxUCol; localMaxUCol.clear();
  std::vector<const PandoraPlus::CaloUnit*> localMaxVCol; localMaxVCol.clear();
  GetLocalMaxBar( m_barUCol, localMaxUCol );
  GetLocalMaxBar( m_barVCol, localMaxVCol );

  std::vector<const PandoraPlus::CaloBarShower*> m_showerColU; m_showerColU.clear();
  std::vector<const PandoraPlus::CaloBarShower*> m_showerColV; m_showerColV.clear();
  for(int j=0; j<localMaxUCol.size(); j++){
    PandoraPlus::CaloBarShower* m_shower = new PandoraPlus::CaloBarShower();
    m_shower->addBar( localMaxUCol[j] );
    m_shower->setSeed( localMaxUCol[j] );
    m_shower->setIDInfo( localMaxUCol[j]->getModule(),
                        localMaxUCol[j]->getStave(),
                        localMaxUCol[j]->getDlayer(),
                        localMaxUCol[j]->getPart(),
                        localMaxUCol[j]->getSlayer() );
    m_showerColU.push_back(m_shower);
  }
  for(int j=0; j<localMaxVCol.size(); j++){
    PandoraPlus::CaloBarShower* m_shower = new PandoraPlus::CaloBarShower();
    m_shower->addBar( localMaxVCol[j] );
    m_shower->setSeed( localMaxVCol[j] );
    m_shower->setIDInfo( localMaxVCol[j]->getModule(),
                        localMaxVCol[j]->getStave(),
                        localMaxVCol[j]->getDlayer(),
                        localMaxVCol[j]->getPart(),
                        localMaxVCol[j]->getSlayer() );
    m_showerColV.push_back(m_shower);
  }

  m_2dClus->setShowerUCol(m_showerColU);
  m_2dClus->setShowerVCol(m_showerColV);

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::GetLocalMaxBar( std::vector<const PandoraPlus::CaloUnit*>& barCol, std::vector<const PandoraPlus::CaloUnit*>& localMaxCol ){
  //std::sort( barCol.begin(), barCol.end(), compBar );

  for(int ib=0; ib<barCol.size(); ib++){
    std::vector<const PandoraPlus::CaloUnit*> m_neighbors = getNeighbors( barCol[ib], barCol );
    if( m_neighbors.size()==0 && barCol[ib]->getEnergy()>settings.map_floatPars["Eth_localMax"] ) { localMaxCol.push_back( barCol[ib] ); continue; }

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
