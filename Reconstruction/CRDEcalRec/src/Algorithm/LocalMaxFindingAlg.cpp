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

  std::vector<PandoraPlus::CaloBlock*> m_blocks = m_datacol.BlockCol;
  std::vector<PandoraPlus::CaloTower*> m_TowerCol; m_TowerCol.clear();

  int Nblocks = m_blocks.size();
  for(int ib=0;ib<Nblocks;ib++){
    GetLocalMax( m_blocks[ib] );

    int b_module = m_blocks[ib]->getModule();
    int b_stave  = m_blocks[ib]->getStave();
    int b_part   = m_blocks[ib]->getPart();

    PandoraPlus::CaloTower *m_tower = nullptr; 
    for(int it=0; it<m_TowerCol.size(); it++){
      if( m_TowerCol[it]->getModule() == b_module && 
          m_TowerCol[it]->getStave()  == b_stave  &&
          m_TowerCol[it]->getPart()   == b_part)  m_tower = m_TowerCol[it];
    }

    if(!m_tower){
      m_tower = new PandoraPlus::CaloTower();
      m_tower->setTowerID( m_blocks.at(ib)->getModule(),
                           m_blocks.at(ib)->getStave(),
                           m_blocks.at(ib)->getPart() );
      m_tower->addBlock( m_blocks[ib] );
      m_TowerCol.push_back(m_tower);
    }
    else{
      m_tower->addBlock( m_blocks[ib] );
    }

  }
//cout<<"  LocalMaxFinding: tower size = "<<m_TowerCol.size();
//cout<<"  #Blocks in tower: ";
//for(int i=0; i<m_TowerCol.size(); i++) cout<<m_TowerCol[i]->getBlocks().size()<<"  ";
//cout<<endl;

  m_datacol.TowerCol = m_TowerCol;
  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::GetLocalMax(PandoraPlus::CaloBlock* m_block){
  if(m_block->getBarXCol().size()==0 && m_block->getBarYCol().size()==0) return StatusCode::SUCCESS;

  std::vector<const PandoraPlus::CaloBar*> m_barXCol = m_block->getBarXCol();
  std::vector<const PandoraPlus::CaloBar*> m_barYCol = m_block->getBarYCol();

  std::vector<const PandoraPlus::CaloBar*> localMaxXCol; localMaxXCol.clear();
  std::vector<const PandoraPlus::CaloBar*> localMaxYCol; localMaxYCol.clear();
  GetLocalMaxBar( m_barXCol, localMaxXCol );
  GetLocalMaxBar( m_barYCol, localMaxYCol );

  std::vector<const PandoraPlus::CaloBarShower*> m_showerColX; m_showerColX.clear();
  std::vector<const PandoraPlus::CaloBarShower*> m_showerColY; m_showerColY.clear();
  for(int j=0; j<localMaxXCol.size(); j++){
    PandoraPlus::CaloBarShower* m_shower = new PandoraPlus::CaloBarShower();
    m_shower->addBar( localMaxXCol[j] );
    m_shower->setSeed( localMaxXCol[j] );
    m_shower->setIDInfo( localMaxXCol[j]->getModule(),
                        localMaxXCol[j]->getStave(),
                        localMaxXCol[j]->getDlayer(),
                        localMaxXCol[j]->getPart(),
                        localMaxXCol[j]->getSlayer() );
    m_showerColX.push_back(m_shower);
  }
  for(int j=0; j<localMaxYCol.size(); j++){
    PandoraPlus::CaloBarShower* m_shower = new PandoraPlus::CaloBarShower();
    m_shower->addBar( localMaxYCol[j] );
    m_shower->setSeed( localMaxYCol[j] );
    m_shower->setIDInfo( localMaxYCol[j]->getModule(),
                        localMaxYCol[j]->getStave(),
                        localMaxYCol[j]->getDlayer(),
                        localMaxYCol[j]->getPart(),
                        localMaxYCol[j]->getSlayer() );
    m_showerColY.push_back(m_shower);
  }

  m_block->setShowerXCol(m_showerColX);
  m_block->setShowerYCol(m_showerColY);

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::GetLocalMaxBar( std::vector<const PandoraPlus::CaloBar*>& barCol, std::vector<const PandoraPlus::CaloBar*>& localMaxCol ){
  std::sort( barCol.begin(), barCol.end(), compBar );

  for(int ib=0; ib<barCol.size(); ib++){
    std::vector<const PandoraPlus::CaloBar*> m_neighbors = getNeighbors( barCol[ib], barCol );
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


std::vector<const PandoraPlus::CaloBar*> LocalMaxFindingAlg::getNeighbors( const PandoraPlus::CaloBar* seed, std::vector<const PandoraPlus::CaloBar*>& barCol){
  std::vector<const PandoraPlus::CaloBar*> m_neighbor;
  for(int i=0;i<barCol.size();i++){
    if( seed->isNeighbor(barCol[i]) ) m_neighbor.push_back(barCol[i]);
  }
  if(m_neighbor.size()>2) std::cout<<"WARNING: more than 2 hits in neighborCol!!"<<std::endl;

  return m_neighbor;
}

#endif
