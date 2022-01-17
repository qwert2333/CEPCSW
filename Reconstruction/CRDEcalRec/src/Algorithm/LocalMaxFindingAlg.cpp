#ifndef _FINDINGLOCALMAX_ALG_C
#define _FINDINGLOCALMAX_ALG_C

#include "Algorithm/LocalMaxFindingAlg.h"
void LocalMaxFindingAlg::Settings::SetInitialValue(){


}

StatusCode LocalMaxFindingAlg::Initialize(){

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::RunAlgorithm( LocalMaxFindingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){
  settings = m_settings; 

  std::vector<CRDEcalEDM::CRDCaloTower> m_TowerCol; m_TowerCol.clear();
  int Nblocks = m_datacol.BlockVec.size();
  for(int ib=0;ib<Nblocks;ib++){
    GetLocalMax( m_datacol.BlockVec[ib] );

    CRDEcalEDM::CRDCaloTower m_tower; m_tower.Clear();
    m_tower.SetTowerID( m_datacol.BlockVec[ib].getModule(),
                        m_datacol.BlockVec[ib].getStave(),
                        m_datacol.BlockVec[ib].getPart() );
    std::vector<CRDEcalEDM::CRDCaloTower>::iterator iter = find(m_TowerCol.begin(), m_TowerCol.end(), m_tower);
    if(iter==m_TowerCol.end()){ m_tower.AddBlock(m_datacol.BlockVec[ib]); m_TowerCol.push_back(m_tower); }
    else iter->AddBlock(m_datacol.BlockVec[ib]);
  }
  m_datacol.TowerCol = m_TowerCol;

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::GetLocalMax(CRDEcalEDM::CRDCaloBlock& m_block){
  if(m_block.getBarXCol().size()==0 && m_block.getBarYCol().size()==0) return StatusCode::SUCCESS;

  std::vector<CRDEcalEDM::CRDCaloBar> m_barXCol = m_block.getBarXCol();
  std::vector<CRDEcalEDM::CRDCaloBar> m_barYCol = m_block.getBarYCol();

  std::vector<CRDEcalEDM::CRDCaloBar> localMaxXCol; localMaxXCol.clear();
  std::vector<CRDEcalEDM::CRDCaloBar> localMaxYCol; localMaxYCol.clear();
  GetLocalMaxBar( m_barXCol, localMaxXCol );
  GetLocalMaxBar( m_barYCol, localMaxYCol );

  std::vector<CRDEcalEDM::CRDCaloBarShower> m_showerColX; m_showerColX.clear();
  std::vector<CRDEcalEDM::CRDCaloBarShower> m_showerColY; m_showerColY.clear();
  for(int j=0; j<localMaxXCol.size(); j++){
    CRDEcalEDM::CRDCaloBarShower m_shower; m_shower.Clear();
    m_shower.addBar( localMaxXCol[j] );
    m_shower.setSeed( localMaxXCol[j] );
    m_shower.setIDInfo( localMaxXCol[j].getModule(),
                        localMaxXCol[j].getStave(),
                        localMaxXCol[j].getDlayer(),
                        localMaxXCol[j].getPart(),
                        localMaxXCol[j].getSlayer() );
    m_showerColX.push_back(m_shower);
  }
  for(int j=0; j<localMaxYCol.size(); j++){
    CRDEcalEDM::CRDCaloBarShower m_shower; m_shower.Clear();
    m_shower.addBar( localMaxYCol[j] );
    m_shower.setSeed( localMaxYCol[j] );
    m_shower.setIDInfo( localMaxYCol[j].getModule(),
                        localMaxYCol[j].getStave(),
                        localMaxYCol[j].getDlayer(),
                        localMaxYCol[j].getPart(),
                        localMaxYCol[j].getSlayer() );
    m_showerColY.push_back(m_shower);
  }

  m_block.setShowerXCol(m_showerColX);
  m_block.setShowerYCol(m_showerColY);

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::GetLocalMaxBar( std::vector<CRDEcalEDM::CRDCaloBar>& barCol, std::vector<CRDEcalEDM::CRDCaloBar>& localMaxCol ){
  std::sort( barCol.begin(), barCol.end() );

  for(int ib=0; ib<barCol.size(); ib++){
    std::vector<CRDEcalEDM::CRDCaloBar> m_neighbors = getNeighbors( barCol[ib], barCol );
    if( m_neighbors.size()==0 && barCol[ib].getEnergy()>settings.Eth_localMax ) { localMaxCol.push_back( barCol[ib] ); continue; }

    bool isLocalMax=true;
    bool isIso;
    double Eneigh=0;
    if(barCol[ib].getEnergy()<settings.Eth_localMax) isLocalMax = false;
    for(int j=0;j<m_neighbors.size();j++){
      if(m_neighbors[j].getEnergy()>barCol[ib].getEnergy()) isLocalMax=false;
      Eneigh += m_neighbors[j].getEnergy();
    }
    isIso = (barCol[ib].getEnergy()/(barCol[ib].getEnergy()+Eneigh))>settings.Eth_MaxWithNeigh;
    if(!isIso) continue;

    if(isLocalMax) localMaxCol.push_back( barCol[ib] );

  }

  return StatusCode::SUCCESS;
}

#endif
