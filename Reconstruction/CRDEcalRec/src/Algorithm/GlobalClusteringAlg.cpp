#ifndef GLOBALCLUSTERING_ALG_C
#define GLOBALCLUSTERING_ALG_C

#include "Algorithm/GlobalClusteringAlg.h"
using namespace PandoraPlus;

StatusCode GlobalClusteringAlg::ReadSettings(PandoraPlus::Settings& m_settings){
  settings = m_settings;
  return StatusCode::SUCCESS;
};

StatusCode GlobalClusteringAlg::Initialize(){

  return StatusCode::SUCCESS;
};

StatusCode GlobalClusteringAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){
  //Readin: m_DataCol.BarCol
  //Output: m_Datacol.BlockCol

  std::vector<PandoraPlus::CaloBar*> m_bars = m_datacol.BarCol; 
  std::vector<PandoraPlus::CaloBlock*> m_blocks; m_blocks.clear();

  for(int ibar=0; ibar<m_bars.size(); ibar++){

    bool fl_foundbl = false;
    for(int ibl=0; ibl<m_blocks.size(); ibl++)
      if( m_bars[ibar]->getModule() == m_blocks[ibl]->getModule() &&
          m_bars[ibar]->getStave()  == m_blocks[ibl]->getStave() &&
          m_bars[ibar]->getPart()   == m_blocks[ibl]->getPart() &&
          m_bars[ibar]->getDlayer() == m_blocks[ibl]->getDlayer() ){
        m_blocks[ibl]->addBar( m_bars[ibar] );
        fl_foundbl=true;
      }

    if(!fl_foundbl){
      PandoraPlus::CaloBlock* m_block = new PandoraPlus::CaloBlock();
      m_datacol.bk_BlockCol.push_back(m_block);

      m_block->addBar( m_bars[ibar] );
      m_block->setIDInfo(m_bars[ibar]->getModule(),  m_bars[ibar]->getStave(), m_bars[ibar]->getDlayer(), m_bars[ibar]->getPart());
      m_blocks.push_back(m_block);
    }
  }
  m_datacol.BlockCol = m_blocks;

  return StatusCode::SUCCESS;
};

StatusCode GlobalClusteringAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
};


#endif

