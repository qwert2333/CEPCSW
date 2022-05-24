#ifndef ECALHIT_CREATOR_H
#define ECALHIT_CREATOR_H

#include "k4FWCore/DataHandle.h"
#include "PandoraPlusDataCol.h"
#include "CRDEcalSvc/ICRDEcalSvc.h"

class EcalHitsCreator{

public: 

  class Settings{
  public: 
    Settings(){};

  };
  
  //initialize a CaloHitCreator
  EcalHitsCreator( Settings& settings ){};
  ~EcalHitsCreator() {};

  StatusCode ReadSettings( Settings& settings ){ return StatusCode::SUCCESS; };

  StatusCode GetEcalBars(PandoraPlusDataCol& dataCol , ICRDEcalSvc& m_svc ){ 
    std::vector<CRDEcalEDM::CRDCaloBar> m_bars; m_bars.clear();
    m_svc.getDigiHits( m_bars );  

    std::vector<CRDEcalEDM::CRDCaloBlock> m_blocks; m_blocks.clear();
    
    for(int ibar=0; ibar<m_bars.size(); ibar++){

      bool fl_foundbl = false;
      for(int ibl=0; ibl<m_blocks.size(); ibl++)
        if( m_bars[ibar].getModule() == m_blocks[ibl].getModule() &&  
            m_bars[ibar].getStave() == m_blocks[ibl].getStave() &&
            m_bars[ibar].getPart() == m_blocks[ibl].getPart() &&
            m_bars[ibar].getDlayer() == m_blocks[ibl].getDlayer() ){
          m_blocks[ibl].addBar( m_bars[ibar] );
          fl_foundbl=true; 
        }

      if(!fl_foundbl){
        CRDEcalEDM::CRDCaloBlock m_block; m_block.Clear();
        m_block.addBar( m_bars[ibar]);
        m_block.setIDInfo(m_bars[ibar].getModule(),  m_bars[ibar].getStave(), m_bars[ibar].getDlayer(), m_bars[ibar].getPart());
        m_blocks.push_back(m_block);
      }
    }

    dataCol.BlockVec = m_blocks; 
    return StatusCode::SUCCESS;
  };

  void Reset(){};

  Settings  settings; 

private: 


};
#endif
