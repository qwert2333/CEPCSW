#ifndef _SHADOWMAKING_ALG_C
#define _SHADOWMAKING_ALG_C

#include "Algorithm/ShadowMakingAlg.h"

void ShadowMakingAlg::Settings::SetInitialValue(){
  fl_UseTrack = true; 
  EndLayer = 15;
  th_GoodLayer = 4; 
}

StatusCode ShadowMakingAlg::Initialize(){

  return StatusCode::SUCCESS;
}

StatusCode ShadowMakingAlg::RunAlgorithm(ShadowMakingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){
  settings = m_settings;

  std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_clusXCol = m_datacol.LongiClusXCol; 
  std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_clusYCol = m_datacol.LongiClusYCol; 

  for(int it=0; it<m_datacol.TowerCol.size(); it++){
    std::vector<CRDEcalEDM::CRDCaloBlock> m_blocksInTower = m_datacol.TowerCol[it].getBlocks();

    for(int ib=0;ib<m_blocksInTower.size(); ib++){ m_blocksInTower[ib].ClearNeuShadowClus(); m_blocksInTower[ib].ClearTrkShadowClus(); }
   
    //Create shadow cluster from track: 
    if(m_datacol.TrackCol.size()!=0 && settings.fl_UseTrack){
      for(int ib=0; ib<m_blocksInTower.size(); ib++){
        CRDEcalEDM::CRDCaloBlock m_block = m_blocksInTower[ib];
        if( m_block.getDlayer()>settings.EndLayer ) continue;
        std::vector<CRDEcalEDM::CRDShadowCluster> m_exptrkvec; m_exptrkvec.clear();
        for(int it=0; it<m_block.getTrkCol().size(); it++){
          CRDEcalEDM::CRDShadowCluster m_trksh; m_trksh.Clear();
          m_trksh.Dlayer = m_block.getDlayer();
          m_trksh.ExpEshower = 0.01;
          m_trksh.ExpDepth = m_block.getDlayer();
          m_trksh.ExpPos = m_block.getTrkCol()[it].getProspectPos(m_block.getDlayer());  //FIXME: only work for (module, stave, part)=(6,6,2). Need full detector extrapolation. 
          m_trksh.Type = 1;
          m_exptrkvec.push_back(m_trksh);
        }
      m_blocksInTower[ib].setTrkShadowClusCol( m_exptrkvec );
    }}  
   
   
    //Identify the cluster quality: 
    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_goodClusX; m_goodClusX.clear(); 
    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_badClusX; m_badClusX.clear(); 
    for(int ic=0; ic<m_clusXCol.size(); ic++){
      if( m_clusXCol[ic].getEndDlayer()-m_clusXCol[ic].getBeginningDlayer()>=settings.th_GoodLayer ) m_goodClusX.push_back(m_clusXCol[ic]);
      else m_badClusX.push_back(m_clusXCol[ic]);
    }
   
    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_goodClusY; m_goodClusY.clear(); 
    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_badClusY; m_badClusY.clear(); 
    for(int ic=0; ic<m_clusYCol.size(); ic++){
      if( m_clusYCol[ic].getEndDlayer()-m_clusYCol[ic].getBeginningDlayer()>=settings.th_GoodLayer ) m_goodClusY.push_back(m_clusYCol[ic]);
      else m_badClusY.push_back(m_clusYCol[ic]);
    }
   
    //Create shadow cluster from good cluster: 
    for(int ic=0; ic<m_goodClusX.size(); ic++){ 
      std::vector<CRDEcalEDM::CRDCaloBlock> m_blocks = GetBlocksNeedModification(m_goodClusX[ic]);  //Only have block ID
   
      for(int ib=0; ib<m_blocksInTower.size(); ib++){
        if(find( m_blocks.begin(), m_blocks.end(), m_blocksInTower[ib] )==m_blocks.end()) continue;
   
        CRDEcalEDM::CRDShadowCluster m_neush; m_neush.Clear();
        m_neush.Dlayer = m_blocksInTower[ib].getDlayer();
        m_neush.Type = 0;
        m_neush.slayer = 0; 
        m_neush.ExpPos = m_goodClusX[ic].getExpPos(m_neush.Dlayer);
   
        //Check if this neutral shadow cluster is overlap with a track shadow cluster
        bool f_overlap = false; 
        CRDEcalEDM::CRDShadowCluster overlappedSdClus; overlappedSdClus.Clear();
        for(int icd=0; icd<m_blocksInTower[ib].getTrkShadowClusCol().size(); icd++){
          CRDEcalEDM::CRDShadowCluster m_trkSdClus = m_blocksInTower[ib].getTrkShadowClusCol()[icd];
          if( (m_neush.ExpPos-m_trkSdClus.ExpPos).Mag()<10 ){ f_overlap=true; overlappedSdClus=m_trkSdClus; break; }
        }
        if(!f_overlap) m_blocksInTower[ib].addNeuShadowClus( m_neush );
      }
    }
   
   
    for(int ic=0; ic<m_goodClusY.size(); ic++){
      std::vector<CRDEcalEDM::CRDCaloBlock> m_blocks = GetBlocksNeedModification(m_goodClusY[ic]);
   
      for(int ib=0; ib<m_blocksInTower.size(); ib++){
        if(find( m_blocks.begin(), m_blocks.end(), m_blocksInTower[ib] )==m_blocks.end()) continue;
   
        CRDEcalEDM::CRDShadowCluster m_neush; m_neush.Clear();
        m_neush.Dlayer = m_blocksInTower[ib].getDlayer();
        m_neush.Type = 0;
        m_neush.ExpPos = m_goodClusY[ic].getExpPos(m_neush.Dlayer);
   
        bool f_overlap = false;
        CRDEcalEDM::CRDShadowCluster overlappedSdClus; overlappedSdClus.Clear();
        for(int icd=0; icd<m_blocksInTower[ib].getTrkShadowClusCol().size(); icd++){
          CRDEcalEDM::CRDShadowCluster m_trkSdClus = m_blocksInTower[ib].getTrkShadowClusCol()[icd];
          if( (m_neush.ExpPos-m_trkSdClus.ExpPos).Mag()<10 ){ f_overlap=true; overlappedSdClus=m_trkSdClus; break; }
        }
        if(!f_overlap) m_blocksInTower[ib].addNeuShadowClus( m_neush );
      }
    }

    dataCol.TowerCol[it].SetBlocks( m_blocksInTower );
  }
  return StatusCode::SUCCESS;
}

StatusCode ShadowMakingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}


std::vector<CRDEcalEDM::CRDCaloBlock> ShadowMakingAlg::GetBlocksNeedModification( CRDEcalEDM::CRDCaloHitLongiCluster& m_clus ){

  std::vector<CRDEcalEDM::CRDCaloBlock> vec_blocks; vec_blocks.clear();

  for(int is=0; is<m_clus.getBarShowers().size(); is++){
    CRDEcalEDM::CRDCaloBlock m_block; m_block.Clear();
    m_block.setIDInfo( m_clus.getBarShowers()[is].getModule(),  
                       m_clus.getBarShowers()[is].getStave(),
                       m_clus.getBarShowers()[is].getDlayer(), 
                       m_clus.getBarShowers()[is].getPart()  );
    vec_blocks.push_back(m_block);
  }

  return vec_blocks; 
}

#endif
