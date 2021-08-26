#ifndef _CANDIDATEMAKING_ALG_C
#define _CANDIDATEMAKING_ALG_C

#include "Algorithm/ShadowMakingAlg.h"
#include <set>

void ShadowMakingAlg::Settings::SetInitialValue(){
  Debug=0;
  UseTrk = false; 
  EndLayer = 15; 
}

StatusCode ShadowMakingAlg::Initialize(){

  //Initialize settings
  return StatusCode::SUCCESS;
}

StatusCode ShadowMakingAlg::RunAlgorithm(ShadowMakingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){
  settings = m_settings;

  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_goodClus = m_datacol.GoodClus3DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_badClus  = m_datacol.BadClus3DCol;
  if(settings.Debug>1) { m_datacol.PrintLayer(); m_datacol.PrintShower(); m_datacol.Print3DClus(); }

  for(int ib=0;ib<m_datacol.BlockVec.size(); ib++){ m_datacol.BlockVec[ib].ClearNeuCandidate(); m_datacol.BlockVec[ib].ClearTrkCandidate(); }

  std::vector<CRDEcalEDM::CRDCaloBlock> m_blockVec = GetBlocksNeedModification( m_goodClus, m_badClus ); //These blocks only have cellID info.

  if(settings.Debug>0){ 
    std::cout<<"DEBUG: Layers need to modify: ";
    for(int il=0;il<m_blockVec.size();il++) 
      printf("  (%d, %d, %d, %d) \n", m_blockVec[il].getModule(), m_blockVec[il].getStave(), m_blockVec[il].getDlayer(), m_blockVec[il].getPart() );
    std::cout<<std::endl;
  }


  //Update CaloBlock with Track
  if(m_datacol.TrackCol.size()!=0 && settings.UseTrk){
    for(int ib=0; ib<m_datacol.BlockVec.size(); ib++){
      CRDEcalEDM::CRDCaloBlock m_block = m_datacol.BlockVec[ib];
      if( m_block.getDlayer()>settings.EndLayer ) continue; 
      std::vector<CRDEcalEDM::CRDShadowCluster> m_exptrkvec; m_exptrkvec.clear();
      for(int it=0; it<m_block.getTrkCol().size(); it++){
        CRDEcalEDM::CRDShadowCluster m_trksh; m_trksh.Clear();
        m_trksh.Dlayer = m_block.getDlayer();
        m_trksh.ExpEshower = 0.01;
        m_trksh.ExpDepth = m_block.getDlayer();
        m_trksh.ExpPos = m_block.getTrkCol()[it].getProspectPos(m_block.getDlayer());
        m_trksh.Type = 1;
        m_exptrkvec.push_back(m_trksh);
      }
    m_datacol.BlockVec[ib].setTrkCandidateCol( m_exptrkvec );
    if(settings.Debug>0)cout<<"  Block #"<<ib<<": TrkCandidate number = "<<m_datacol.BlockVec[ib].getTrkCandidateCol().size()<<endl;
  }}



  //Update CaloBlock with Shower Candidate
  for(int ib=0; ib<m_datacol.BlockVec.size(); ib++){

    std::vector<CRDEcalEDM::CRDCaloBlock>::iterator iter = find(m_blockVec.begin(), m_blockVec.end(), m_datacol.BlockVec[ib]);
    if(iter == m_blockVec.end() ) continue; //This layer doesn't have candidate.

    std::vector<CRDEcalEDM::CRDShadowCluster> m_expshvec; m_expshvec.clear(); 
    for(int icl=0; icl<m_goodClus.size(); icl++){
      CRDEcalEDM::CRDShadowCluster m_neush; m_neush.Clear();
      m_neush.Dlayer = m_datacol.BlockVec[ib].getDlayer();
      m_neush.Type = 0; 
      m_neush.ExpEshower = m_goodClus[icl].getExpEnergy(m_neush.Dlayer);
      if(m_neush.ExpEshower <=0 ) continue;
      m_neush.ExpEseed = m_neush.ExpEshower*0.8;
      m_neush.ExpDepth = m_neush.Dlayer;
      m_neush.ExpPos = m_goodClus[icl].getExpPos(m_neush.Dlayer);

      //Overlap removal: remove the neutral candidate overlapped with track candidate.
      bool f_overlap = false; 
      CRDEcalEDM::CRDShadowCluster overlappedCandi; overlappedCandi.Clear(); 
      for(int icd=0; icd<m_datacol.BlockVec[ib].getTrkCandidateCol().size(); icd++){
        CRDEcalEDM::CRDShadowCluster m_trkCandi = m_datacol.BlockVec[ib].getTrkCandidateCol()[icd];
        if( (m_neush.ExpPos-m_trkCandi.ExpPos).Mag()<10 ){ f_overlap=true; overlappedCandi=m_trkCandi; break; }
      }
      if(f_overlap){ 
        if(settings.Debug>1) printf("One neutral candidate is overlapped with a trk candidate: NeuCandidate(%.2f, %.2f, %.2f), TrkCandidate(%.2f, %.2f, %.2f). \n", 
                                    m_neush.ExpPos.x(), m_neush.ExpPos.y(), m_neush.ExpPos.z(), 
                                    overlappedCandi.ExpPos.x(), overlappedCandi.ExpPos.y(), overlappedCandi.ExpPos.z() );
        continue; 
      }
      m_expshvec.push_back(m_neush);
    }

    m_datacol.BlockVec[ib].setNeuCandidateCol(m_expshvec);
    if(settings.Debug>0) cout<<"  Block #"<<ib<<": NeuCandidate number = "<<m_datacol.BlockVec[ib].getNeuCandidateCol().size()<<endl;
  }


  return StatusCode::SUCCESS;

}



std::vector<CRDEcalEDM::CRDCaloBlock> ShadowMakingAlg::GetBlocksNeedModification( std::vector<CRDEcalEDM::CRDCaloHit3DCluster>& m_goodClus, std::vector<CRDEcalEDM::CRDCaloHit3DCluster>& m_badClus ){
  
  std::vector<CRDEcalEDM::CRDCaloBlock> vec_blocks; vec_blocks.clear();
  //Case0: no good cluster: Do nothing, skip this event. 
  if(m_goodClus.size()==0) return vec_blocks;

  //All blocks in good cluster
  for(int icl=0; icl<m_goodClus.size(); icl++){
  for(int ish=0; ish<m_goodClus[icl].get2DShowers().size(); ish++){
    CRDEcalEDM::CRDCaloBlock m_block; m_block.Clear();
    m_block.setIDInfo( m_goodClus[icl].get2DShowers()[ish].getModule(), 
                       m_goodClus[icl].get2DShowers()[ish].getStave(), 
                       m_goodClus[icl].get2DShowers()[ish].getDlayer(), 
                       m_goodClus[icl].get2DShowers()[ish].getPart() );
    vec_blocks.push_back(m_block);
  }}
  //All blocks in bad cluster
  for(int icl=0; icl<m_badClus.size(); icl++){
  for(int ish=0; ish<m_badClus[icl].get2DShowers().size(); ish++){  
    CRDEcalEDM::CRDCaloBlock m_block; m_block.Clear();
    m_block.setIDInfo( m_badClus[icl].get2DShowers()[ish].getModule(),
                       m_badClus[icl].get2DShowers()[ish].getStave(),
                       m_badClus[icl].get2DShowers()[ish].getDlayer(),
                       m_badClus[icl].get2DShowers()[ish].getPart()    );
    vec_blocks.push_back(m_block);
  }}

  //Clear the duplicated blocks
  std::vector<CRDEcalEDM::CRDCaloBlock> vec_outblock; vec_outblock.clear();
  for(int ib=0; ib<vec_blocks.size(); ib++){
    int m_Dlayer = vec_blocks[ib].getDlayer(); 
    std::vector<CRDEcalEDM::CRDCaloBlock>::iterator iter = find(vec_outblock.begin(), vec_outblock.end(), vec_blocks[ib]);
    if( iter==vec_outblock.end() && m_Dlayer<=settings.EndLayer ) vec_outblock.push_back(vec_blocks[ib]);
  }

  return vec_outblock;


}


StatusCode ShadowMakingAlg::ClearAlgorithm(){


  return StatusCode::SUCCESS;

}

#endif
