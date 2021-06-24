#ifndef _CANDIDATEMAKING_ALG_C
#define _CANDIDATEMAKING_ALG_C

#include "Algorithm/CandidateMakingAlg.h"
#include <set>

void CandidateMakingAlg::Settings::SetInitialValue(){
  Debug=0;
  UseTrk = false; 
}

StatusCode CandidateMakingAlg::Initialize(){

  //Initialize settings
  return StatusCode::SUCCESS;
}

StatusCode CandidateMakingAlg::RunAlgorithm(CandidateMakingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){
  settings = m_settings;

  std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_goodClus = m_datacol.GoodClus3DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_badClus  = m_datacol.BadClus3DCol;
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
      std::vector<CRDEcalEDM::CRDShowerCandidate> m_exptrkvec; m_exptrkvec.clear();
      for(int it=0; it<m_block.getTrkCol().size(); it++){
        CRDEcalEDM::CRDShowerCandidate m_trksh; m_trksh.Clear();
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

    bool f_match = false; 
    for(int jb=0; jb<m_blockVec.size(); jb++) if(m_datacol.BlockVec[ib]==m_blockVec[jb]){ f_match=true; break; } 
    if(!f_match) continue;  //This layer doesn't have candidate. 

    std::vector<CRDEcalEDM::CRDShowerCandidate> m_expshvec; m_expshvec.clear(); 
    for(int icl=0; icl<m_goodClus.size(); icl++){
      CRDEcalEDM::CRDShowerCandidate m_neush; m_neush.Clear();
      m_neush.Dlayer = m_datacol.BlockVec[ib].getDlayer();
      m_neush.Type = 0; 
      m_neush.ExpEshower = m_goodClus[icl].getExpEnergy(m_neush.Dlayer);
      if(m_neush.ExpEshower <=0 ) continue;
      m_neush.ExpEseed = m_neush.ExpEshower*0.8;
      m_neush.ExpDepth = m_neush.Dlayer;
      m_neush.ExpPos = m_goodClus[icl].getExpPos(m_neush.Dlayer);

      //Overlap removal: remove the neutral candidate overlapped with track candidate.
      bool f_overlap = false; 
      for(int icd=0; icd<m_datacol.BlockVec[ib].getTrkCandidateCol().size(); icd++){
        CRDEcalEDM::CRDShowerCandidate m_trkCandi = m_datacol.BlockVec[ib].getTrkCandidateCol()[icd];
        if( (m_neush.ExpPos-m_trkCandi.ExpPos).Mag()<10 ){ f_overlap=true; break; }
      }
      if(f_overlap) continue; 

      m_expshvec.push_back(m_neush);
    }

    m_datacol.BlockVec[ib].setNeuCandidateCol(m_expshvec);
    if(settings.Debug>0)cout<<"  Block #"<<ib<<": NeuCandidate number = "<<m_datacol.BlockVec[ib].getNeuCandidateCol().size()<<endl;
  }


  return StatusCode::SUCCESS;

}



std::vector<CRDEcalEDM::CRDCaloBlock> CandidateMakingAlg::GetBlocksNeedModification( std::vector<CRDEcalEDM::CRDCaloHit3DShower>& m_goodClus, std::vector<CRDEcalEDM::CRDCaloHit3DShower>& m_badClus ){
  
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
    bool f_exist = false;
    for(int jb=0; jb<vec_outblock.size(); jb++)
      if(vec_blocks[ib]==vec_outblock[jb]) { f_exist=true; break; }

    if(!f_exist) vec_outblock.push_back(vec_blocks[ib]);
  }

  return vec_outblock;


}


StatusCode CandidateMakingAlg::ClearAlgorithm(){


  return StatusCode::SUCCESS;

}

#endif
