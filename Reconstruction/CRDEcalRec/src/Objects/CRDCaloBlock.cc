#ifndef _CRD_CALOBLOCK_C
#define _CRD_CALOBLOCK_C

#include "Objects/CRDCaloBlock.h"
#include "TMath.h"
namespace CRDEcalEDM{

  bool CRDCaloBlock::MatchTrk( CRDEcalEDM::Track _trk ) const {
    //For line: TODO: write this code
    //TVector3 vtx = _trk.getVertex(); 
    //TVector3 axis = _trk.getMomentum(); axis.SetMag(1.); 

    if(module==6 && stave==6 && part==2) return true;
    return false;
  }

  std::vector<CRDEcalEDM::CRDShadowCluster> CRDCaloBlock::getAllCandidateCol() const{
    std::vector<CRDEcalEDM::CRDShadowCluster> m_allcandi; m_allcandi.clear();
    m_allcandi = NeuCandidateCol; 
    m_allcandi.insert(m_allcandi.end(), TrkCandidateCol.begin(), TrkCandidateCol.end());
    return m_allcandi; 
  }


  void CRDCaloBlock::PrintCandidates() const{
    std::cout<<"    Track Candidate: # "<<TrkCandidateCol.size()<<std::endl;
    for(int ic=0; ic<TrkCandidateCol.size(); ic++)
    printf("\t  #%d: (%.3f, %.3f, %.3f), type = %d \n", 
           ic, TrkCandidateCol[ic].ExpPos.x(), TrkCandidateCol[ic].ExpPos.y(), TrkCandidateCol[ic].ExpPos.z(), TrkCandidateCol[ic].Type);
    std::cout<<"    Neutral Candidate: # "<<NeuCandidateCol.size()<<std::endl;
    for(int ic=0; ic<NeuCandidateCol.size(); ic++)
    printf("\t  #%d: (%.3f, %.3f, %.3f), type = %d \n",
           ic, NeuCandidateCol[ic].ExpPos.x(), NeuCandidateCol[ic].ExpPos.y(), NeuCandidateCol[ic].ExpPos.z(), NeuCandidateCol[ic].Type);
    
  }

};
#endif
