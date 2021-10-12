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

  std::vector<CRDEcalEDM::CRDShadowCluster> CRDCaloBlock::getAllShadowClusCol() const{
    std::vector<CRDEcalEDM::CRDShadowCluster> m_allcandi; m_allcandi.clear();
    m_allcandi = NeuShadowClusCol; 
    m_allcandi.insert(m_allcandi.end(), TrkShadowClusCol.begin(), TrkShadowClusCol.end());
    return m_allcandi; 
  }


  void CRDCaloBlock::PrintShadowClus() const{
    std::cout<<"    Track ShadowClus: # "<<TrkShadowClusCol.size()<<std::endl;
    for(int ic=0; ic<TrkShadowClusCol.size(); ic++)
    printf("\t  #%d: (%.3f, %.3f, %.3f), type = %d \n", 
           ic, TrkShadowClusCol[ic].ExpPos.x(), TrkShadowClusCol[ic].ExpPos.y(), TrkShadowClusCol[ic].ExpPos.z(), TrkShadowClusCol[ic].Type);
    std::cout<<"    Neutral ShadowClus: # "<<NeuShadowClusCol.size()<<std::endl;
    for(int ic=0; ic<NeuShadowClusCol.size(); ic++)
    printf("\t  #%d: (%.3f, %.3f, %.3f), type = %d \n",
           ic, NeuShadowClusCol[ic].ExpPos.x(), NeuShadowClusCol[ic].ExpPos.y(), NeuShadowClusCol[ic].ExpPos.z(), NeuShadowClusCol[ic].Type);
    
  }

};
#endif
