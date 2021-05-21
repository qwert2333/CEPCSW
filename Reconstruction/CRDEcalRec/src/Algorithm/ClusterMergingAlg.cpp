#ifndef _CLUSTERMERGING_ALG_C
#define _CLUSTERMERGING_ALG_C

#include "Algorithm/ClusterMergingAlg.h"

void ClusterMergingAlg::Settings::SetInitialValue(){
  Debug=0;
}

StatusCode ClusterMergingAlg::Initialize(){

  //Initialize settings
  return StatusCode::SUCCESS;
}

StatusCode ClusterMergingAlg::RunAlgorithm(ClusterMergingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){
  settings = m_settings;

  std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_goodClus = m_datacol.GoodClus3DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_badClus  = m_datacol.BadClus3DCol;
  if(settings.Debug>0) { m_datacol.PrintLayer(); m_datacol.PrintShower(); m_datacol.Print3DClus(); }

  for(int ib=0;ib<m_datacol.BlockVec.size(); ib++) m_datacol.BlockVec[ib].ClearCandidate();

  std::vector<int> m_DlayerVec = GetLayerNeedModification( m_goodClus, m_badClus );
  if(settings.Debug>0){ 
    std::cout<<"DEBUG: Layers need to modify: ";
    for(int il=0;il<m_DlayerVec.size();il++) std::cout<<m_DlayerVec[il]<<"   ";
    std::cout<<std::endl;
  }

  //Update CaloBlock with EMCandidate
  for(int i=0;i<m_DlayerVec.size();i++){
    int dlayer = m_DlayerVec[i];
    //std::cout<<"DEBUG: layer with ghost: "<<dlayer<<std::endl;
    std::vector<CRDEcalEDM::CRDShowerCandidate> m_expshvec; m_expshvec.clear();
    for(int icl=0; icl<m_goodClus.size(); icl++){
      CRDEcalEDM::CRDShowerCandidate m_expsh; m_expsh.Clear(); 
      m_expsh.Dlayer = dlayer; 
      m_expsh.ExpEshower = m_goodClus[icl].getExpEnergy(dlayer); 
      if(m_expsh.ExpEshower <0 ) continue; 
      m_expsh.ExpEseed = m_expsh.ExpEshower*0.8; 
      m_expsh.ExpDepth = dlayer; 
      m_expsh.ExpPos = m_goodClus[icl].getExpPos(dlayer); 
      m_expshvec.push_back(m_expsh);
    }

    if(settings.Debug>0){
      std::cout<<"DEBUG: #EMCandidate: "<<m_expshvec.size()<<std::endl;
      for(int iexp=0;iexp<m_expshvec.size(); iexp++) 
        printf("\t DEBUG: EMCandidate position (%.2f, %.2f, %.2f), Energy %.3f \n", 
               m_expshvec[iexp].ExpPos.x(), m_expshvec[iexp].ExpPos.y(),m_expshvec[iexp].ExpPos.z(), m_expshvec[iexp].ExpEshower );

    }

    for(int ib=0;ib<m_datacol.BlockVec.size(); ib++){
      if(m_datacol.BlockVec[ib].getDlayer()==dlayer){
        m_datacol.BlockVec[ib].setForced(true);
        m_datacol.BlockVec[ib].setCandidateCol(m_expshvec);
    }}
  }

  //m_datacol.Flag_Iter++;
  return StatusCode::SUCCESS;

}


std::vector<int> ClusterMergingAlg::GetGhostHitsLayer( std::vector<CRDEcalEDM::CRDCaloHit3DShower>& m_goodClus, std::vector<CRDEcalEDM::CRDCaloHit3DShower>& m_badClus ){

  std::vector<int> vec_dlayer; vec_dlayer.clear();

  std::vector<CRDEcalEDM::CRDCaloHit3DShower> ghostCandidates;
  for(int i=0;i<m_badClus.size();i++)  if(m_badClus[i].get2DShowers().size()==1) ghostCandidates.push_back(m_badClus[i]);

  //map<layer, 2DshowerCol>
  std::map<int, std::vector<CRDEcalEDM::CRDCaloHit2DShower> > map_goodshower; map_goodshower.clear();
  std::map<int, std::vector<CRDEcalEDM::CRDCaloHit2DShower> > map_ghostshower;  map_ghostshower.clear();

  for(int i=0;i<m_goodClus.size();i++){
    for(int j=0;j<m_goodClus[i].get2DShowers().size();j++){
      int dlayer = m_goodClus[i].get2DShowers()[j].getDlayer();
      map_goodshower[dlayer].push_back(m_goodClus[i].get2DShowers()[j]);
  }}
  for(int i=0;i<ghostCandidates.size();i++){
    int dlayer = ghostCandidates[i].get2DShowers()[0].getDlayer();
    map_ghostshower[dlayer].push_back(ghostCandidates[i].get2DShowers()[0]);
  }

  std::map<int, std::vector<CRDEcalEDM::CRDCaloHit2DShower> >::iterator iter = map_ghostshower.begin(); 

  for(iter; iter!=map_ghostshower.end(); iter++){
    int dlayer = iter->first;
    std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_showerinlayer = iter->second;
    if( (m_showerinlayer.size()+map_goodshower[dlayer].size()) <4) continue;
    vec_dlayer.push_back(dlayer);
  }

  return vec_dlayer;
}


std::vector<int> ClusterMergingAlg::GetLayerNeedModification( std::vector<CRDEcalEDM::CRDCaloHit3DShower>& m_goodClus, std::vector<CRDEcalEDM::CRDCaloHit3DShower>& m_badClus ){
  
  std::vector<int> vec_dlayer; vec_dlayer.clear();
  //Case0: no good cluster: Do nothing, skip this event. 
  if(m_goodClus.size()==0) return vec_dlayer;

  //Case1: only 1 good cluster, i.e. 1 photon: tune all fired layers with the good cluster. 
  else if(m_goodClus.size()==1){
    for(int ish=0; ish<m_goodClus[0].get2DShowers().size(); ish++) vec_dlayer.push_back( m_goodClus[0].get2DShowers()[ish].getDlayer() );
    for(int icl=0; icl<m_badClus.size(); icl++){
    for(int ish=0; ish<m_badClus[icl].get2DShowers().size(); ish++) vec_dlayer.push_back( m_badClus[icl].get2DShowers()[ish].getDlayer() );
    }
    std::sort(vec_dlayer.begin(), vec_dlayer.end());
    vec_dlayer.erase( unique(vec_dlayer.begin(), vec_dlayer.end()), vec_dlayer.end() );

    return vec_dlayer;
  }


  //Case2: multiple clusters: tune in all layers. 
  std::vector<int> layer_start; layer_start.clear(); 
  std::vector<int> layer_end; layer_end.clear(); 
  for(int icl=0; icl<m_goodClus.size(); icl++){ 
    layer_start.push_back( m_goodClus[icl].getBeginningDlayer() ); 
    layer_end.push_back( m_goodClus[icl].getEndDlayer() );
  }
  std::sort( layer_start.begin(), layer_start.end() );
  std::sort( layer_end.begin(), layer_end.end() );

  int start = layer_start[0];                
  int end   = layer_end[layer_end.size()-1]; 

  for(int i=start; i<end+1; i++) vec_dlayer.push_back(i);
  return vec_dlayer; 

}


StatusCode ClusterMergingAlg::ClearAlgorithm(){


  return StatusCode::SUCCESS;

}

#endif
