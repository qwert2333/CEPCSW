#ifndef _CLUSTERMERGING_ALG_C
#define _CLUSTERMERGING_ALG_C

#include "Algorithm/ClusterMergingAlg.h"

void ClusterMergingAlg::Settings::SetInitialValue(){

}

StatusCode ClusterMergingAlg::Initialize(){

  //Initialize settings
/*  m_ESAlgSettings = new EnergySplittingAlg::Settings();
  m_ETAlgSettings = new EnergyTimeMatchingAlg::Settings();
  m_CCAlgSettings = new ConeClusteringAlg::Settings();

  //Initialize Algorthm
  m_energysplittingAlg = new EnergySplittingAlg();
  m_etmatchingAlg = new EnergyTimeMatchingAlg();
  m_coneclusterAlg = new ConeClusteringAlg();
*/

  return StatusCode::SUCCESS;
}

StatusCode ClusterMergingAlg::RunAlgorithm(ClusterMergingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){

//std::cout<<"DEBUG: Run ClusterMergingAlg"<<std::endl; 
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_goodClus = m_datacol.GoodClus3DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_badClus  = m_datacol.BadClus3DCol;

  //std::cout<<"#GoodClus and #BadClus: "<<m_goodClus.size()<<'\t'<<m_badClus.size()<<endl;

/*
//  std::vector<int> m_DlayerVec = GetGhostHitsLayer( m_goodClus, m_badClus );
  //std::cout<<"#GhostHitsLayer: "<<m_DlayerVec.size()<<std::endl;

  // No ghost hit layer, merge the bad cluster into closest good cluster
  if(m_DlayerVec.size()==0){ 
    m_datacol.ClearCluster(); 
    m_datacol.GoodClus3DCol = m_goodClus; 
    m_datacol.BadClus3DCol = m_badClus; 
    std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_mergedClus; m_mergedClus.clear(); 
    m_mergedClus = m_goodClus; 
    for(int icl=0;icl<m_badClus.size();icl++) if(!MergeToGoodCluster(m_mergedClus, m_badClus[icl], true))  cout<<"WARNING: Cluster merging fail!"<<endl;

    m_datacol.Clus3DCol = m_mergedClus;
    //std::cout<<"Cluster Merged! Go back to PandoraPlusPFAlg"<<std::endl;
    return StatusCode::SUCCESS;
  }
*/

  std::vector<int> m_DlayerVec = GetLayerNeedModification( m_goodClus, m_badClus );

//std::cout<<"Layers need modification: ";
//for(int i=0; i<m_DlayerVec.size(); i++) std::cout<<m_DlayerVec[i]<<'\t';
//std::cout<<std::endl;

  //Update CaloBlock with ExpShower
  for(int i=0;i<m_DlayerVec.size();i++){
    int dlayer = m_DlayerVec[i];
    //std::cout<<"DEBUG: layer with ghost: "<<dlayer<<std::endl;
    std::vector<CRDEcalEDM::CRDExpEMShower> m_expshvec; m_expshvec.clear();
    for(int icl=0; icl<m_goodClus.size(); icl++){
      CRDEcalEDM::CRDExpEMShower m_expsh; m_expsh.Clear(); 
      m_expsh.Dlayer = dlayer; 
      m_expsh.ExpEshower = m_goodClus[icl].getExpEnergy(dlayer); 
      if(m_expsh.ExpEshower <0 ) continue; 
      m_expsh.ExpEseed = m_expsh.ExpEshower*0.8; 
      m_expsh.ExpDepth = dlayer; 
      m_expsh.ExpPos = m_goodClus[icl].getExpPos(dlayer); 
      m_expshvec.push_back(m_expsh);
    }

//std::cout<<"DEBUG: #ExpShower: "<<m_expshvec.size()<<std::endl;
//for(int iexp=0;iexp<m_expshvec.size(); iexp++) printf("\t DEBUG: ExpShower position (%.2f, %.2f, %.2f), Energy %.3f \n", m_expshvec[iexp].ExpPos.x(), m_expshvec[iexp].ExpPos.y(),m_expshvec[iexp].ExpPos.z(), m_expshvec[iexp].ExpEshower );

    for(int ib=0;ib<m_datacol.blockVec.size(); ib++){
      if(m_datacol.blockVec[ib].getDlayer()==dlayer){
//std::cout<<"DEBUG: Found corresponding block!"<<std::endl;
        m_datacol.blockVec[ib].setForced(true);
        m_datacol.blockVec[ib].setMode(1);
        m_datacol.blockVec[ib].setExpShowerCol(m_expshvec);
//std::cout<<"DEBUG: Added ExpShowers into block!"<<std::endl;
    }}
  }


  //Do the iteration1
/*  m_datacol.ClearLayer(); 
  m_datacol.ClearShower(); 
  m_datacol.ClearCluster();
  m_datacol.ClearPFO(); 

  m_ESAlgSettings->SetInitialValue();
  m_ETAlgSettings->SetInitialValue();
  m_CCAlgSettings->SetInitialValue(); 
  m_ESAlgSettings->SetValues();
  m_CCAlgSettings->SetMergeBadClus(true);

  m_energysplittingAlg->RunAlgorithm( *m_ESAlgSettings, m_datacol );
  m_etmatchingAlg->     RunAlgorithm( *m_ETAlgSettings, m_datacol );
  m_coneclusterAlg->    RunAlgorithm( *m_CCAlgSettings, m_datacol); 
*/
  m_datacol.Flag_Iter = true;
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


  //Case2: multiple clusters: tune in overlapping layers. 
  std::vector<int> layer_start; layer_start.clear(); 
  std::vector<int> layer_end; layer_end.clear(); 
  for(int icl=0; icl<m_goodClus.size(); icl++){ 
    layer_start.push_back( m_goodClus[icl].getBeginningDlayer() ); 
    layer_end.push_back( m_goodClus[icl].getEndDlayer() );
  }
  std::sort( layer_start.begin(), layer_start.end() );
  std::sort( layer_end.begin(), layer_end.end() );

  int start = layer_start[1];                //second 
  int end   = layer_end[layer_end.size()-2]; //second last. Dlayer in [start, end] would have overlapping clusters.  

  for(int i=start-1; i<end+1; i++) vec_dlayer.push_back(i);
  return vec_dlayer; 

}

bool ClusterMergingAlg::MergeToGoodCluster( std::vector<CRDEcalEDM::CRDCaloHit3DShower>& goodClusCol,  CRDEcalEDM::CRDCaloHit3DShower& badClus, bool ForceMerging){

   CRDEcalEDM::CRDCaloHit3DShower goodClus = GetClosestGoodCluster(goodClusCol, badClus);

   //WIP: check the chi2 before/after megering.
   //double chi2_before;
   //double chi2_after;
   //if(!ForceMerging && chi2_after>chi2_before){
   // std::cout<<"WARNING: chi2 increases after merging. Skip this merging!"<<std<<endl;
   // return false;
   //}

   std::vector<CRDEcalEDM::CRDCaloHit3DShower>::iterator iter = find(goodClusCol.begin(), goodClusCol.end(), goodClus);
   if(iter==goodClusCol.end()) return false;
   else{
      iter->MergeCluster(badClus);
   }
   return true;
}


CRDEcalEDM::CRDCaloHit3DShower ClusterMergingAlg::GetClosestGoodCluster( std::vector<CRDEcalEDM::CRDCaloHit3DShower>& goodClusCol,  CRDEcalEDM::CRDCaloHit3DShower& badClus ){
   TVector3 m_clusCent = badClus.getShowerCenter();
   CRDEcalEDM::CRDCaloHit3DShower m_clusCandi;
   double minTheta=999;
   for(int i=0;i<goodClusCol.size();i++){
      TVector3 vec_goodClus = goodClusCol[i].getShowerCenter();
      TVector3 vec_goodAxis = goodClusCol[i].getAxis();
      double theta = vec_goodAxis.Cross(m_clusCent-vec_goodClus).Mag();
      if(theta<minTheta){
         minTheta=theta;
         m_clusCandi = goodClusCol[i];
      }
   }

   return m_clusCandi;
}


StatusCode ClusterMergingAlg::ClearAlgorithm(){

/*  delete m_energysplittingAlg;
  delete m_etmatchingAlg;
  delete m_coneclusterAlg;

  delete m_ESAlgSettings;
  delete m_ETAlgSettings;
  delete m_CCAlgSettings;
*/
  return StatusCode::SUCCESS;

}

#endif
