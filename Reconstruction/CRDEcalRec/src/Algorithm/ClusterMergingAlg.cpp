#ifndef _CLUSTERMERGING_ALG_C
#define _CLUSTERMERGING_ALG_C

#include "Algorithm/ClusterMergingAlg.h"

void ClusterMergingAlg::Settings::SetInitialValue(){
  axis_Angle = PI/6.; 
  relP_Angle = PI/5.; 
  skipLayer = 3;
  fl_MergeGoodClus = true;
  fl_MergeBadClus = true;
  fl_MergeEMTail = true; 
  Debug = 0;
};

StatusCode ClusterMergingAlg::Initialize(){

  return StatusCode::SUCCESS;
};

StatusCode ClusterMergingAlg::RunAlgorithm( ClusterMergingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){
  settings = m_settings; 

  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_goodClusCol = m_datacol.GoodClus3DCol; 
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_badClusCol  = m_datacol.BadClus3DCol; 
  if(settings.Debug>0){ 
    int Ntots = 0; 
    for(int is=0; is<m_goodClusCol.size(); is++) Ntots += m_goodClusCol[is].get2DShowers().size(); 
    std::cout<<"GoodCluster number: "<<m_goodClusCol.size()<<"  Total showers in good cluster: "<<Ntots<<std::endl;
    Ntots = 0; 
    for(int is=0; is<m_badClusCol.size(); is++) Ntots += m_badClusCol[is].get2DShowers().size(); 
    std::cout<<"BadCluster number: "<<m_badClusCol.size()<<"  Total showers in bad cluster: "<<Ntots<<std::endl;
  }

  //std::vector<CRDEcalEDM::Track>              m_trkCol  = dataCol.TrackCol;

  //Merge 2 good clusters with several criteria
  if(settings.fl_MergeGoodClus){
  std::sort(m_goodClusCol.begin(), m_goodClusCol.end(), compBegin); 
  for(int ic=0; ic<m_goodClusCol.size() && m_goodClusCol.size()>1; ic++){
    CRDEcalEDM::CRDCaloHit3DCluster m_clus = m_goodClusCol[ic]; 
    for(int jc=ic+1; jc<m_goodClusCol.size(); jc++){
      CRDEcalEDM::CRDCaloHit3DCluster p_clus = m_goodClusCol[jc];
      double minRelAngle = min( sin((m_clus.getShowerCenter()-p_clus.getShowerCenter()).Angle(m_clus.getAxis())), 
                                sin((m_clus.getShowerCenter()-p_clus.getShowerCenter()).Angle(p_clus.getAxis())) );
      if( !(p_clus.getBeginningDlayer()>m_clus.getEndDlayer() || m_clus.getBeginningDlayer()>p_clus.getEndDlayer()) ) continue;

      if( sin( m_clus.getAxis().Angle(p_clus.getAxis()) ) < sin(settings.axis_Angle) && 
          minRelAngle < sin(settings.relP_Angle) &&
          p_clus.getBeginningDlayer()-m_clus.getEndDlayer() <= settings.skipLayer 
        ) {
        m_goodClusCol[ic].MergeCluster( m_goodClusCol[jc] ); 
        m_goodClusCol.erase(m_goodClusCol.begin()+jc);
        ic--; jc--;
        break; 
      }
    }
  }}
  

  //TODO: consider different cluster type cases.
  //Divide good clusters with cluster type
  //std::map< int, std::vector<CRDEcalEDM::CRDCaloHit3DCluster>> map_goodClusType; map_goodClusType.clear();
  //for(int ic=0; ic<m_goodClusCol.size(); ic++)
  //  map_goodClusType[ m_goodClusCol[ic].getType() ].push_back( m_goodClusCol[ic] );

  //Case1: MIP tracks
  //Case2: EM cluster
  //Case3: others

  //Merge bad clusters into closest good cluster
  if(settings.fl_MergeBadClus){
    for(int icl=0;icl<m_badClusCol.size();icl++) if(!MergeToGoodCluster(m_goodClusCol, m_badClusCol[icl], true))  cout<<"WARNING: Cluster merging fail!"<<endl;
  }

  //Merge 2 good clusters with several criteria
  if(settings.fl_MergeGoodClus){
  std::sort(m_goodClusCol.begin(), m_goodClusCol.end(), compBegin);
  for(int ic=0; ic<m_goodClusCol.size() && m_goodClusCol.size()>1; ic++){
    CRDEcalEDM::CRDCaloHit3DCluster m_clus = m_goodClusCol[ic];
    for(int jc=ic+1; jc<m_goodClusCol.size(); jc++){
      CRDEcalEDM::CRDCaloHit3DCluster p_clus = m_goodClusCol[jc];
      if( !(p_clus.getBeginningDlayer()>m_clus.getEndDlayer() || m_clus.getBeginningDlayer()>p_clus.getEndDlayer()) ) continue;

      double minRelAngle = min( sin((m_clus.getShowerCenter()-p_clus.getShowerCenter()).Angle(m_clus.getAxis())), 
                                sin((m_clus.getShowerCenter()-p_clus.getShowerCenter()).Angle(p_clus.getAxis())) );
      if( sin( m_clus.getAxis().Angle(p_clus.getAxis()) ) < sin(settings.axis_Angle) &&
          minRelAngle < sin(settings.relP_Angle) &&
          p_clus.getBeginningDlayer()-m_clus.getEndDlayer() <= settings.skipLayer
        ) {
        m_goodClusCol[ic].MergeCluster( m_goodClusCol[jc] );
        m_goodClusCol.erase(m_goodClusCol.begin()+jc);
        ic--; jc--;
        break;
      }
    }
  }}


  //Merge EM cluster tail into EM core by profile
  if(settings.fl_MergeEMTail){
  for(int ic=0; ic<m_goodClusCol.size() && m_goodClusCol.size()>1; ic++){
    m_goodClusCol[ic].IdentifyCluster();
    if(m_goodClusCol[ic].getType()!=1) continue; 
    //m_goodClusCol[ic].FitProfile(); 
    for(int jc=ic+1; jc<m_goodClusCol.size(); jc++){
      if( fabs(m_goodClusCol[ic].getFitExpEn()-m_goodClusCol[ic].getShowerE()-m_goodClusCol[jc].getShowerE()) < fabs(m_goodClusCol[ic].getFitExpEn()-m_goodClusCol[ic].getShowerE()) ){
        m_goodClusCol[ic].MergeCluster( m_goodClusCol[jc] );
        m_goodClusCol.erase(m_goodClusCol.begin()+jc);
        ic--; jc--;
        break;
      }
    }  
  }}

  m_datacol.GoodClus3DCol = m_goodClusCol;
  if(settings.fl_MergeBadClus || settings.fl_MergeGoodClus) m_datacol.Clus3DCol = m_goodClusCol;


  return StatusCode::SUCCESS;
};

bool ClusterMergingAlg::MergeToGoodCluster( std::vector<CRDEcalEDM::CRDCaloHit3DCluster>& goodClusCol,  CRDEcalEDM::CRDCaloHit3DCluster& badClus, bool ForceMerging){

   CRDEcalEDM::CRDCaloHit3DCluster m_goodClus = GetClosestGoodCluster(goodClusCol, badClus);

   //WIP: check the chi2 before/after megering.
   //double chi2_before;
   //double chi2_after;
   //if(!ForceMerging && chi2_after>chi2_before){
   // std::cout<<"WARNING: chi2 increases after merging. Skip this merging!"<<std<<endl;
   // return false;
   //}

   std::vector<CRDEcalEDM::CRDCaloHit3DCluster>::iterator iter = find(goodClusCol.begin(), goodClusCol.end(), m_goodClus);

   if(iter==goodClusCol.end()) return false;
   else{
      iter->MergeCluster(badClus);
   }
   return true;
};


CRDEcalEDM::CRDCaloHit3DCluster ClusterMergingAlg::GetClosestGoodCluster( std::vector<CRDEcalEDM::CRDCaloHit3DCluster>& goodClusCol,  CRDEcalEDM::CRDCaloHit3DCluster& badClus ){

   TVector3 m_clusCent = badClus.getShowerCenter();
   CRDEcalEDM::CRDCaloHit3DCluster m_clusCandi; m_clusCandi.Clear();
   double minTheta=999;
   for(int i=0;i<goodClusCol.size();i++){
      TVector3 vec_goodClus = goodClusCol[i].getShowerCenter();
      TVector3 vec_goodAxis = goodClusCol[i].getAxis();
      double theta = vec_goodAxis.Cross(m_clusCent-vec_goodClus).Mag();

      if(theta<minTheta){
         minTheta=theta;
         m_clusCandi.Clear();
         m_clusCandi = goodClusCol[i];
      }
   }

   return m_clusCandi;
}

StatusCode ClusterMergingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}

#endif
