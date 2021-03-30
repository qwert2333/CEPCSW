#ifndef _CONECLUSTERING_ALG_C
#define _CONECLUSTERING_ALG_C

#include "Algorithm/ConeClusteringAlg.h"

ConeClusteringAlg::ConeClusteringAlg(){

}


void ConeClusteringAlg::Settings::SetInitialValue(){
  th_ConeTheta = PI/6.; 
  th_ConeR = 50.; //50mm
  th_ClusChi2 = 10e17;
  fl_MergeBadClus=false; 
}


StatusCode ConeClusteringAlg::Initialize(){
  return StatusCode::SUCCESS;
}


StatusCode ConeClusteringAlg::RunAlgorithm( ConeClusteringAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){
  settings = m_settings;


  std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_ClusterCol;  m_ClusterCol.clear();
  std::map<int, std::vector<CRDEcalEDM::CRDCaloHit2DShower> > m_orderedShower;  //map<layer, showers>
  m_orderedShower.clear();

  std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_2DshowerCol = m_datacol.shower2DCol;
  for(int is=0;is<m_2DshowerCol.size();is++){
     m_orderedShower[m_2DshowerCol[is].getDlayer()].push_back(m_2DshowerCol[is]);
  }
 

 std::map<int, std::vector<CRDEcalEDM::CRDCaloHit2DShower>>::iterator iter = m_orderedShower.begin();

  //In first layer: initial clusters. All showers in the first layer are regarded as cluster seed.
  //cluster initial direction = R.
  std::vector<CRDEcalEDM::CRDCaloHit2DShower> ShowersinFirstLayer;  ShowersinFirstLayer.clear();
  ShowersinFirstLayer = iter->second;
  for(int i=0;i<ShowersinFirstLayer.size(); i++){
    CRDEcalEDM::CRDCaloHit3DShower m_clus;
    m_clus.AddShower(ShowersinFirstLayer[i]);
    m_clus.FitAxis();
    m_ClusterCol.push_back(m_clus);
  }
  iter++;

  //Loop later layers
  for(iter;iter!=m_orderedShower.end();iter++){
    std::vector<CRDEcalEDM::CRDCaloHit2DShower> ShowersinLayer = iter->second;
    for(int is=0; is<ShowersinLayer.size(); is++){
      CRDEcalEDM::CRDCaloHit2DShower m_shower = ShowersinLayer[is];

      for(int ic=0; ic<m_ClusterCol.size(); ic++ ){
        CRDEcalEDM::CRDCaloHit2DShower shower_in_clus = m_ClusterCol[ic].get2DShowers().back();
        dd4hep::Position relR = m_shower.getPos()-shower_in_clus.getPos();
        TVector3 relR_vec(relR.x(), relR.y(), relR.z());
        if(relR_vec.Angle(m_ClusterCol[ic].getAxis())< settings.th_ConeTheta && relR_vec.Mag()< settings.th_ConeR ){
          m_ClusterCol[ic].AddShower(m_shower);
          ShowersinLayer.erase(ShowersinLayer.begin()+is);
          is--;
          break;
        }
      }

    }//end loop showers in layer.
      if(ShowersinLayer.size()>0){
      for(int i=0;i<ShowersinLayer.size(); i++){
        CRDEcalEDM::CRDCaloHit3DShower m_clus;
        m_clus.AddShower(ShowersinLayer[i]);
        m_clus.FitAxis();
        m_ClusterCol.push_back(m_clus);
      }}//end new cluster
  }//end loop layers.


  //Check cluster quality, merge bad clusters to good clusters.
  std::vector< CRDEcalEDM::CRDCaloHit3DShower >  goodClus;
  std::vector< CRDEcalEDM::CRDCaloHit3DShower >  badClus;
  for(int icl=0; icl<m_ClusterCol.size(); icl++){
    if(CheckClusterQuality(m_ClusterCol[icl])) goodClus.push_back(m_ClusterCol[icl]);
    else badClus.push_back(m_ClusterCol[icl]);
  }
//std::cout<<"Number of good cluster: "<<goodClus.size()<<std::endl;
//std::cout<<"Number of bad cluster:  "<<badClus.size()<<std::endl;

  m_datacol.GoodClus3DCol = goodClus;
  m_datacol.BadClus3DCol = badClus;

  //Merge the bad clusters into closest good cluster
  if(settings.fl_MergeBadClus) 
    for(int icl=0;icl<badClus.size();icl++) if(!MergeToGoodCluster(goodClus, badClus[icl], true))  cout<<"WARNING: Cluster merging fail!"<<endl;
  
  m_datacol.Clus3DCol = goodClus;

  return StatusCode::SUCCESS;
}


bool ConeClusteringAlg::CheckClusterQuality(CRDEcalEDM::CRDCaloHit3DShower& clus){

  //Criterion1: at least 3 layers. 
  if(clus.get2DShowers().size()<=2 ) { return false; }

  //Criterion2: first 3 layser should satisfy: (E2/E1>3 && E3/E2>1.2) || (E3/E2>2). 
  bool isSeed;
  std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_2Dshowers = clus.get2DShowers();
  std::map<int, CRDEcalEDM::CRDCaloHit2DShower > m_orderedShower;  //map<layer, showers>
  for(int is=0;is<m_2Dshowers.size();is++){
     m_orderedShower[m_2Dshowers[is].getDlayer()] = m_2Dshowers[is];
  }
  std::map<int, CRDEcalEDM::CRDCaloHit2DShower>::iterator iter = m_orderedShower.begin();
  CRDEcalEDM::CRDCaloHit2DShower showerbegin[3];
  showerbegin[0] = iter->second;  iter++; 
  showerbegin[1] = iter->second;  iter++; 
  showerbegin[2] = iter->second;
  if( (showerbegin[1].getShowerE()/showerbegin[0].getShowerE()>3 && showerbegin[2].getShowerE()/showerbegin[1].getShowerE()>1.2) ||
       showerbegin[2].getShowerE()/showerbegin[1].getShowerE()>2 )
  isSeed = true; 
  if(!isSeed) return false;

  //Criterion3: profile fit chi2 < th_chi2
  if( clus.getChi2() > settings.th_ClusChi2 ) return false;

  return true;
}


bool ConeClusteringAlg::MergeToGoodCluster( std::vector<CRDEcalEDM::CRDCaloHit3DShower>& goodClusCol,  CRDEcalEDM::CRDCaloHit3DShower& badClus, bool ForceMerging){

   CRDEcalEDM::CRDCaloHit3DShower m_goodClus = GetClosestGoodCluster(goodClusCol, badClus);

   //WIP: check the chi2 before/after megering.
   //double chi2_before;
   //double chi2_after;
   //if(!ForceMerging && chi2_after>chi2_before){
   // std::cout<<"WARNING: chi2 increases after merging. Skip this merging!"<<std<<endl;
   // return false;
   //}

   std::vector<CRDEcalEDM::CRDCaloHit3DShower>::iterator iter = find(goodClusCol.begin(), goodClusCol.end(), m_goodClus);

   if(iter==goodClusCol.end()) return false;
   else{
      iter->MergeCluster(badClus);
   }
   return true;
}


CRDEcalEDM::CRDCaloHit3DShower ConeClusteringAlg::GetClosestGoodCluster( std::vector<CRDEcalEDM::CRDCaloHit3DShower>& goodClusCol,  CRDEcalEDM::CRDCaloHit3DShower& badClus ){

   TVector3 m_clusCent = badClus.getShowerCenter();
   CRDEcalEDM::CRDCaloHit3DShower m_clusCandi; m_clusCandi.Clear();
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

#endif
