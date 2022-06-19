#ifndef _CONECLUSTERING_ALG_C
#define _CONECLUSTERING_ALG_C

#include "Algorithm/ConeClusteringAlg.h"

ConeClusteringAlg::ConeClusteringAlg(){

}


void ConeClusteringAlg::Settings::SetInitialValue(){
  th_ConeTheta_l1 = PI/6.; 
  th_ConeR_l1 = 30.; //30mm
  th_ConeTheta_l2 = PI/10.;  
  th_ConeR_l2 = 30.; //30mm
  th_ClusChi2 = 10e17;
  fl_GoodClusLevel = 2;
  fl_UseCandidate=0; 
  clusType = ""; 
}

void ConeClusteringAlg::Settings::SetConeValue( double _coneTheta_l1, double _coneR_l1, double _coneTheta_l2, double _coneR_l2){
  th_ConeTheta_l1 = _coneTheta_l1; 
  th_ConeR_l1     = _coneR_l1;
  th_ConeTheta_l2 = _coneTheta_l2;
  th_ConeR_l2     = _coneR_l2;
}

StatusCode ConeClusteringAlg::Initialize(){
  return StatusCode::SUCCESS;
}


StatusCode ConeClusteringAlg::RunAlgorithm( ConeClusteringAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){
  settings = m_settings;

  std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_2DshowerCol; m_2DshowerCol.clear();
  if( settings.clusType=="MIP" )     m_2DshowerCol = m_datacol.MIPShower2DCol;
  else if( settings.clusType=="EM" ) m_2DshowerCol = m_datacol.EMShower2DCol;
  else m_2DshowerCol = m_datacol.Shower2DCol;

  if(m_2DshowerCol.size()==0){ 
    std::cout<<"Warning: Empty input in ConeClusteringAlg. Please check previous algorithm!"<<endl;  
    if(settings.fl_overwrite) m_datacol.ClearCluster();
    return StatusCode::SUCCESS; 
  }

//for(int is=0; is<m_2DshowerCol.size(); is++)
//  printf("      Shower %d: (%.1f, %.1f, %.1f, %.3f) \n", is,
//                                                         m_2DshowerCol[is].getPos().x(),
//                                                         m_2DshowerCol[is].getPos().y(),
//                                                         m_2DshowerCol[is].getPos().z(),
//                                                         m_2DshowerCol[is].getShowerE()  );


  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_trkClusterCol;  m_trkClusterCol.clear();
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_neuClusterCol;  m_neuClusterCol.clear();

  std::map<int, std::vector<CRDEcalEDM::CRDCaloHit2DShower> > m_orderedTrkShower;  m_orderedTrkShower.clear(); //map<layer, showers>
  std::map<int, std::vector<CRDEcalEDM::CRDCaloHit2DShower> > m_orderedNeuShower;  m_orderedNeuShower.clear(); //map<layer, showers>
  for(int is=0;is<m_2DshowerCol.size();is++){
    if(settings.fl_UseCandidate==0) m_orderedNeuShower[m_2DshowerCol[is].getDlayer()].push_back(m_2DshowerCol[is]);
    else{
      if(m_2DshowerCol[is].isTrkShower()) m_orderedTrkShower[m_2DshowerCol[is].getDlayer()].push_back(m_2DshowerCol[is]);
      else m_orderedNeuShower[m_2DshowerCol[is].getDlayer()].push_back(m_2DshowerCol[is]);
    }
  }

//cout<<"  Ordered NeuShower size: "<<m_orderedNeuShower.size()<<endl;
//cout<<"  Ordered TrkShower size: "<<m_orderedTrkShower.size()<<endl;

  //Longitudinal linking
  LongiConeLinking( m_orderedNeuShower, m_neuClusterCol );
  if(settings.fl_UseCandidate>0 && m_orderedTrkShower.size()>0) LongiConeLinking( m_orderedTrkShower, m_trkClusterCol ); 
//cout<<"  NeuCluster size: "<<m_neuClusterCol.size(); 
//cout<<"  TrkCluster size: "<<m_trkClusterCol.size()<<endl; 

//for(int icl=0; icl<m_neuClusterCol.size(); icl++){
//for(int a=0; a<m_neuClusterCol[icl].get2DShowers().size(); a++)
//  printf("      Shower %d: (%.1f, %.1f, %.1f, %.1f) \n", a,
//                                                         m_neuClusterCol[icl].get2DShowers()[a].getPos().x(),
//                                                         m_neuClusterCol[icl].get2DShowers()[a].getPos().y(),
//                                                         m_neuClusterCol[icl].get2DShowers()[a].getPos().z(),
//                                                         m_neuClusterCol[icl].get2DShowers()[a].getShowerE()  );
//}

  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_ClusterCol;  m_ClusterCol.clear();
  m_ClusterCol.insert(m_ClusterCol.end(), m_neuClusterCol.begin(), m_neuClusterCol.end() );
  if(m_trkClusterCol.size()>0) m_ClusterCol.insert(m_ClusterCol.end(), m_trkClusterCol.begin(), m_trkClusterCol.end() );

//cout<<"  Cluster size: "<<m_ClusterCol.size()<<endl;

  //Check cluster quality.
  std::vector< CRDEcalEDM::CRDCaloHit3DCluster >  goodClus; goodClus.clear();
  std::vector< CRDEcalEDM::CRDCaloHit3DCluster >  badClus; badClus.clear(); 
  for(int icl=0; icl<m_ClusterCol.size(); icl++){
    //if(CheckClusterQuality(m_ClusterCol[icl]) >= settings.fl_GoodClusLevel) goodClus.push_back(m_ClusterCol[icl]);
    //if( m_ClusterCol[icl].get2DShowers().size() >= settings.fl_GoodClusLevel) goodClus.push_back(m_ClusterCol[icl]);
    int m_Nlayer = m_ClusterCol[icl].getEndDlayer()-m_ClusterCol[icl].getBeginningDlayer();
    if( m_Nlayer >= settings.fl_GoodClusLevel) goodClus.push_back(m_ClusterCol[icl]);
    else badClus.push_back(m_ClusterCol[icl]);
  }

//cout<<"  Good cluster size: "<<goodClus.size()<<endl;
//cout<<"  Bad cluster size: "<<badClus.size()<<endl;

  if(settings.fl_overwrite) m_datacol.ClearCluster();
  m_datacol.GoodClus3DCol.insert(m_datacol.GoodClus3DCol.end(), goodClus.begin(), goodClus.end() );
  m_datacol.BadClus3DCol.insert( m_datacol.BadClus3DCol.end(),  badClus.begin(), badClus.end() );
  m_datacol.Clus3DCol.insert( m_datacol.Clus3DCol.end(), m_ClusterCol.begin(), m_ClusterCol.end() );

  return StatusCode::SUCCESS;
}


StatusCode ConeClusteringAlg::LongiConeLinking(
  std::map<int, std::vector<CRDEcalEDM::CRDCaloHit2DShower> >& orderedShower, 
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster>& ClusterCol)
{
  if(orderedShower.size()==0) return StatusCode::SUCCESS;
//cout<<"Start LongiConeLinking: orderedShower.size = "<<orderedShower.size()<<endl;

   std::map<int, std::vector<CRDEcalEDM::CRDCaloHit2DShower>>::iterator iter = orderedShower.begin();
  //In first layer: initial clusters. All showers in the first layer are regarded as cluster seed.
  //cluster initial direction = R.
  std::vector<CRDEcalEDM::CRDCaloHit2DShower> ShowersinFirstLayer;  ShowersinFirstLayer.clear();
  ShowersinFirstLayer = iter->second;
//cout<<" First layer: "<<iter->first<<", shower size: "<<ShowersinFirstLayer.size()<<endl;
  for(int i=0;i<ShowersinFirstLayer.size(); i++){
    CRDEcalEDM::CRDCaloHit3DCluster m_clus; m_clus.Clear(); 
    m_clus.AddShower(ShowersinFirstLayer[i]);
    ClusterCol.push_back(m_clus);
  }
  iter++;


  //Use different cone angle for 1->2/2->3 and 3->n case
  //Loop later layers
  for(iter;iter!=orderedShower.end();iter++){
    std::vector<CRDEcalEDM::CRDCaloHit2DShower> ShowersinLayer = iter->second;
    for(int is=0; is<ShowersinLayer.size(); is++){
      CRDEcalEDM::CRDCaloHit2DShower m_shower = ShowersinLayer[is];
      for(int ic=0; ic<ClusterCol.size(); ic++ ){
        int m_N2dshs = ClusterCol[ic].get2DShowers().size();
        CRDEcalEDM::CRDCaloHit2DShower shower_in_clus = ClusterCol[ic].get2DShowers().back();
        dd4hep::Position relR = m_shower.getPos()-shower_in_clus.getPos();
        TVector3 relR_vec(relR.x(), relR.y(), relR.z());
        if(  (m_N2dshs<3  && m_N2dshs>0 && relR_vec.Angle(ClusterCol[ic].getAxis())< settings.th_ConeTheta_l1 && relR_vec.Mag()< settings.th_ConeR_l1) ||
             (m_N2dshs>=3               && relR_vec.Angle(ClusterCol[ic].getAxis())< settings.th_ConeTheta_l2 && relR_vec.Mag()< settings.th_ConeR_l2)  ){

//printf("      Merged shower: (%.2f, %.2f, %.2f), Cluster last: (%.2f, %.2f, %.2f), Cluster axis: (%.2f, %.2f, %.2f) \n", 
//    m_shower.getPos().x(), m_shower.getPos().y(),m_shower.getPos().z(), 
//    shower_in_clus.getPos().x(), shower_in_clus.getPos().y(), shower_in_clus.getPos().z(),
//    ClusterCol[ic].getAxis().x(), ClusterCol[ic].getAxis().y(), ClusterCol[ic].getAxis().z() );

          ClusterCol[ic].AddShower(m_shower);
          ShowersinLayer.erase(ShowersinLayer.begin()+is);
          is--;
          break;
        }
      }
    }//end loop showers in layer.
    if(ShowersinLayer.size()>0){
      for(int i=0;i<ShowersinLayer.size(); i++){
        CRDEcalEDM::CRDCaloHit3DCluster m_clus;
        m_clus.AddShower(ShowersinLayer[i]);
        ClusterCol.push_back(m_clus);
    }}//end new cluster
  }//end loop layers.


  return StatusCode::SUCCESS;
}


int ConeClusteringAlg::CheckClusterQuality(CRDEcalEDM::CRDCaloHit3DCluster& clus){
  //Quality in different level: 
  //L0: return 0. <=1 layers(2DShowers), fluctuations in layer. 
  //L1: return 1. 2 layers(2DShower), little hadronic shower, wrong combination.
  //L2: return 2. At least 3 layers, but not pass energy ladder(E1, E2, E3). Hadronic shower/mip. 
  //L3: return 3. Pass L0 and L1 and energy ladder. EM shower/hadronic shower. 

  return clus.get2DShowers().size(); 
/*  if(clus.get2DShowers().size()<=1 ) { return 0; }

  //Criterion1: at least 3 layers. 
  if(clus.get2DShowers().size()<=3 ) { return 1; }

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
  if(!isSeed) return 2;

  return 3;
*/
}



#endif
