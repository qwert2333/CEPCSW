#ifndef _CONECLUSTERING_ALG_C
#define _CONECLUSTERING_ALG_C

#include "Algorithm/ConeClusteringAlg.h"

StatusCode ConeClusteringAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Option: readin information
  if(settings.map_stringPars.find("ReadinHit")==settings.map_stringPars.end())        settings.map_stringPars["ReadinHit"] = "EcalTransShower";
  if(settings.map_stringPars.find("OutputCluster")==settings.map_stringPars.end())    settings.map_stringPars["OutputCluster"] = "EcalBarCluster";


  //Set initial values
  if(settings.map_floatPars.find("th_ConeTheta_l1")==settings.map_floatPars.end())    settings.map_floatPars["th_ConeTheta_l1"] = TMath::Pi()/4.;
  if(settings.map_floatPars.find("th_ConeR_l1")==settings.map_floatPars.end())        settings.map_floatPars["th_ConeR_l1"] = 30.;
  if(settings.map_floatPars.find("th_ConeTheta_l2")==settings.map_floatPars.end())    settings.map_floatPars["th_ConeTheta_l2"] = TMath::Pi()/6.;
  if(settings.map_floatPars.find("th_ConeR_l2")==settings.map_floatPars.end())        settings.map_floatPars["th_ConeR_l2"] = 30.;
  if(settings.map_floatPars.find("th_ClusChi2")==settings.map_floatPars.end())        settings.map_floatPars["th_ClusChi2"] = 10e17;
  if(settings.map_floatPars.find("fl_GoodClusLevel")==settings.map_floatPars.end())   settings.map_floatPars["fl_GoodClusLevel"] = 4;
  //For Cluster merging 
  if(settings.map_floatPars.find("axis_Angle")==settings.map_floatPars.end())         settings.map_floatPars["axis_Angle"] = TMath::Pi()/6.;
  if(settings.map_floatPars.find("relP_Angle")==settings.map_floatPars.end())         settings.map_floatPars["relP_Angle"] = TMath::Pi()/5.;
  if(settings.map_floatPars.find("skipLayer")==settings.map_floatPars.end())          settings.map_floatPars["skipLayer"] = 3;
  if(settings.map_floatPars.find("fl_MergeGoodClus")==settings.map_floatPars.end())   settings.map_floatPars["fl_MergeGoodClus"] = 1;
  if(settings.map_floatPars.find("fl_MergeBadClus")==settings.map_floatPars.end())    settings.map_floatPars["fl_MergeBadClus"] = 1;
  if(settings.map_floatPars.find("fl_MergeEMTail")==settings.map_floatPars.end())     settings.map_floatPars["fl_MergeEMTail"] = 1;

  return StatusCode::SUCCESS;
};

StatusCode ConeClusteringAlg::Initialize(){

  return StatusCode::SUCCESS;
}


StatusCode ConeClusteringAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol){
cout<<"  ConeClustering: Readin "<<settings.map_stringPars["ReadinHit"]<<", Output: "<<settings.map_stringPars["OutputCluster"]<<endl;

  std::vector<PandoraPlus::CaloHit*> m_hitCol; m_hitCol.clear(); 
  //Get readin info
  if(settings.map_stringPars["ReadinHit"] == "EcalTransShower"){

    std::vector<PandoraPlus::TransShower*>* p_Tshowers = &(m_datacol.TransShowerCol);
    if(p_Tshowers->size()==0){ 
      std::cout<<"Warning: Empty input in ConeClusteringAlg. Please check previous algorithm!"<<endl;  
      return StatusCode::SUCCESS; 
    }
    //Convert transShower to hit
    for(int is=0; is<p_Tshowers->size(); is++){
      PandoraPlus::CaloHit* m_hit = new PandoraPlus::CaloHit();
      m_hit->setEnergy( p_Tshowers->at(is)->getShowerE() );
      m_hit->setPosition( p_Tshowers->at(is)->getPos() );
      m_hit->setLayer( p_Tshowers->at(is)->getDlayer() );
      //m_hit->setParentShower( p_Tshowers->at(is) );
      m_hitCol.push_back(m_hit);
      m_datacol.bk_HitCol.push_back(m_hit);
    }
    p_Tshowers = nullptr;
  }
  else m_hitCol = m_datacol.map_CaloHit[settings.map_stringPars["ReadinHit"]]; 
  if(m_hitCol.size()==0) { std::cout<<"Warning: Empty input in ConeClusteringAlg. Please check previous algorithm!"<<endl; return StatusCode::SUCCESS; }



cout<<"ConeClustering: Print Input Hits"<<endl;
for(int i=0; i<m_hitCol.size(); i++)
  printf("  Hit #%d: pos/E (%.2f, %.2f, %.2f, %.3f), Layer %d. \n", i,  
        m_hitCol[i]->getPosition().x(), m_hitCol[i]->getPosition().y(), m_hitCol[i]->getPosition().z(), m_hitCol[i]->getEnergy(), m_hitCol[i]->getLayer() );
cout<<"ConeClustering: End Print"<<endl;







  //Store the ordered hits.
  std::map<int, std::vector<PandoraPlus::CaloHit*> > m_orderedHit;  m_orderedHit.clear(); //map<layer, showers>
  for(int ih=0;ih<m_hitCol.size();ih++)
    m_orderedHit[m_hitCol[ih]->getLayer()].push_back(m_hitCol[ih]);


  //Longitudinal linking
  std::vector<PandoraPlus::CaloCluster*> m_clusterCol;  m_clusterCol.clear();
  LongiConeLinking( m_orderedHit, m_clusterCol );
  m_datacol.bk_ClusterCol.insert(  m_datacol.bk_ClusterCol.end(), m_clusterCol.begin(), m_clusterCol.end() );
cout<<"  Cluster size: "<<m_clusterCol.size()<<endl;


  //Check cluster quality.
  std::vector<PandoraPlus::CaloCluster*>  goodClus;
  std::vector<PandoraPlus::CaloCluster*>  badClus;
  for(int icl=0; icl<m_clusterCol.size(); icl++){
    if( m_clusterCol[icl]->getCaloHits().size() >= settings.map_floatPars["fl_GoodClusLevel"]) goodClus.push_back(m_clusterCol[icl]);
    else badClus.push_back(m_clusterCol[icl]);
  }

cout<<"  Good cluster: "<<goodClus.size()<<endl;
cout<<"  Bad cluster: "<<badClus.size()<<endl;

  //Merge clusters
  //MergeGoodClusters( goodClus );
  //for(int icl=0;icl<badClus.size();icl++) MergeBadToGoodCluster(goodClus, badClus[icl] ); 
  //MergeGoodClusters( goodClus );
  //MergeEMTail( goodClus );

  //m_datacol.GoodClus3DCol.insert(m_datacol.GoodClus3DCol.end(), goodClus.begin(), goodClus.end() );
  //m_datacol.BadClus3DCol.insert( m_datacol.BadClus3DCol.end(),  badClus.begin(), badClus.end() );
  //m_datacol.Clus3DCol.insert( m_datacol.Clus3DCol.end(), m_clusterCol.begin(), m_clusterCol.end() );

  m_datacol.map_CaloCluster[settings.map_stringPars["OutputCluster"]] = m_clusterCol;

  return StatusCode::SUCCESS;
}

StatusCode ConeClusteringAlg::ClearAlgorithm(){
  settings.Clear();

  return StatusCode::SUCCESS;
}

StatusCode ConeClusteringAlg::LongiConeLinking(  const std::map<int, std::vector<PandoraPlus::CaloHit*> >& orderedHit, 
                                                 std::vector<PandoraPlus::CaloCluster*>& ClusterCol)
{

  if(orderedHit.size()==0) return StatusCode::SUCCESS;

  auto iter = orderedHit.begin();
  //In first layer: initial clusters. All showers in the first layer are regarded as cluster seed.
  //cluster initial direction = R.
  std::vector<PandoraPlus::CaloHit*> HitsinFirstLayer = iter->second;
  for(int i=0;i<HitsinFirstLayer.size(); i++){
    PandoraPlus::CaloCluster* m_clus = new PandoraPlus::CaloCluster();
    m_clus->addHit(HitsinFirstLayer[i]);
    ClusterCol.push_back(m_clus);
  }
  iter++;

cout<<"    LongiConeLinking: Cluster seed in first layer: "<<ClusterCol.size()<<endl;

  //Use different cone angle for 1->2/2->3 and 3->n case
  //Loop later layers
  for(iter; iter!=orderedHit.end(); iter++){
    std::vector<PandoraPlus::CaloHit*> HitsinLayer = iter->second;
cout<<"    In Layer: "<<iter->first<<"  Hit size: "<<HitsinLayer.size()<<endl;

    for(int is=0; is<HitsinLayer.size(); is++){
      PandoraPlus::CaloHit* m_hit = HitsinLayer[is];
printf("     New Hit: (%.3f, %.3f, %.3f), Layer %d \n", m_hit->getPosition().x(), m_hit->getPosition().y(), m_hit->getPosition().z(),  m_hit->getLayer() );
cout<<"     Cluster size: "<<ClusterCol.size()<<endl;

      for(int ic=0; ic<ClusterCol.size(); ic++ ){
        int m_Nhits = ClusterCol[ic]->getCaloHits().size();
        const PandoraPlus::CaloHit* hit_in_clus = ClusterCol[ic]->getCaloHits().back();
        TVector3 relR_vec = m_hit->getPosition() - hit_in_clus->getPosition();
printf("      New hit: (%.2f, %.2f, %.2f), Cluster last: (%.2f, %.2f, %.2f, %d), Cluster axis: (%.2f, %.2f, %.2f) \n",
    m_hit->getPosition().x(), m_hit->getPosition().y(),m_hit->getPosition().z(),
    hit_in_clus->getPosition().x(), hit_in_clus->getPosition().y(), hit_in_clus->getPosition().z(), hit_in_clus->getLayer(), 
    ClusterCol[ic]->getAxis().x(), ClusterCol[ic]->getAxis().y(), ClusterCol[ic]->getAxis().z() );

        if(  (m_Nhits<3 && m_Nhits>0 && relR_vec.Angle(ClusterCol[ic]->getAxis())< settings.map_floatPars["th_ConeTheta_l1"] && relR_vec.Mag()< settings.map_floatPars["th_ConeR_l1"]) ||
             (m_Nhits>=3             && relR_vec.Angle(ClusterCol[ic]->getAxis())< settings.map_floatPars["th_ConeTheta_l2"] && relR_vec.Mag()< settings.map_floatPars["th_ConeR_l2"])  ){

          ClusterCol[ic]->addHit(m_hit);
          HitsinLayer.erase(HitsinLayer.begin()+is);
          is--;
          break;
        }
      }
    }//end loop showers in layer.
    if(HitsinLayer.size()>0){
      for(int i=0;i<HitsinLayer.size(); i++){
        PandoraPlus::CaloCluster* m_clus = new PandoraPlus::CaloCluster();
        m_clus->addHit(HitsinLayer[i]);
        ClusterCol.push_back(m_clus);
    }}//end new cluster
  }//end loop layers.


  return StatusCode::SUCCESS;
}

StatusCode ConeClusteringAlg::MergeGoodClusters( std::vector<PandoraPlus::CaloCluster*>& m_clusCol ){

  return StatusCode::SUCCESS;
}

StatusCode ConeClusteringAlg::MergeBadToGoodCluster( std::vector<PandoraPlus::CaloCluster*>& m_goodClusCol, PandoraPlus::CaloCluster* m_badClus ){
  PandoraPlus::CaloCluster* m_clus = GetClosestGoodCluster( m_goodClusCol, m_badClus );
  if(!m_clus) return StatusCode::FAILURE;

  auto iter = find( m_goodClusCol.begin(), m_goodClusCol.end(), m_clus );
  if(iter==m_goodClusCol.end()) return StatusCode::FAILURE;
  else m_clus->MergeCluster( m_badClus );
    
  return StatusCode::SUCCESS;
}



PandoraPlus::CaloCluster* ConeClusteringAlg::GetClosestGoodCluster( std::vector< PandoraPlus::CaloCluster* >& m_goodClusCol, PandoraPlus::CaloCluster* m_badClus ){

  TVector3 m_clusCent = m_badClus->getShowerCenter();
  PandoraPlus::CaloCluster* m_clus = nullptr;
  double minTheta=999;
  for(int i=0;i<m_goodClusCol.size();i++){
     TVector3 vec_goodClus = m_goodClusCol[i]->getShowerCenter();
     TVector3 vec_goodAxis = m_goodClusCol[i]->getAxis();
     double theta = vec_goodAxis.Cross(m_clusCent-vec_goodClus).Mag();

     if(theta<minTheta){
        minTheta=theta;
        m_clus = m_goodClusCol[i];
     }
  }

  return m_clus;
}

#endif
