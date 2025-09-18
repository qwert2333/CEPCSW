#ifndef HCALCLUSTERING_ALG_C
#define HCALCLUSTERING_ALG_C

#include "Algorithm/HcalClusteringAlg.h"
using namespace Cyber;

StatusCode HcalClusteringAlg::ReadSettings(Cyber::Settings& m_settings){
  settings = m_settings;

  if(settings.map_stringVecPars.find("InputHCALHits")==settings.map_stringVecPars.end())         settings.map_stringVecPars["InputHCALHits"] = {"HCALBarrel", "HCALEndcaps"};
  if(settings.map_stringPars.find("OutputHCALClusters")==settings.map_stringPars.end())    settings.map_stringPars["OutputHCALClusters"] = "HCALCluster";

  //Set initial values
  if(settings.map_floatPars.find("th_ConeTheta_l1")==settings.map_floatPars.end())    settings.map_floatPars["th_ConeTheta_l1"] = TMath::Pi()/2.;
  if(settings.map_floatPars.find("th_ConeR_l1")==settings.map_floatPars.end())        settings.map_floatPars["th_ConeR_l1"] = 70.;
  if(settings.map_floatPars.find("th_ConeTheta_l2")==settings.map_floatPars.end())    settings.map_floatPars["th_ConeTheta_l2"] = TMath::Pi()/3.;
  if(settings.map_floatPars.find("th_ConeR_l2")==settings.map_floatPars.end())        settings.map_floatPars["th_ConeR_l2"] = 120.;
  if(settings.map_floatPars.find("th_ClusChi2")==settings.map_floatPars.end())        settings.map_floatPars["th_ClusChi2"] = 10e17;
  if(settings.map_floatPars.find("fl_GoodClusLevel")==settings.map_floatPars.end())   settings.map_floatPars["fl_GoodClusLevel"] = 10;
  if(settings.map_intPars.find("minNhit")==settings.map_intPars.end())                settings.map_intPars["minNhit"] = 3;
  if(settings.map_intPars.find("minMergeLayer")==settings.map_intPars.end())          settings.map_intPars["minMergeLayer"] = 5;
  if(settings.map_floatPars.find("maxMergeR")==settings.map_floatPars.end())          settings.map_floatPars["maxMergeR"] = 500;
  
  return StatusCode::SUCCESS;
};

StatusCode HcalClusteringAlg::Initialize( CyberDataCol& m_datacol ){
  return StatusCode::SUCCESS;
};

StatusCode HcalClusteringAlg::RunAlgorithm( CyberDataCol& m_datacol ){
  //Readin data from DataCol: 
//   std::vector<Cyber::CaloHit*> m_hcalHits;
//   m_hcalHits.clear();
//   for(int ih=0; ih<m_datacol.map_CaloHit["HCALBarrel"].size(); ih++)
//     m_hcalHits.push_back( m_datacol.map_CaloHit["HCALBarrel"][ih].get() );
//   //ordered hits by layer
//   std::map<int, std::vector<Cyber::CaloHit*> > m_orderedHit;  
//   m_orderedHit.clear();
//   for(int ih=0;ih<m_hcalHits.size();ih++)
//     m_orderedHit[m_hcalHits[ih]->getLayer()].push_back(m_hcalHits[ih]);

  std::vector<Cyber::CaloHit*> m_hcalHits;
  m_hcalHits.clear();
  for(int icl=0; icl<settings.map_stringVecPars["InputHCALHits"].size(); icl++){
    for(int ih=0; ih<m_datacol.map_CaloHit[settings.map_stringVecPars["InputHCALHits"][icl]].size(); ih++)
      m_hcalHits.push_back( m_datacol.map_CaloHit[settings.map_stringVecPars["InputHCALHits"][icl]][ih].get() );
  }

   //ordered hits by layer. 
	 //TODO: this did not seperate barrel and endcap, need to do cone-clustering separately. 
   std::map<int, std::vector<Cyber::CaloHit*> > m_orderedHit;
   m_orderedHit.clear();
   for(int ih=0;ih<m_hcalHits.size();ih++)
     m_orderedHit[m_hcalHits[ih]->getLayer()].push_back(m_hcalHits[ih]);

  std::vector<std::shared_ptr<Cyber::Calo3DCluster>> m_clusterCol;  
  m_clusterCol.clear();

  //LongiConeLinking( m_orderedHit, m_clusterCol );
  Clustering(m_hcalHits, m_clusterCol);


//   m_datacol.bk_Cluster3DCol.insert(  m_datacol.bk_Cluster3DCol.end(), m_clusterCol.begin(), m_clusterCol.end() );

  //Merge isolated hits. 
  //Do not merge hits at HCAL front 10 layers, they may need to be linked to ECAL clusters. 
  std::vector<const Cyber::CaloHit*> m_isoHits;
  for(int ic=0; ic<m_clusterCol.size(); ic++){
    if(m_clusterCol[ic]->getCaloHits().size()<=settings.map_intPars["minNhit"]){ //TODO: arbitrory criterion: cluster with <=3 hits are removed.
      int minLayer = 999;
      for(int ihit=0; ihit<m_clusterCol[ic]->getCaloHits().size(); ihit++){
        if(m_clusterCol[ic]->getCaloHits()[ihit]->getLayer()<minLayer) 
          minLayer = m_clusterCol[ic]->getCaloHits()[ihit]->getLayer();
      }

      if(minLayer<=settings.map_intPars["minMergeLayer"]) continue;

      for(int ihit=0; ihit<m_clusterCol[ic]->getCaloHits().size(); ihit++)      
        m_isoHits.push_back(m_clusterCol[ic]->getCaloHits()[ihit]);

      m_clusterCol.erase(m_clusterCol.begin()+ic);      
    }
  }

  MergeToCluster(m_isoHits, m_clusterCol);
 
  //cout<<"  After HCAL clustering: Cluster size "<<m_clusterCol.size()<<endl;
  //int Nisohits = 0;
  //for(int ic=0; ic<m_clusterCol.size(); ic++)
  //{
  //  cout<<"    Cluster #"<<ic<<": hit size "<<m_clusterCol[ic]->getCaloHits().size()<<endl;
  //  if(m_clusterCol[ic]->getCaloHits().size()==1) Nisohits++;
  //}
  //cout<<"  Isolated hit size: "<<Nisohits<<endl;
 
//   m_datacol.map_CaloCluster[settings.map_stringPars["OutputCluster"]] = m_clusterCol;
  m_datacol.map_CaloCluster[settings.map_stringPars["OutputHCALClusters"]]= m_clusterCol;
  return StatusCode::SUCCESS;
};

StatusCode HcalClusteringAlg::ClearAlgorithm(){

  //Clear local memory
//   m_hcalHits.clear();

  return StatusCode::SUCCESS;
};

StatusCode HcalClusteringAlg::LongiConeLinking(  const std::map<int, std::vector<Cyber::CaloHit*> >& orderedHit, 
                                                 std::vector<std::shared_ptr<Cyber::Calo3DCluster> >& ClusterCol)
{

  if(orderedHit.size()==0) return StatusCode::SUCCESS;

  auto iter = orderedHit.begin();
  //In first layer: initial clusters. All showers in the first layer are regarded as cluster seed.
  //cluster initial direction = R.
  std::vector<Cyber::CaloHit*> HitsinFirstLayer = iter->second;
  for(int i=0;i<HitsinFirstLayer.size(); i++){
    std::shared_ptr<Cyber::Calo3DCluster> m_clus = std::make_shared<Cyber::Calo3DCluster>();
    m_clus->addHit(HitsinFirstLayer[i]);
    ClusterCol.push_back(m_clus);
  }
  iter++;

//cout<<"    LongiConeLinking: Cluster seed in first layer: "<<ClusterCol.size()<<endl;

  //Use different cone angle for 1->2/2->3 and 3->n case
  //Loop later layers
  for(iter; iter!=orderedHit.end(); iter++){
    std::vector<Cyber::CaloHit*> HitsinLayer = iter->second;
//cout<<"    In Layer: "<<iter->first<<"  Hit size: "<<HitsinLayer.size()<<endl;

    for(int is=0; is<HitsinLayer.size(); is++){
      Cyber::CaloHit* m_hit = HitsinLayer[is];
//printf("     New Hit: (%.3f, %.3f, %.3f), Layer %d \n", m_hit->getPosition().x(), m_hit->getPosition().y(), m_hit->getPosition().z(),  m_hit->getLayer() );
//cout<<"     Cluster size: "<<ClusterCol.size()<<endl;

      for(int ic=0; ic<ClusterCol.size(); ic++ ){
        int m_Nhits = ClusterCol[ic]->getCaloHits().size();
        const Cyber::CaloHit* hit_in_clus = ClusterCol[ic]->getCaloHits().back();
        TVector3 relR_vec = m_hit->getPosition() - hit_in_clus->getPosition();
//printf("      New hit: (%.2f, %.2f, %.2f), Cluster last: (%.2f, %.2f, %.2f, %d), Cluster axis: (%.2f, %.2f, %.2f) \n",
//    m_hit->getPosition().x(), m_hit->getPosition().y(),m_hit->getPosition().z(),
//    hit_in_clus->getPosition().x(), hit_in_clus->getPosition().y(), hit_in_clus->getPosition().z(), hit_in_clus->getLayer(), 
//    ClusterCol[ic]->getAxis().x(), ClusterCol[ic]->getAxis().y(), ClusterCol[ic]->getAxis().z() );

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
        std::shared_ptr<Cyber::Calo3DCluster> m_clus = std::make_shared<Cyber::Calo3DCluster>();
        m_clus->addHit(HitsinLayer[i]);
        ClusterCol.push_back(m_clus);
    }}//end new cluster
  }//end loop layers.


  return StatusCode::SUCCESS;
}

template<typename T1, typename T2> StatusCode HcalClusteringAlg::Clustering(std::vector<T1*> &m_input, std::vector<std::shared_ptr<T2>> &m_output) 
{
  std::vector<std::shared_ptr<T2>> record;
  record.clear();

  for(int i=0; i<m_input.size(); i++)
  {
    T1* lowlevelcluster = m_input.at(i);
    for(int j=0; j<m_output.size(); j++)
    {
      if(m_output.at(j).get()->isHCALNeighbor(lowlevelcluster)) // //m_output.at(j).isNeighbor(lowlevelcluster)
        record.push_back(m_output.at(j));

    }
    if(record.size()>0)
    {
      record.at(0).get()->addHit(lowlevelcluster);
      for(int k=1; k<record.size(); k++)
      {
        for(int l=0; l<record.at(k).get()->getCaloHits().size(); l++)
        {
          record.at(0).get()->addHit(record.at(k).get()->getCaloHits().at(l));
        }
      }
      for(int m=1; m<record.size(); m++)
      {
        m_output.erase(find(m_output.begin(),m_output.end(),record.at(m)));
        //delete record.at(m); 
        record.at(m) = nullptr;
      }
      record.clear();
      lowlevelcluster = nullptr;
      continue;
    }
    //T2* highlevelcluster = new T2(); //first new
    std::shared_ptr<T2> highlevelcluster = std::make_shared<T2>();

    highlevelcluster.get()->addHit(lowlevelcluster);
    m_output.push_back(highlevelcluster);
  }
  return StatusCode::SUCCESS;
	//cout<<"how many neighbors: "<<number<<endl;
}

StatusCode HcalClusteringAlg::MergeToCluster( std::vector<const Cyber::CaloHit*>& m_isohits, std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_clusters ){

  if(m_isohits.size()==0) return StatusCode::SUCCESS;

  for(int ihit=0; ihit<m_isohits.size(); ihit++){
    int index_closest = -1;
    double minDistance = 1e6;
    for(int icl=0; icl<m_clusters.size(); icl++){
      for(auto jhit : m_clusters[icl]->getCaloHits()){
        double tmp_R = sqrt( pow(jhit->getPosition().x()-m_isohits[ihit]->getPosition().x(),2) + pow(jhit->getPosition().y()-m_isohits[ihit]->getPosition().y(),2) + pow(jhit->getPosition().z()-m_isohits[ihit]->getPosition().z(),2) );
        if( tmp_R>settings.map_floatPars["maxMergeR"] && tmp_R<minDistance)
          { minDistance = tmp_R; index_closest = icl; }
      }
    }

    if(index_closest>=0){
      m_clusters[index_closest]->addHit( m_isohits[ihit] );
      m_isohits.erase(m_isohits.begin()+ihit);
      ihit--;
    }  
    else{
      std::shared_ptr<Cyber::Calo3DCluster> m_clus = std::make_shared<Cyber::Calo3DCluster>();
      m_clus->addHit(m_isohits[ihit]);
      m_clusters.push_back(m_clus);
      m_isohits.erase(m_isohits.begin()+ihit);
      ihit--;
    }

  }

  return StatusCode::SUCCESS;
}

#endif
