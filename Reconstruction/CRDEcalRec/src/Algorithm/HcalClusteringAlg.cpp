#ifndef HCALCLUSTERING_ALG_C
#define HCALCLUSTERING_ALG_C

#include "Algorithm/HcalClusteringAlg.h"
using namespace PandoraPlus;

StatusCode HcalClusteringAlg::ReadSettings(PandoraPlus::Settings& m_settings){
  settings = m_settings;

  if(settings.map_stringPars.find("HCALCluster")==settings.map_stringPars.end())    settings.map_stringPars["HCALCluster"] = "HCALCluster";

  //Set initial values
  if(settings.map_floatPars.find("th_ConeTheta_l1")==settings.map_floatPars.end())    settings.map_floatPars["th_ConeTheta_l1"] = TMath::Pi()/2.;
  if(settings.map_floatPars.find("th_ConeR_l1")==settings.map_floatPars.end())        settings.map_floatPars["th_ConeR_l1"] = 70.;
  if(settings.map_floatPars.find("th_ConeTheta_l2")==settings.map_floatPars.end())    settings.map_floatPars["th_ConeTheta_l2"] = TMath::Pi()/3.;
  if(settings.map_floatPars.find("th_ConeR_l2")==settings.map_floatPars.end())        settings.map_floatPars["th_ConeR_l2"] = 120.;
  if(settings.map_floatPars.find("th_ClusChi2")==settings.map_floatPars.end())        settings.map_floatPars["th_ClusChi2"] = 10e17;
  if(settings.map_floatPars.find("fl_GoodClusLevel")==settings.map_floatPars.end())   settings.map_floatPars["fl_GoodClusLevel"] = 10;

  return StatusCode::SUCCESS;
};

StatusCode HcalClusteringAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  return StatusCode::SUCCESS;
};

StatusCode HcalClusteringAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){
  //Readin data from DataCol: 
//   std::vector<PandoraPlus::CaloHit*> m_hcalHits;
//   m_hcalHits.clear();
//   for(int ih=0; ih<m_datacol.map_CaloHit["HCALBarrel"].size(); ih++)
//     m_hcalHits.push_back( m_datacol.map_CaloHit["HCALBarrel"][ih].get() );
//   //ordered hits by layer
//   std::map<int, std::vector<PandoraPlus::CaloHit*> > m_orderedHit;  
//   m_orderedHit.clear();
//   for(int ih=0;ih<m_hcalHits.size();ih++)
//     m_orderedHit[m_hcalHits[ih]->getLayer()].push_back(m_hcalHits[ih]);

  std::vector<PandoraPlus::CaloHit*> m_hcalHits;
  m_hcalHits.clear();
  for(int ih=0; ih<m_datacol.map_CaloHit["HCALBarrel"].size(); ih++)
    m_hcalHits.push_back( m_datacol.map_CaloHit["HCALBarrel"][ih].get() );

  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> m_clusterCol;  
  m_clusterCol.clear();

//   LongiConeLinking( m_orderedHit, m_clusterCol );
  Clustering(m_hcalHits, m_clusterCol);
//   m_datacol.bk_Cluster3DCol.insert(  m_datacol.bk_Cluster3DCol.end(), m_clusterCol.begin(), m_clusterCol.end() );
  // cout<<"  Cluster size: "<<m_clusterCol.size()<<endl;
  // for(int ic=0; ic<m_clusterCol.size(); ic++)
  // {
  //   cout<<"    Cluster "<<ic<<":"<<m_clusterCol[ic]->getCaloHits().size()<<endl;
  // }

//   m_datacol.map_CaloCluster[settings.map_stringPars["OutputCluster"]] = m_clusterCol;
  m_datacol.map_CaloCluster["HCALCluster"]= m_clusterCol;
  return StatusCode::SUCCESS;
};

StatusCode HcalClusteringAlg::ClearAlgorithm(){

  //Clear local memory
//   m_hcalHits.clear();

  return StatusCode::SUCCESS;
};

StatusCode HcalClusteringAlg::LongiConeLinking(  const std::map<int, std::vector<PandoraPlus::CaloHit*> >& orderedHit, 
                                                 std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster> >& ClusterCol)
{

  if(orderedHit.size()==0) return StatusCode::SUCCESS;

  auto iter = orderedHit.begin();
  //In first layer: initial clusters. All showers in the first layer are regarded as cluster seed.
  //cluster initial direction = R.
  std::vector<PandoraPlus::CaloHit*> HitsinFirstLayer = iter->second;
  for(int i=0;i<HitsinFirstLayer.size(); i++){
    std::shared_ptr<PandoraPlus::Calo3DCluster> m_clus = std::make_shared<PandoraPlus::Calo3DCluster>();
    m_clus->addHit(HitsinFirstLayer[i]);
    ClusterCol.push_back(m_clus);
  }
  iter++;

//cout<<"    LongiConeLinking: Cluster seed in first layer: "<<ClusterCol.size()<<endl;

  //Use different cone angle for 1->2/2->3 and 3->n case
  //Loop later layers
  for(iter; iter!=orderedHit.end(); iter++){
    std::vector<PandoraPlus::CaloHit*> HitsinLayer = iter->second;
//cout<<"    In Layer: "<<iter->first<<"  Hit size: "<<HitsinLayer.size()<<endl;

    for(int is=0; is<HitsinLayer.size(); is++){
      PandoraPlus::CaloHit* m_hit = HitsinLayer[is];
//printf("     New Hit: (%.3f, %.3f, %.3f), Layer %d \n", m_hit->getPosition().x(), m_hit->getPosition().y(), m_hit->getPosition().z(),  m_hit->getLayer() );
//cout<<"     Cluster size: "<<ClusterCol.size()<<endl;

      for(int ic=0; ic<ClusterCol.size(); ic++ ){
        int m_Nhits = ClusterCol[ic]->getCaloHits().size();
        const PandoraPlus::CaloHit* hit_in_clus = ClusterCol[ic]->getCaloHits().back();
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
        std::shared_ptr<PandoraPlus::Calo3DCluster> m_clus = std::make_shared<PandoraPlus::Calo3DCluster>();
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

#endif
