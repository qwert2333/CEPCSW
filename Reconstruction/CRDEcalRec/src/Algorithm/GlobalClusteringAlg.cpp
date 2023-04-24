#ifndef GLOBALCLUSTERING_ALG_C
#define GLOBALCLUSTERING_ALG_C

#include "Algorithm/GlobalClusteringAlg.h"
using namespace PandoraPlus;

StatusCode GlobalClusteringAlg::ReadSettings(PandoraPlus::Settings& m_settings){
  settings = m_settings;
  if(settings.map_floatPars.find("unit_threshold")==settings.map_floatPars.end())  settings.map_floatPars["unit_threshold"] = 0.001;

  return StatusCode::SUCCESS;
};

StatusCode GlobalClusteringAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  m_bars.clear();
  m_processbars.clear(); 
  m_restbars.clear();
  m_1dclusters.clear();
  m_halfclusters.clear();

  //Readin data from DataCol: 
  m_bars = m_datacol.map_BarCol["BarCol"];


  return StatusCode::SUCCESS;
};

StatusCode GlobalClusteringAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){
  //system("/cefs/higgs/songwz/winter22/CEPCSW/workarea/memory/memory_test.sh cluster_begin");
  //time_t time_cb;  
  //time(&time_cb);  
  //cout<<" When begin clustering: "<<ctime(&time_cb)<<endl;
cout<<"  GlobalClusteringAlg: Readin Bar size: "<<m_bars.size()<<endl;

  //Threshold and scale factor (Todo) for bars. 
  for(int ibar=0; ibar<m_bars.size(); ibar++)
  {
    if(m_bars.at(ibar)->getEnergy()>settings.map_floatPars["unit_threshold"])
    {
      m_processbars.push_back(m_bars.at(ibar));
    }
    else
    {
      m_restbars.push_back(m_bars.at(ibar));
    }
  }

  //Clustering
  Clustering(m_processbars, m_1dclusters);
  Clustering(m_1dclusters, m_halfclusters);

  //Store created objects to backup col. 
  //for(int irest=0; irest<m_restbars.size(); irest++)       m_datacol.bk_BarCol.push_back(m_restbars.at(irest));
  //for(int i1d=0; i1d<m_1dclusters.size(); i1d++)           m_datacol.bk_Cluster1DCol.push_back(m_1dclusters.at(i1d));
  //for(int ihalf=0; ihalf<m_halfclusters.size(); ihalf++)   m_datacol.bk_ClusterHalfCol.push_back(m_halfclusters.at(ihalf));


  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_HalfClusterV; m_HalfClusterV.clear();
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_HalfClusterU; m_HalfClusterU.clear();
  for(int i=0; i<m_halfclusters.size() && m_halfclusters[i]; i++){
    if(m_halfclusters[i]->getSlayer()==0)
      m_HalfClusterU.push_back(m_halfclusters[i]);
    else if(m_halfclusters[i]->getSlayer()==1)
      m_HalfClusterV.push_back(m_halfclusters[i]);
  }

printf("  GlobalClustering: RestBarCol size %d, 1DCluster size %d, HalfCluster size (%d, %d) \n", m_restbars.size(), m_1dclusters.size(), m_HalfClusterU.size(), m_HalfClusterV.size());
//for(int ic=0; ic<m_HalfClusterU.size(); ic++) cout<<m_HalfClusterU[ic]->getEnergy()<<'\t';
//cout<<endl;
//for(int ic=0; ic<m_HalfClusterV.size(); ic++) cout<<m_HalfClusterV[ic]->getEnergy()<<'\t';
//cout<<endl;


  //Write results into DataCol.
  //m_datacol.map_BarCol["RestBarCol"] = m_restbars;
  m_datacol.map_1DCluster["Cluster1DCol"] = m_1dclusters;
  m_datacol.map_HalfCluster["HalfClusterColU"] = m_HalfClusterU;
  m_datacol.map_HalfCluster["HalfClusterColV"] = m_HalfClusterV;

  //time_t time_ce;  
  //time(&time_ce); 
  //cout<<" When end clustering: "<<ctime(&time_ce)<<endl;
  //system("/cefs/higgs/songwz/winter22/CEPCSW/workarea/memory/memory_test.sh cluster_end");
  return StatusCode::SUCCESS;
};

StatusCode GlobalClusteringAlg::ClearAlgorithm(){

  //Clear local memory
  m_bars.clear();
  m_processbars.clear();
  m_restbars.clear();
  m_1dclusters.clear();
  m_halfclusters.clear();

  return StatusCode::SUCCESS;
};

template<typename T1, typename T2> StatusCode GlobalClusteringAlg::Clustering(std::vector<std::shared_ptr<T1>> &m_input, std::vector<std::shared_ptr<T2>> &m_output) 
{
  std::vector<std::shared_ptr<T2>> record;
  record.clear();

  for(int i=0; i<m_input.size(); i++)
  {
    T1* lowlevelcluster = m_input.at(i).get();
    for(int j=0; j<m_output.size(); j++)
    {
      if(m_output.at(j).get()->isNeighbor(lowlevelcluster)) // //m_output.at(j).isNeighbor(lowlevelcluster)
        record.push_back(m_output.at(j));

    }
    if(record.size()>0)
    {
      record.at(0).get()->addUnit(lowlevelcluster);
      for(int k=1; k<record.size(); k++)
      {
        for(int l=0; l<record.at(k).get()->getCluster().size(); l++)
        {
          record.at(0).get()->addUnit(record.at(k).get()->getCluster().at(l));
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

    highlevelcluster.get()->addUnit(lowlevelcluster);
    m_output.push_back(highlevelcluster);
  }
	return StatusCode::SUCCESS;
	//cout<<"how many neighbors: "<<number<<endl;
}


#endif
