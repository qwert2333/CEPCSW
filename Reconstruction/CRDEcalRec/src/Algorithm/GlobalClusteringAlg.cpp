#ifndef GLOBALCLUSTERING_ALG_C
#define GLOBALCLUSTERING_ALG_C

#include "Algorithm/GlobalClusteringAlg.h"
using namespace PandoraPlus;

StatusCode GlobalClusteringAlg::ReadSettings(PandoraPlus::Settings& m_settings){
  settings = m_settings;
  if(settings.map_floatPars.find("unit_threshold")==settings.map_floatPars.end())  settings.map_floatPars["unit_threshold"] = 0.001;

  return StatusCode::SUCCESS;
};

StatusCode GlobalClusteringAlg::Initialize(){

  return StatusCode::SUCCESS;
};

StatusCode GlobalClusteringAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){
  system("/cefs/higgs/songwz/winter22/CEPCSW/workarea/memory/memory_test.sh cluster_begin");
  time_t time_cb;  
  time(&time_cb);  
  cout<<" When begin clustering: "<<ctime(&time_cb)<<endl;

  std::vector<PandoraPlus::CaloUnit*> m_bars = m_datacol.BarCol; 
  std::vector<PandoraPlus::CaloUnit*> m_processbars;         m_processbars.clear();
  std::vector<PandoraPlus::CaloUnit*> m_restbars;            m_restbars.clear();
  std::vector<PandoraPlus::Calo1DCluster*> m_1dclusters;     m_1dclusters.clear();
  std::vector<PandoraPlus::Calo2DCluster*> m_2dclusters;     m_2dclusters.clear();
  std::vector<PandoraPlus::CaloHalfCluster*> m_halfclusters; m_halfclusters.clear();
  std::vector<PandoraPlus::Calo3DCluster*> m_3dclusters;     m_3dclusters.clear();

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
  cout<<"  How many bars: "<<m_bars.size()<<endl;
  cout<<"  How many bars over threshold: "<<m_processbars.size()<<endl;
  cout<<"  Clustering bars to 1DClusters: "<<endl;
  Clustering(m_processbars, m_1dclusters);
  cout<<"  1DCluster size: "<<m_1dclusters.size()<<".  Clustering 1DClusters to HalfClusters: "<<endl;
  Clustering(m_1dclusters, m_halfclusters);
  cout<<"  HalfCluster size: "<<m_halfclusters.size()<<endl;
  // Clustering(m_2dclusters, m_3dclusters);
  // cout<<"  3DCluster size: "<<m_3dclusters.size()<<endl;
  


/*
cout<<endl;
cout<<"  Check 3DClusters"<<endl;
for(int i3d=0; i3d<m_3dclusters.size(); i3d++){
  printf("    3DClus #%d: energy %.5f, 2DClus size %d, tower size %d, towerID: ", i3d, m_3dclusters[i3d]->getEnergy(), m_3dclusters[i3d]->getCluster().size(), m_3dclusters[i3d]->getTowerID().size() );
  for(int it=0; it<m_3dclusters[i3d]->getTowerID().size(); it++) printf("[%d, %d, %d], ", m_3dclusters[i3d]->getTowerID()[it][0],  m_3dclusters[i3d]->getTowerID()[it][1], m_3dclusters[i3d]->getTowerID()[it][2] );
  cout<<endl;
  cout<<"    Check 2DClus in this 3D: "<<endl;
  for(int i2d=0; i2d<m_3dclusters[i3d]->getCluster().size(); i2d++){
    const Calo2DCluster* p_clus = m_3dclusters[i3d]->getCluster()[i2d];
    printf("      2DClus #%d: Layer %d, energy %.5f, 1DClus size (%d, %d), tower size %d, towerID: ", i2d, p_clus->getDlayer(), p_clus->getEnergy(), p_clus->getClusterU().size(), p_clus->getClusterV().size(), p_clus->getTowerID().size() );
    for(int it=0; it<p_clus->getTowerID().size(); it++) printf("[%d, %d, %d],", p_clus->getTowerID()[it][0], p_clus->getTowerID()[it][1], p_clus->getTowerID()[it][2]);
    cout<<endl;
    p_clus = nullptr;
  }
cout<<endl;
}
cout<<endl;
*/

  for(int irest=0; irest<m_restbars.size(); irest++) m_datacol.bk_RestBarCol.push_back(m_restbars.at(irest));
  for(int i1d=0; i1d<m_1dclusters.size(); i1d++) m_datacol.bk_Cluster1DCol.push_back(m_1dclusters.at(i1d));
  for(int i2d=0; i2d<m_2dclusters.size(); i2d++) m_datacol.bk_Cluster2DCol.push_back(m_2dclusters.at(i2d));
  for(int ihalf=0; ihalf<m_halfclusters.size(); ihalf++) m_datacol.bk_ClusterHalfCol.push_back(m_halfclusters.at(ihalf));
  for(int i3d=0; i3d<m_3dclusters.size(); i3d++) m_datacol.bk_Cluster3DCol.push_back(m_3dclusters.at(i3d));

  m_datacol.RestBarCol = m_restbars;
	m_datacol.Cluster1DCol = m_1dclusters;
	m_datacol.Cluster2DCol = m_2dclusters;
  m_datacol.ClusterHalfCol = m_halfclusters;
	m_datacol.Cluster3DCol = m_3dclusters;

  time_t time_ce;  
  time(&time_ce); 
  cout<<" When end clustering: "<<ctime(&time_ce)<<endl;
  system("/cefs/higgs/songwz/winter22/CEPCSW/workarea/memory/memory_test.sh cluster_end");
  return StatusCode::SUCCESS;
};

StatusCode GlobalClusteringAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
};

template<typename T1, typename T2> StatusCode GlobalClusteringAlg::Clustering(std::vector<T1*> &m_input, std::vector<T2*> &m_output) 
{
  std::vector<T2*> record;
  record.clear();
  for(int i=0; i<m_input.size(); i++)
  {
    T1* lowlevelcluster = m_input.at(i);
    for(int j=0; j<m_output.size(); j++)
    {
      if(m_output.at(j)->isNeighbor(lowlevelcluster)) // //m_output.at(j).isNeighbor(lowlevelcluster)
        record.push_back(m_output.at(j));

    }
    if(record.size()>0)
    {
      record.at(0)->addUnit(lowlevelcluster);
      for(int k=1; k<record.size(); k++)
      {
        for(int l=0; l<record.at(k)->getCluster().size(); l++)
        {
          record.at(0)->addUnit(record.at(k)->getCluster().at(l));
        }
      }
      for(int m=1; m<record.size(); m++)
      {
        m_output.erase(find(m_output.begin(),m_output.end(),record.at(m)));
        delete record.at(m); 
        record.at(m) = nullptr;
      }
      record.clear();
      lowlevelcluster = nullptr;
      continue;
    }
  
    T2* highlevelcluster = new T2(); //first new
    highlevelcluster->addUnit(lowlevelcluster);
    lowlevelcluster = nullptr;
    m_output.push_back(highlevelcluster);
  }
	return StatusCode::SUCCESS;
	//cout<<"how many neighbors: "<<number<<endl;
}


// StatusCode GlobalClusteringAlg::Towering(std::vector<PandoraPlus::Calo3DCluster*>& m_3dcluster,std::vector<PandoraPlus::CaloTower*>& m_tower)
// {
// 	int record = 0;
// 	for(int i=0; i<m_3dcluster.size(); i++)
// 	{
// 		PandoraPlus::Calo3DCluster* cluster3d = m_3dcluster.at(i);
// 		std::vector<const PandoraPlus::CaloUnit*> m_bar = cluster3d->getBars();
// 		std::vector<PandoraPlus::Calo1DCluster*> m_1dcluster; m_1dcluster.clear();
// 		std::vector<PandoraPlus::Calo2DCluster*> m_2dcluster; m_2dcluster.clear();
// 		//
// 		for(int j=0; j<m_bar.size(); j++)
// 		{
// 			const PandoraPlus::CaloUnit* bar = m_bar.at(j);

// 			for(int k=0; k<m_1dcluster.size(); k++)
// 			{
// 				if(m_1dcluster.at(k)->getBars().at(0)->getModule() == bar->getModule() && m_1dcluster.at(k)->getBars().at(0)->getPart() == bar->getPart() && m_1dcluster.at(k)->getBars().at(0)->getStave() == bar->getStave()
// 				&& m_1dcluster.at(k)->getBars().at(0)->getDlayer() == bar->getDlayer() && m_1dcluster.at(k)->getBars().at(0)->getSlayer() == bar->getSlayer())
// 				{
// 					m_1dcluster.at(k)->addUnit(bar);
// 					record = 1;
// 					break;
// 				}
// 			}

// 			if(record==1)
// 			{
// 				record = 0;
// 				continue;
// 			}

// 			PandoraPlus::Calo1DCluster* cluster1d = new PandoraPlus::Calo1DCluster(); //first new
// 			cluster1d->addUnit(bar);
// 			bar = nullptr;
// 			m_1dcluster.push_back(cluster1d);	
// 		}
// 		//
// 		record = 0;
// 		for(int j=0; j<m_1dcluster.size(); j++)
// 		{
// 			PandoraPlus::Calo1DCluster* cluster1d = m_1dcluster.at(j);
// 			for(int k=0; k<m_2dcluster.size(); k++)
// 			{
// 				if(cluster1d->getTowerID().at(0)==m_2dcluster.at(k)->getTowerID().at(0))
// 				{
// 					m_2dcluster.at(k)->addUnit(cluster1d);
// 					record = 1;
// 					break;
// 				}
// 			}
// 			if(record==1)
// 			{
// 				record = 0;
// 				continue;
// 			}
// 			PandoraPlus::Calo2DCluster* cluster2d = new PandoraPlus::Calo2DCluster(); //first new
// 			cluster2d->addUnit(cluster1d);
// 			cluster1d = nullptr;
// 			m_2dcluster.push_back(cluster2d);
// 		}
// 		//
// 		record = 0;
// 		for(int j=0; j<m_2dcluster.size(); j++)
// 		{
// 			PandoraPlus::Calo2DCluster* cluster2d = m_2dcluster.at(j);
// 			for(int k=0; k<m_tower.size(); k++)
// 			{
// 				if(cluster2d->getTowerID().at(0)==(m_tower.at(k)->getModule()*16*16 + m_tower.at(k)->getPart()*16 + m_tower.at(k)->getStave()))  //(m_tower.at(k)->getModule()*16*16 + m_tower.at(k)->getPart()*16 + m_tower.at(k)->getStave())
// 				{
// 					m_tower.at(k)->addUnit(cluster2d);
// 					record = 1;
// 					break;
// 				}
// 			}
// 			if(record==1)
// 			{
// 				record = 0;
// 				continue;
// 			}
// 			PandoraPlus::CaloTower* tower = new PandoraPlus::CaloTower(); //first new
// 			tower->addUnit(cluster2d);
// 			cluster2d = nullptr;
// 			m_tower.push_back(tower);
// 		}

// 	}
// 	return StatusCode::SUCCESS;
// }

#endif

