#ifndef GLOBALCLUSTERING_ALG_C
#define GLOBALCLUSTERING_ALG_C

#include "Algorithm/GlobalClusteringAlg.h"
using namespace PandoraPlus;

StatusCode GlobalClusteringAlg::ReadSettings(PandoraPlus::Settings& m_settings){
  settings = m_settings;

  return StatusCode::SUCCESS;
};

StatusCode GlobalClusteringAlg::Initialize(){

  return StatusCode::SUCCESS;
};

StatusCode GlobalClusteringAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

/*
  std::vector<PandoraPlus::CaloUnit*> m_bars = m_datacol.BarCol; 
  std::vector<PandoraPlus::CaloBlock*> m_blocks; m_blocks.clear();

  for(int ibar=0; ibar<m_bars.size(); ibar++){

    bool fl_foundbl = false;
    for(int ibl=0; ibl<m_blocks.size(); ibl++)
      if( m_bars[ibar]->getModule() == m_blocks[ibl]->getModule() &&
          m_bars[ibar]->getStave()  == m_blocks[ibl]->getStave() &&
          m_bars[ibar]->getPart()   == m_blocks[ibl]->getPart() &&
          m_bars[ibar]->getDlayer() == m_blocks[ibl]->getDlayer() ){
        m_blocks[ibl]->addBar( m_bars[ibar] );
        fl_foundbl=true;
      }

    if(!fl_foundbl){
      PandoraPlus::CaloBlock* m_block = new PandoraPlus::CaloBlock();
      m_datacol.bk_BlockCol.push_back(m_block);

      m_block->addBar( m_bars[ibar] );
      m_block->setIDInfo(m_bars[ibar]->getModule(),  m_bars[ibar]->getStave(), m_bars[ibar]->getDlayer(), m_bars[ibar]->getPart());
      m_blocks.push_back(m_block);
    }
  }
  m_datacol.BlockCol = m_blocks;
*/

   std::vector<PandoraPlus::CaloUnit*> m_bars = m_datacol.BarCol; 
	std::vector<PandoraPlus::Calo1DCluster*> m_1dclusters; m_1dclusters.clear();
	std::vector<PandoraPlus::Calo2DCluster*> m_2dclusters; m_2dclusters.clear();
	std::vector<PandoraPlus::Calo3DCluster*> m_3dclusters; m_3dclusters.clear();
	
cout<<"check bar data: "<<m_bars.size()<<endl;
	Clustering(m_bars, m_1dclusters);
	for(int i1d=0; i1d<m_1dclusters.size(); i1d++)
		m_datacol.bk_Cluster1DCol.push_back(m_1dclusters.at(i1d));


cout<<"check 1d data: "<<m_1dclusters.size()<<endl;
	Clustering(m_1dclusters, m_2dclusters);
	for(int i2d=0; i2d<m_2dclusters.size(); i2d++)
		m_datacol.bk_Cluster2DCol.push_back(m_2dclusters.at(i2d));

cout<<"check 2d data: "<<m_2dclusters.size()<<endl;
	Clustering(m_2dclusters, m_3dclusters);
	for(int i3d=0; i3d<m_3dclusters.size(); i3d++)
		m_datacol.bk_Cluster3DCol.push_back(m_3dclusters.at(i3d));

cout<<"check 3d data: "<<m_3dclusters.size()<<endl;

	m_datacol.Cluster1DCol = m_1dclusters;
	m_datacol.Cluster2DCol = m_2dclusters;
	m_datacol.Cluster3DCol = m_3dclusters;

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
      record.at(0)->addCluster(lowlevelcluster);
      for(int k=1; k<record.size(); k++)
      {
        for(int l=0; l<record.at(k)->getCluster().size(); l++)
        {
          record.at(0)->addCluster(record.at(k)->getCluster().at(l));
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
    highlevelcluster->addCluster(lowlevelcluster);
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
// 					m_1dcluster.at(k)->addCluster(bar);
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
// 			cluster1d->addCluster(bar);
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
// 					m_2dcluster.at(k)->addCluster(cluster1d);
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
// 			cluster2d->addCluster(cluster1d);
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
// 					m_tower.at(k)->addCluster(cluster2d);
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
// 			tower->addCluster(cluster2d);
// 			cluster2d = nullptr;
// 			m_tower.push_back(tower);
// 		}

// 	}
// 	return StatusCode::SUCCESS;
// }

#endif

