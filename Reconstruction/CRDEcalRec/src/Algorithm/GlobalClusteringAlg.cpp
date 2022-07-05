#ifndef GLOBALCLUSTERING_ALG_C
#define GLOBALCLUSTERING_ALG_C

#include "Algorithm/GlobalClusteringAlg.h"
using namespace PandoraPlus;

StatusCode GlobalClusteringAlg::ReadSettings(PandoraPlus::Settings& m_settings){
  settings = m_settings;
/*   if(settings.map_intPars.find("m_module")==settings.map_intPars.end()) settings.map_intPars["m_module"] = 7;
  if(settings.map_intPars.find("m_modulestart")==settings.map_intPars.end()) settings.map_intPars["m_modulestart"] = 0;
  if(settings.map_intPars.find("m_part")==settings.map_intPars.end()) settings.map_intPars["m_part"] = 4;
  if(settings.map_intPars.find("m_stave")==settings.map_intPars.end()) settings.map_intPars["m_stave"] = 11;
  if(settings.map_intPars.find("m_superlayer")==settings.map_intPars.end()) settings.map_intPars["m_superlayer"] = 14;
  if(settings.map_intPars.find("m_startnumber")==settings.map_intPars.end()) settings.map_intPars["m_startnumber"] = 1;
  if(settings.map_intPars.find("m_phibarnumber")==settings.map_intPars.end()) settings.map_intPars["m_phibarnumber"] = 60;
  if(settings.map_intPars.find("m_zbarnumber")==settings.map_intPars.end()) settings.map_intPars["m_zbarnumber"] = 47; */
  return StatusCode::SUCCESS;
};

StatusCode GlobalClusteringAlg::Initialize(){

  return StatusCode::SUCCESS;
};

StatusCode GlobalClusteringAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){
  // //Readin: m_DataCol.BarCol
  // //Output: m_Datacol.BlockCol

  std::vector<PandoraPlus::CaloBar*> m_bars = m_datacol.BarCol; 
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

//   std::vector<PandoraPlus::CaloBar*> m_bars = m_datacol.BarCol; 
	std::vector<PandoraPlus::Calo1DCluster*> m_1dclusters; m_1dclusters.clear();
	std::vector<PandoraPlus::Calo2DCluster*> m_2dclusters; m_2dclusters.clear();
	std::vector<PandoraPlus::Calo3DCluster*> m_3dclusters; m_3dclusters.clear();
	
	cout<<"check bar data: "<<m_bars.size()<<endl;
	Clustering(m_bars, m_1dclusters);
	for(int i1d=0; i1d<m_1dclusters.size(); i1d++)
	{
		m_datacol.bk_Cluster1DCol.push_back(m_1dclusters.at(i1d));
	}
	cout<<"check 1d data: "<<m_1dclusters.size()<<endl;
	Clustering(m_1dclusters, m_2dclusters);
	for(int i2d=0; i2d<m_2dclusters.size(); i2d++)
	{
		m_datacol.bk_Cluster2DCol.push_back(m_2dclusters.at(i2d));
	}
	cout<<"check 2d data: "<<m_2dclusters.size()<<endl;
	Clustering(m_2dclusters, m_3dclusters);
	for(int i3d=0; i3d<m_3dclusters.size(); i3d++)
	{
		m_datacol.bk_Cluster3DCol.push_back(m_3dclusters.at(i3d));
	}
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
    std::vector<int> record;
    record.clear();
	// int number = 0;
	// cout<<"1: "<<record.size();
    for(int i=0; i<m_input.size(); i++)
    {
        T1* lowlevelcluster = m_input.at(i);

        for(int j=0; j<m_output.size(); j++)
        {
            if(ifAdjacent(lowlevelcluster,m_output.at(j))) // //m_output.at(j).isNeighbor(lowlevelcluster)
            {
                record.push_back(j);
				// number = number + 1;
            }
        }
		// cout<<"2: "<<record.size();

        if(record.size()>0)
        {
            sort(record.begin(), record.end(), greater<int>());
            m_output.at(record.at(0))->addCluster(lowlevelcluster);
            for(int k=1; k<record.size(); k++)
            {
                for(int l=0; l<m_output.at(record.at(k))->getCluster().size(); l++)
                {
                    m_output.at(record.at(0))->addCluster(m_output.at(record.at(k))->getCluster().at(l));
                }
            }
            for(int m=1; m<record.size(); m++)
            {
				delete m_output[record.at(m)]; 
				m_output[record.at(m)]=NULL;
                m_output.erase(m_output.begin()+record.at(m)); //memory
            }
            record.clear();
			lowlevelcluster = nullptr;
			// cout<<"4: "<<m_output.size()<<endl;
            continue;
        }
		
		T2* highlevelcluster = new T2(); //first new
        highlevelcluster->addCluster(lowlevelcluster);
		lowlevelcluster = nullptr;
        m_output.push_back(highlevelcluster);
		// cout<<"3: "<<m_output.size()<<endl;
    }
	return StatusCode::SUCCESS;
	//cout<<"how many neighbors: "<<number<<endl;
}

//bar->1d
bool GlobalClusteringAlg::ifAdjacent(PandoraPlus::CaloBar* m_bar, PandoraPlus::Calo1DCluster* cluster1d) //there are six situations
{
	std::vector<const PandoraPlus::CaloBar*> m_1dcluster = cluster1d->getCluster();
	
	for(int i1d = 0; i1d<m_1dcluster.size(); i1d++)
	{
		if(m_bar->getModule()==m_1dcluster[i1d]->getModule() && m_bar->getDlayer()==m_1dcluster[i1d]->getDlayer() && m_bar->getSlayer()==m_1dcluster[i1d]->getSlayer())
		{
			if(m_1dcluster[i1d]->getSlayer()==0)
			{
				if(m_1dcluster[i1d]->getBar()==1)
				{
					if(
					(m_bar->getBar()==m_phibarnumber && m_bar->getStave()==(m_1dcluster[i1d]->getStave()-1) && abs(m_bar->getPart()-m_1dcluster[i1d]->getPart())<=1)||
					((m_bar->getBar()==1 || m_bar->getBar()==2)  && m_bar->getStave()==m_1dcluster[i1d]->getStave() && abs(m_bar->getPart()-m_1dcluster[i1d]->getPart())<=1)
					)
					{
						return true;
					}
				}
				else if(m_1dcluster[i1d]->getBar()==m_phibarnumber)
				{
					if(
					(m_bar->getBar()==1 && m_bar->getStave()==(m_1dcluster[i1d]->getStave()+1) && abs(m_bar->getPart()-m_1dcluster[i1d]->getPart())<=1)||
					((m_bar->getBar()==m_phibarnumber || m_bar->getBar()==(m_phibarnumber-1))  && m_bar->getStave()==m_1dcluster[i1d]->getStave() && abs(m_bar->getPart()-m_1dcluster[i1d]->getPart())<=1)
					)
					{
						return true;
					}
				}
				else
				{
					if(abs(m_bar->getBar()-m_1dcluster[i1d]->getBar())<=1 && m_bar->getStave()==m_1dcluster[i1d]->getStave() && abs(m_bar->getPart()-m_1dcluster[i1d]->getPart())<=1)
					{
						return true;
					}
				}
			}
			else
			{
				if(m_1dcluster[i1d]->getBar()==1)
				{
					if(
					(m_bar->getBar()==(m_zbarnumber-2*(m_bar->getDlayer()-1)) && m_bar->getPart()==(m_1dcluster[i1d]->getPart()-1) && abs(m_bar->getStave()-m_1dcluster[i1d]->getStave())<=1)||
					((m_bar->getBar()==1 || m_bar->getBar()==2)  && m_bar->getPart()==m_1dcluster[i1d]->getPart() && abs(m_bar->getStave()-m_1dcluster[i1d]->getStave())<=1)
					)
					{
						return true;
					}
				}
				else if(m_1dcluster[i1d]->getBar()==(m_zbarnumber-2*(m_bar->getDlayer()-1)))
				{
					if(
					(m_bar->getBar()==1 && m_bar->getPart()==(m_1dcluster[i1d]->getPart()+1) && abs(m_bar->getStave()-m_1dcluster[i1d]->getStave())<=1)||
					((m_bar->getBar()==(m_zbarnumber-2*(m_bar->getDlayer()-1)) || m_bar->getBar()==((m_zbarnumber-2*(m_bar->getDlayer()-1))-1))  && m_bar->getPart()==m_1dcluster[i1d]->getPart() && abs(m_bar->getStave()-m_1dcluster[i1d]->getStave())<=1)
					)
					{
						return true;
					}
				}
				else
				{
					if(abs(m_bar->getBar()-m_1dcluster[i1d]->getBar())<=1 && m_bar->getPart()==m_1dcluster[i1d]->getPart() && abs(m_bar->getStave()-m_1dcluster[i1d]->getStave())<=1)
					{
						return true;
					}
				}
			}
		}		
	}

	
	return false;
}

//1d->2d
bool GlobalClusteringAlg::ifAdjacent(PandoraPlus::Calo1DCluster* m_1dcluster, PandoraPlus::Calo2DCluster* m_2dcluster)
{
	assert(m_1dcluster->getBars().size() > 0 && m_2dcluster->getCluster().at(0)->getCluster().size()>0);
	return (m_1dcluster->getBars().at(0)->getDlayer() == m_2dcluster->getCluster().at(0)->getCluster().at(0)->getDlayer() && ifSameTower(m_1dcluster,m_2dcluster));
}

//2d->3d
bool GlobalClusteringAlg::ifAdjacent(PandoraPlus::Calo2DCluster* m_2dcluster, PandoraPlus::Calo3DCluster* m_3dcluster)
{
	if(ifSameTower(m_2dcluster,m_3dcluster))
	{
		return true;
	}
	
	std::vector<const PandoraPlus::CaloBar*> bars_2d = m_2dcluster->getBars();
	std::vector<const PandoraPlus::CaloBar*> bars_3d = m_3dcluster->getBars();
	
	for(int i=0; i<bars_2d.size(); i++)
	{
		const PandoraPlus::CaloBar* bob = bars_2d.at(i);
		
		for(int j=0; j<bars_3d.size(); j++)
		{
			const PandoraPlus::CaloBar* alice = bars_3d.at(j);
			// cout<<i<<j<<endl;
			if(ifModuleAdjacent(bob,alice))
			{
				// cout<<"ifModuleAdjacent"<<endl;
				return true;
			}
		}
	}

	return false;

	
}

//some other templates and functions
template<typename T1, typename T2> bool GlobalClusteringAlg::ifSameTower(T1* m_uncluster, T2* m_incluster) //1d-->2d and 2d-->3d(normal situation)
{
	std::vector<int> untower = m_uncluster->getTowerID();
	std::vector<int> intower = m_incluster->getTowerID();
    for(int i=0; i<untower.size(); i++)
    {
        for(int j=0; j<intower.size(); j++)
        {
            if(untower.at(i) == intower.at(j))
            {
                return true;
            }
        }
    }
    return false;
}

bool GlobalClusteringAlg::ifModuleAdjacent(const PandoraPlus::CaloBar* bar_2d, const PandoraPlus::CaloBar* bar_3d)
{
	PandoraPlus::CaloBar bob; //just first layer
	PandoraPlus::CaloBar alice; //second new
	// m_datacol.bk_BarCol.push_back(bob);
	// m_datacol.bk_BarCol.push_back(alice);
	if(bar_2d->getModule()==m_module && bar_3d->getModule()==m_modulestart)
	{
		bob = *bar_2d;
		alice = *bar_3d;
		return ifAdjacent(bob, alice);
		
	}
	else if(bar_2d->getModule()==m_modulestart && bar_3d->getModule()==m_module)
	{
		bob = *bar_3d;
		alice = *bar_2d;
		return ifAdjacent(bob, alice);
	}
	else if(bar_2d->getModule() == bar_3d->getModule())
	{
		return false;
	}
	else if((bar_2d->getModule() < bar_3d->getModule()) && (bar_2d->getModule() - bar_3d->getModule())==-1)
	{
		bob = *bar_2d;
		alice = *bar_3d;
		return ifAdjacent(bob, alice);
	}
	else if((bar_2d->getModule() > bar_3d->getModule()) && (bar_2d->getModule() - bar_3d->getModule())==1)
	{
		bob = *bar_3d;
		alice = *bar_2d;
		return ifAdjacent(bob, alice);
	}	
	else
	{
		return false;
	}	
}

//
bool GlobalClusteringAlg::ifAdjacent(PandoraPlus::CaloBar &bob, PandoraPlus::CaloBar &alice)
{
	if( (bob.getDlayer()==1 && bob.getSlayer()==0 && bob.getPart()==m_part && alice.getSlayer()==0 && alice.getPart()==1 && bob.getStave()==alice.getStave() && abs(bob.getBar()-alice.getBar())<=1 ) ||
		(bob.getDlayer()==1 && bob.getSlayer()==0 && bob.getPart()==m_part && alice.getSlayer()==1 && alice.getPart()==1 && bob.getStave()==alice.getStave() && alice.getBar()==1) ||
		(bob.getDlayer()==1 && bob.getSlayer()==0 && bob.getPart()==m_part && alice.getSlayer()==0 && alice.getPart()==1 && abs(bob.getStave()-alice.getStave())<=1 && ((bob.getBar()==1&&alice.getBar()==60)||(bob.getBar()==60&&alice.getBar()==1)))
	) 
	{
		return true;
	}
	else
	{
		return false;
	}
}
//
StatusCode GlobalClusteringAlg::Towering(std::vector<PandoraPlus::Calo3DCluster*>& m_3dcluster,std::vector<PandoraPlus::CaloTower*>& m_tower)
{
	int record = 0;
	for(int i=0; i<m_3dcluster.size(); i++)
	{
		PandoraPlus::Calo3DCluster* cluster3d = m_3dcluster.at(i);
		std::vector<const PandoraPlus::CaloBar*> m_bar = cluster3d->getBars();
		std::vector<PandoraPlus::Calo1DCluster*> m_1dcluster; m_1dcluster.clear();
		std::vector<PandoraPlus::Calo2DCluster*> m_2dcluster; m_2dcluster.clear();
		//
		for(int j=0; j<m_bar.size(); j++)
		{
			const PandoraPlus::CaloBar* bar = m_bar.at(j);

			for(int k=0; k<m_1dcluster.size(); k++)
			{
				if(m_1dcluster.at(k)->getBars().at(0)->getModule() == bar->getModule() && m_1dcluster.at(k)->getBars().at(0)->getPart() == bar->getPart() && m_1dcluster.at(k)->getBars().at(0)->getStave() == bar->getStave()
				&& m_1dcluster.at(k)->getBars().at(0)->getDlayer() == bar->getDlayer() && m_1dcluster.at(k)->getBars().at(0)->getSlayer() == bar->getSlayer())
				{
					m_1dcluster.at(k)->addCluster(bar);
					record = 1;
					break;
				}
			}

			if(record==1)
			{
				record = 0;
				continue;
			}

			PandoraPlus::Calo1DCluster* cluster1d = new PandoraPlus::Calo1DCluster(); //first new
			cluster1d->addCluster(bar);
			bar = nullptr;
			m_1dcluster.push_back(cluster1d);	
		}
		//
		record = 0;
		for(int j=0; j<m_1dcluster.size(); j++)
		{
			PandoraPlus::Calo1DCluster* cluster1d = m_1dcluster.at(j);
			for(int k=0; k<m_2dcluster.size(); k++)
			{
				if(cluster1d->getTowerID().at(0)==m_2dcluster.at(k)->getTowerID().at(0))
				{
					m_2dcluster.at(k)->addCluster(cluster1d);
					record = 1;
					break;
				}
			}
			if(record==1)
			{
				record = 0;
				continue;
			}
			PandoraPlus::Calo2DCluster* cluster2d = new PandoraPlus::Calo2DCluster(); //first new
			cluster2d->addCluster(cluster1d);
			cluster1d = nullptr;
			m_2dcluster.push_back(cluster2d);
		}
		//
		record = 0;
		for(int j=0; j<m_2dcluster.size(); j++)
		{
			PandoraPlus::Calo2DCluster* cluster2d = m_2dcluster.at(j);
			for(int k=0; k<m_tower.size(); k++)
			{
				if(cluster2d->getTowerID().at(0)==(m_tower.at(k)->getModule()*16*16 + m_tower.at(k)->getPart()*16 + m_tower.at(k)->getStave()))  //(m_tower.at(k)->getModule()*16*16 + m_tower.at(k)->getPart()*16 + m_tower.at(k)->getStave())
				{
					m_tower.at(k)->addCluster(cluster2d);
					record = 1;
					break;
				}
			}
			if(record==1)
			{
				record = 0;
				continue;
			}
			PandoraPlus::CaloTower* tower = new PandoraPlus::CaloTower(); //first new
			tower->addCluster(cluster2d);
			cluster2d = nullptr;
			m_tower.push_back(tower);
		}

	}
	return StatusCode::SUCCESS;
}

#endif

