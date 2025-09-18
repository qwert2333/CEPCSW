#ifndef _LOCALMAXFINDING_ALG_C
#define _LOCALMAXFINDING_ALG_C

#include "Algorithm/LocalMaxFindingAlg.h"
using namespace Cyber; 

StatusCode LocalMaxFindingAlg::ReadSettings(Cyber::Settings& m_settings){
  settings = m_settings;

  if(settings.map_floatPars.find("Eth_localMax")==settings.map_floatPars.end()) settings.map_floatPars["Eth_localMax"] = 0.0075;
  if(settings.map_floatPars.find("Eth_MaxWithNeigh")==settings.map_floatPars.end()) settings.map_floatPars["Eth_MaxWithNeigh"] = 0.;
  if(settings.map_stringPars.find("OutputLocalMaxName")==settings.map_stringPars.end()) settings.map_stringPars["OutputLocalMaxName"] = "AllLocalMax";
  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::Initialize( CyberDataCol& m_datacol ){
  p_HalfClusU = NULL;
  p_HalfClusV = NULL; 

  p_HalfClusU = &(m_datacol.map_HalfCluster["HalfClusterColU"]);
  p_HalfClusV = &(m_datacol.map_HalfCluster["HalfClusterColV"]);


  return StatusCode::SUCCESS;
}


StatusCode LocalMaxFindingAlg::RunAlgorithm( CyberDataCol& m_datacol){

  if(!p_HalfClusU || !p_HalfClusV) {std::cout<<"ERROR: No HalfClusters in present data collection! "<<std::endl; return StatusCode::FAILURE; }

  for(int iu = 0; iu<p_HalfClusU->size(); iu++){
    std::vector<const Calo1DCluster*> m_1dClusCol = p_HalfClusU->at(iu).get()->getCluster();

    std::vector<std::shared_ptr<Cyber::Calo1DCluster>> ptr_localMax; ptr_localMax.clear();
    for(int i1d=0; i1d<m_1dClusCol.size(); i1d++) GetLocalMax(m_1dClusCol[i1d], ptr_localMax);

    std::vector<const Cyber::Calo1DCluster*> tmp_localMax; tmp_localMax.clear();
    for(auto iter : ptr_localMax){
      iter->getLinkedMCPfromUnit();
      tmp_localMax.push_back(iter.get());
      m_datacol.map_1DCluster["bk1DCluster"].push_back(iter);
    }
    p_HalfClusU->at(iu).get()->setLocalMax(settings.map_stringPars["OutputLocalMaxName"], tmp_localMax);
printf("  HalfClusterU #%d: energy %.3f, localMax size %d \n", iu, p_HalfClusU->at(iu).get()->getEnergy(), tmp_localMax.size());
    m_1dClusCol.clear(); 
    ptr_localMax.clear();
  }

  for(int iv=0; iv<p_HalfClusV->size(); iv++){
    std::vector<const Calo1DCluster*> m_1dClusCol = p_HalfClusV->at(iv).get()->getCluster();

    std::vector<std::shared_ptr<Cyber::Calo1DCluster>> ptr_localMax; ptr_localMax.clear();
    for(int i1d=0; i1d<m_1dClusCol.size(); i1d++) GetLocalMax(m_1dClusCol[i1d], ptr_localMax);

    std::vector<const Cyber::Calo1DCluster*> tmp_localMax; tmp_localMax.clear();
    for(auto iter : ptr_localMax){
      iter->getLinkedMCPfromUnit();
      tmp_localMax.push_back(iter.get());
      m_datacol.map_1DCluster["bk1DCluster"].push_back(iter);
    }
    p_HalfClusV->at(iv).get()->setLocalMax(settings.map_stringPars["OutputLocalMaxName"], tmp_localMax);
printf("  HalfClusterV #%d: energy %.3f, localMax size %d \n", iv, p_HalfClusV->at(iv).get()->getEnergy(), tmp_localMax.size());
    m_1dClusCol.clear(); 
    ptr_localMax.clear();
  }

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::ClearAlgorithm(){
  p_HalfClusU = nullptr;
  p_HalfClusV = nullptr;

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::GetLocalMax( const Cyber::Calo1DCluster* m_1dClus, 
                                            std::vector<std::shared_ptr<Cyber::Calo1DCluster>>& m_output){

  if(m_1dClus->getBars().size()==0) return StatusCode::SUCCESS;

  std::vector<const Cyber::CaloUnit*> m_barCol = m_1dClus->getBars();

//cout<<"  LocalMaxFindingAlg::GetLocalMax: Input bar collection size: "<<m_barCol.size()<<endl;

  std::vector<const Cyber::CaloUnit*> localMaxCol; localMaxCol.clear();

  GetLocalMaxBar( m_barCol, localMaxCol );

//cout<<"  LocalMaxFindingAlg::GetLocalMax: Found local max bar size: "<<localMaxCol.size()<<endl;
//cout<<"  Transfer bar to barShower"<<endl;

  for(int j=0; j<localMaxCol.size(); j++){
    std::shared_ptr<Cyber::Calo1DCluster> m_shower = std::make_shared<Cyber::Calo1DCluster>();
    m_shower->addUnit( localMaxCol[j] );
    m_output.push_back(m_shower);
  }

//cout<<"  Output bar shower size: "<<m_output.size()<<endl;
//cout<<endl;

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::GetLocalMaxBar( std::vector<const Cyber::CaloUnit*>& barCol, std::vector<const Cyber::CaloUnit*>& localMaxCol ){
  //std::sort( barCol.begin(), barCol.end(), compBar );
  for(int ib=0; ib<barCol.size(); ib++){
    std::vector<const Cyber::CaloUnit*> m_neighbors = getNeighbors( barCol[ib], barCol );
    if( m_neighbors.size()==0 && barCol[ib]->getEnergy()>settings.map_floatPars["Eth_localMax"] ) { 
      localMaxCol.push_back( barCol[ib] ); continue; 
    }

    bool isLocalMax=true;
    double Eneigh=0;
    if(barCol[ib]->getEnergy()<settings.map_floatPars["Eth_localMax"]) isLocalMax = false;
    for(int j=0;j<m_neighbors.size();j++){
      if(m_neighbors[j]->getEnergy() > barCol[ib]->getEnergy()) isLocalMax=false;
      Eneigh += m_neighbors[j]->getEnergy();
    }
    if( (barCol[ib]->getEnergy()/(barCol[ib]->getEnergy()+Eneigh))<settings.map_floatPars["Eth_MaxWithNeigh"] ) isLocalMax = false;
    if(isLocalMax) localMaxCol.push_back( barCol[ib] );
  }

  return StatusCode::SUCCESS;
}


std::vector<const Cyber::CaloUnit*> LocalMaxFindingAlg::getNeighbors( const Cyber::CaloUnit* seed, std::vector<const Cyber::CaloUnit*>& barCol){
  std::vector<const Cyber::CaloUnit*> m_neighbor; m_neighbor.clear();
  for(int i=0;i<barCol.size();i++){
    //bool fl_neighbor = false; 
    //if( seed->getSystem()==barCol[i]->getSystem() &&
    //    seed->getModule()==barCol[i]->getModule() && 
    //    seed->getStave()==barCol[i]->getStave() && 
    //    seed->getDlayer()==barCol[i]->getDlayer() &&
    //    seed->getSlayer()==barCol[i]->getSlayer() &&
    //    abs( seed->getBar()-barCol[i]->getBar() )==1 ) fl_neighbor=true;
    //else if( seed->getSystem()==barCol[i]->getSystem() && seed->getStave()==barCol[i]->getStave() && 
    //         ( ( seed->getModule()-barCol[i]->getModule()==1 && seed->isAtLowerEdgePhi() && barCol[i]->isAtUpperEdgePhi() ) ||
    //           ( barCol[i]->getModule()-seed->getModule()==1 && seed->isAtUpperEdgePhi() && barCol[i]->isAtLowerEdgePhi() ) ) ) fl_neighbor=true;
    //else if( seed->getSystem()==barCol[i]->getSystem() && seed->getModule()==barCol[i]->getModule() && 
    //         ( ( seed->getStave()-barCol[i]->getStave()==1 && seed->isAtLowerEdgeZ() && barCol[i]->isAtUpperEdgeZ() ) || 
    //           ( barCol[i]->getStave()-seed->getStave()==1 && seed->isAtUpperEdgeZ() && barCol[i]->isAtLowerEdgeZ() ) ) ) fl_neighbor=true;

    bool fl_neighbor = seed->isNeighbor( barCol[i] );
    if(fl_neighbor){ 
      if(seed->getSystem()==CaloUnit::System_Barrel && seed->getSlayer()==0 && seed->getModule()!=barCol[i]->getModule()) continue;
      if(seed->getSystem()==CaloUnit::System_Barrel && seed->getSlayer()==1 && seed->getStave()!=barCol[i]->getStave())   continue;
      if(seed->getSystem()==CaloUnit::System_Endcap && seed->getSlayer()==0 && seed->getPart()!=barCol[i]->getPart() )    continue;
      if(seed->getSystem()==CaloUnit::System_Endcap && seed->getSlayer()==1 && seed->getStave()!=barCol[i]->getStave() )  continue;
      m_neighbor.push_back(barCol[i]);
    }
  }
  if(m_neighbor.size()>2){ 
    std::cout<<"WARNING: "<<m_neighbor.size()<<" hits in neighborCol!!"<<std::endl;
    printf("Seed cellID: ( %d, %d, %d, %d, %d, %d, %d), position (%.3f, %.3f, %.3f) \n", seed->getSystem(), seed->getModule(), seed->getStave(), seed->getPart(), seed->getDlayer(), seed->getSlayer(), seed->getBar(), seed->getPosition().x(), seed->getPosition().y(), seed->getPosition().z());
    for(int i=0; i<m_neighbor.size(); i++)
      printf("Bar #%d cellID: ( %d, %d, %d, %d, %d, %d, %d), position (%.3f, %.3f, %.3f) \n", i, m_neighbor[i]->getSystem(), m_neighbor[i]->getModule(), m_neighbor[i]->getStave(), m_neighbor[i]->getPart(), m_neighbor[i]->getDlayer(), m_neighbor[i]->getSlayer(), m_neighbor[i]->getBar(), m_neighbor[i]->getPosition().x(), m_neighbor[i]->getPosition().y(), m_neighbor[i]->getPosition().z());
  
  }

  return m_neighbor;
}

#endif
