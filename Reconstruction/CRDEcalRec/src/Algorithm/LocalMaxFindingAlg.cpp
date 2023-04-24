#ifndef _LOCALMAXFINDING_ALG_C
#define _LOCALMAXFINDING_ALG_C

#include "Algorithm/LocalMaxFindingAlg.h"
using namespace PandoraPlus; 

StatusCode LocalMaxFindingAlg::ReadSettings(PandoraPlus::Settings& m_settings){
  settings = m_settings;

  if(settings.map_floatPars.find("Eth_localMax")==settings.map_floatPars.end()) settings.map_floatPars["Eth_localMax"] = 0.005;
  if(settings.map_floatPars.find("Eth_MaxWithNeigh")==settings.map_floatPars.end()) settings.map_floatPars["Eth_MaxWithNeigh"] = 0.;
  if(settings.map_stringPars.find("OutputLocalMaxName")==settings.map_stringPars.end()) settings.map_stringPars["OutputLocalMaxName"] = "AllLocalMax";
  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  p_HalfClusU = NULL;
  p_HalfClusV = NULL; 

  p_HalfClusU = &(m_datacol.map_HalfCluster["HalfClusterColU"]);
  p_HalfClusV = &(m_datacol.map_HalfCluster["HalfClusterColV"]);


  return StatusCode::SUCCESS;
}


StatusCode LocalMaxFindingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol){

  if(!p_HalfClusU || !p_HalfClusV) {std::cout<<"ERROR: No HalfClusters in present data collection! "<<std::endl; return StatusCode::FAILURE; }
printf("  LocalMaxFindingAlg: Readin HalfCluster size (%d, %d) \n", p_HalfClusU->size(), p_HalfClusV->size());

  for(int iu = 0; iu<p_HalfClusU->size(); iu++){
    std::vector<const Calo1DCluster*> m_1dClusCol = p_HalfClusU->at(iu).get()->getCluster();

    std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>> ptr_localMax; ptr_localMax.clear();
    for(int i1d=0; i1d<m_1dClusCol.size(); i1d++) GetLocalMax(m_1dClusCol[i1d], ptr_localMax);

    std::vector<const PandoraPlus::Calo1DCluster*> tmp_localMax; tmp_localMax.clear();
    for(auto iter : ptr_localMax){
      tmp_localMax.push_back(iter.get());
      m_datacol.map_1DCluster["bk1DCluster"].push_back(iter);
    }
    p_HalfClusU->at(iu).get()->setLocalMax(settings.map_stringPars["OutputLocalMaxName"], tmp_localMax);
    m_1dClusCol.clear(); 
    ptr_localMax.clear();
  }

  for(int iv=0; iv<p_HalfClusV->size(); iv++){
    std::vector<const Calo1DCluster*> m_1dClusCol = p_HalfClusV->at(iv).get()->getCluster();

    std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>> ptr_localMax; ptr_localMax.clear();
    for(int i1d=0; i1d<m_1dClusCol.size(); i1d++) GetLocalMax(m_1dClusCol[i1d], ptr_localMax);

    std::vector<const PandoraPlus::Calo1DCluster*> tmp_localMax; tmp_localMax.clear();
    for(auto iter : ptr_localMax){
      tmp_localMax.push_back(iter.get());
      m_datacol.map_1DCluster["bk1DCluster"].push_back(iter);
    }
    p_HalfClusV->at(iv).get()->setLocalMax(settings.map_stringPars["OutputLocalMaxName"], tmp_localMax);
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

StatusCode LocalMaxFindingAlg::GetLocalMax( const PandoraPlus::Calo1DCluster* m_1dClus, 
                                            std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>>& m_output){

  if(m_1dClus->getBars().size()==0) return StatusCode::SUCCESS;

  std::vector<const PandoraPlus::CaloUnit*> m_barCol = m_1dClus->getBars();

//cout<<"  LocalMaxFindingAlg::GetLocalMax: Input bar collection size: "<<m_barUCol.size()<<"  "<<m_barVCol.size()<<endl;

  std::vector<const PandoraPlus::CaloUnit*> localMaxCol; localMaxCol.clear();

  GetLocalMaxBar( m_barCol, localMaxCol );

//cout<<"  LocalMaxFindingAlg::GetLocalMax: Found local max bar size: "<<localMaxUCol.size()<<"  "<<localMaxVCol.size()<<endl;
//cout<<"  Transfer bar to barShower"<<endl;

  for(int j=0; j<localMaxCol.size(); j++){
    //PandoraPlus::Calo1DCluster* m_shower = new PandoraPlus::Calo1DCluster();
    std::shared_ptr<PandoraPlus::Calo1DCluster> m_shower = std::make_shared<PandoraPlus::Calo1DCluster>();
    m_shower->addUnit( localMaxCol[j] );
    m_output.push_back(m_shower);
  }

//cout<<"  Output bar shower size: "<<m_outputU.size()<<"  "<<m_outputV.size()<<endl;
//cout<<endl;

  return StatusCode::SUCCESS;
}

StatusCode LocalMaxFindingAlg::GetLocalMaxBar( std::vector<const PandoraPlus::CaloUnit*>& barCol, std::vector<const PandoraPlus::CaloUnit*>& localMaxCol ){
  //std::sort( barCol.begin(), barCol.end(), compBar );

  for(int ib=0; ib<barCol.size(); ib++){
    std::vector<const PandoraPlus::CaloUnit*> m_neighbors = getNeighbors( barCol[ib], barCol );
    if( m_neighbors.size()==0 && barCol[ib]->getEnergy()>settings.map_floatPars["Eth_localMax"] ) { 
      localMaxCol.push_back( barCol[ib] ); continue; 
    }

    bool isLocalMax=true;
    double Eneigh=0;
    if(barCol[ib]->getEnergy()<settings.map_floatPars["Eth_localMax"]) isLocalMax = false;
    for(int j=0;j<m_neighbors.size();j++){
      if(m_neighbors[j]->getEnergy()>barCol[ib]->getEnergy()) isLocalMax=false;
      Eneigh += m_neighbors[j]->getEnergy();
    }
    if( (barCol[ib]->getEnergy()/(barCol[ib]->getEnergy()+Eneigh))<settings.map_floatPars["Eth_MaxWithNeigh"] ) isLocalMax = false;

    if(isLocalMax) localMaxCol.push_back( barCol[ib] );
  }

  return StatusCode::SUCCESS;
}


std::vector<const PandoraPlus::CaloUnit*> LocalMaxFindingAlg::getNeighbors( const PandoraPlus::CaloUnit* seed, std::vector<const PandoraPlus::CaloUnit*>& barCol){
  std::vector<const PandoraPlus::CaloUnit*> m_neighbor; m_neighbor.clear();
  for(int i=0;i<barCol.size();i++){
    bool fl_neighbor = false; 
    if( seed->getModule()==barCol[i]->getModule() && 
        seed->getPart()==barCol[i]->getPart() && 
        seed->getStave()==barCol[i]->getStave() && 
        seed->getDlayer()==barCol[i]->getDlayer() &&
        seed->getSlayer()==barCol[i]->getSlayer() &&
        abs( seed->getBar()-barCol[i]->getBar() )==1 ) fl_neighbor=true;
    else if( seed->getModule()==barCol[i]->getModule() && 
             seed->getStave()==barCol[i]->getStave() && 
             ( ( seed->getPart()-barCol[i]->getPart()==1 && seed->isAtLowerEdgePhi() && barCol[i]->isAtUpperEdgePhi() ) ||
               ( barCol[i]->getPart()-seed->getPart()==1 && seed->isAtUpperEdgePhi() && barCol[i]->isAtLowerEdgePhi() ) ) ) fl_neighbor=true;
    else if( seed->getModule()==barCol[i]->getModule() && 
             seed->getPart()==barCol[i]->getPart() &&
             ( ( seed->getStave()-barCol[i]->getStave()==1 && seed->isAtLowerEdgeZ() && barCol[i]->isAtUpperEdgeZ() ) || 
               ( barCol[i]->getStave()-seed->getStave()==1 && seed->isAtUpperEdgeZ() && barCol[i]->isAtLowerEdgeZ() ) ) ) fl_neighbor=true;


    if(fl_neighbor) m_neighbor.push_back( barCol[i] );

  }
  if(m_neighbor.size()>2) std::cout<<"WARNING: more than 2 hits in neighborCol!!"<<std::endl;

  return m_neighbor;
}

#endif
