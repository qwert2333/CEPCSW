#ifndef _AXISMERGING_ALG_C
#define _AXISMERGING_ALG_C

#include "Algorithm/AxisMergingAlg.h"

StatusCode AxisMergingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_stringPars.find("OutputAxisName")==settings.map_stringPars.end()) settings.map_stringPars["OutputAxisName"] = "MergedAxis";
  if(settings.map_floatPars.find("axis_Angle")==settings.map_floatPars.end()) settings.map_floatPars["axis_Angle"] = TMath::Pi()/4.;
  if(settings.map_floatPars.find("relP_Angle")==settings.map_floatPars.end()) settings.map_floatPars["relP_Angle"] = TMath::Pi()/4.;
  if(settings.map_floatPars.find("skipLayer")==settings.map_floatPars.end()) settings.map_floatPars["skipLayer"] = 3;

  return StatusCode::SUCCESS;
};

StatusCode AxisMergingAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  m_axisUCol.clear();
  m_axisVCol.clear();
  m_newAxisUCol.clear();
  m_newAxisVCol.clear();

  p_HalfClustersU = &(m_datacol.map_HalfCluster["HalfClusterColU"]);
  p_HalfClustersV = &(m_datacol.map_HalfCluster["HalfClusterColV"]);

  return StatusCode::SUCCESS;
};

StatusCode AxisMergingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  if( !p_HalfClustersU || !p_HalfClustersV ){
    std::cout<<"AxisMergingAlg: Input HalfCluster is NULL! "<<std::endl;
    return StatusCode::SUCCESS;
  }
  if( p_HalfClustersU->size() + p_HalfClustersV->size()==0 ) {
    std::cout<<"AxisMergingAlg: No HalfCluster input"<<std::endl;
    return StatusCode::SUCCESS;
  }

cout<<"AxisMergingAlg: Readin halfcluster size: "<<p_HalfClustersU->size()<<", "<<p_HalfClustersV->size()<<endl;


  //Merge axis in HalfClusterU: 
  for(int ih=0; ih<p_HalfClustersU->size(); ih++){
    m_axisUCol.clear(); m_newAxisUCol.clear();
    m_axisUCol = p_HalfClustersU->at(ih)->getAllHalfClusterCol();

    for(int ic=0; ic<m_axisUCol.size(); ic++) m_newAxisUCol.push_back( m_axisUCol[ic]->Clone() );

    if(m_newAxisUCol.size()<2){ 
      std::vector<const PandoraPlus::CaloHalfCluster*> tmp_axisCol; tmp_axisCol.clear();
      for(int ic=0; ic<m_newAxisUCol.size(); ic++) tmp_axisCol.push_back(m_newAxisUCol[ic]);
      p_HalfClustersU->at(ih)->setHalfClusters(settings.map_stringPars["OutputAxisName"], tmp_axisCol);
      continue; 
    }

printf("  In HalfClusterU #%d: readin axis size %d \n", ih, m_newAxisUCol.size());

    //Case1: Merge track axes and other axes.


    //Case2: Merge each type of axes. 
    for(int iax=0; iax<m_newAxisUCol.size(); iax++){
      const PandoraPlus::CaloHalfCluster* m_axis = m_newAxisUCol[iax];
      for(int jax=iax+1; jax<m_newAxisUCol.size(); jax++){
        const PandoraPlus::CaloHalfCluster* p_axis = m_newAxisUCol[jax];
 
        double minRelAngle = min( sin( (m_axis->getPos()-p_axis->getPos()).Angle(m_axis->getAxis())), 
                                  sin( (m_axis->getPos()-p_axis->getPos()).Angle(p_axis->getAxis())) );
        int skipLayer = m_axis->getBeginningDlayer()<p_axis->getBeginningDlayer() ? 
                        (p_axis->getBeginningDlayer()-m_axis->getEndDlayer()) : (m_axis->getBeginningDlayer()-p_axis->getEndDlayer());

printf("Loop check pair: (%d, %d), criteria: theta1 = %.3f, theta2 = %.3f, gap = %d \n", 
iax, jax, m_axis->getAxis().Angle(p_axis->getAxis()), minRelAngle, skipLayer);

        if( sin( m_axis->getAxis().Angle(p_axis->getAxis()) ) < sin(settings.map_floatPars["axis_Angle"]) &&
            sin(minRelAngle) < sin(settings.map_floatPars["relP_Angle"]) && 
            skipLayer <= settings.map_floatPars["skipLayer"] ){

          m_newAxisUCol[iax]->mergeHalfCluster( m_newAxisUCol[jax] );
          delete m_newAxisUCol[jax]; m_newAxisUCol[jax]=NULL;
          m_newAxisUCol.erase(m_newAxisUCol.begin()+jax);
          jax--;
          if(iax>jax+1) iax--;
        }
        p_axis = nullptr;
      }
      m_axis = nullptr;
    }
cout<<"  After merging: axis size "<<m_newAxisUCol.size()<<", Check the overlap"<<endl;

    //convert to constant object and save in HalfCluster: 
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_axisCol; tmp_axisCol.clear();
    for(int ic=0; ic<m_newAxisUCol.size(); ic++) tmp_axisCol.push_back(m_newAxisUCol[ic]);
    p_HalfClustersU->at(ih)->setHalfClusters(settings.map_stringPars["OutputAxisName"], tmp_axisCol);
  }


  //Merge axis in HalfClusterV:
  for(int ih=0; ih<p_HalfClustersV->size(); ih++){
    m_axisVCol.clear(); m_newAxisVCol.clear();
    m_axisVCol = p_HalfClustersV->at(ih)->getAllHalfClusterCol();

    for(int ic=0; ic<m_axisVCol.size(); ic++) m_newAxisVCol.push_back( m_axisVCol[ic]->Clone() );

    if(m_newAxisVCol.size()<2){
      std::vector<const PandoraPlus::CaloHalfCluster*> tmp_axisCol; tmp_axisCol.clear();
      for(int ic=0; ic<m_newAxisVCol.size(); ic++) tmp_axisCol.push_back(m_newAxisVCol[ic]);
      p_HalfClustersV->at(ih)->setHalfClusters(settings.map_stringPars["OutputAxisName"], tmp_axisCol);
      continue;
    }

printf("  In HalfClusterV #%d: readin axis size %d \n", ih, m_newAxisVCol.size());

    //Case1: Merge track axes and other axes.


    //Case2: Merge each type of axes.
    for(int iax=0; iax<m_newAxisVCol.size(); iax++){
      const PandoraPlus::CaloHalfCluster* m_axis = m_newAxisVCol[iax];
      for(int jax=iax+1; jax<m_newAxisVCol.size(); jax++){
        const PandoraPlus::CaloHalfCluster* p_axis = m_newAxisVCol[jax];

        double minRelAngle = min( sin( (m_axis->getPos()-p_axis->getPos()).Angle(m_axis->getAxis())),
                                  sin( (m_axis->getPos()-p_axis->getPos()).Angle(p_axis->getAxis())) );
        int skipLayer = m_axis->getBeginningDlayer()<p_axis->getBeginningDlayer() ?
                        (p_axis->getBeginningDlayer()-m_axis->getEndDlayer()) : (m_axis->getBeginningDlayer()-p_axis->getEndDlayer());
        if( sin( m_axis->getAxis().Angle(p_axis->getAxis()) ) < sin(settings.map_floatPars["axis_Angle"]) &&
            sin(minRelAngle) < sin(settings.map_floatPars["relP_Angle"]) &&
            skipLayer <= settings.map_floatPars["skipLayer"] ){

          m_newAxisVCol[iax]->mergeHalfCluster( m_newAxisVCol[jax] );
          delete m_newAxisVCol[jax]; m_newAxisVCol[jax]=NULL;
          m_newAxisVCol.erase(m_newAxisVCol.begin()+jax);
          jax--;
          if(iax>jax+1) iax--;
        }
        p_axis = nullptr;
      }
      m_axis = nullptr;
    }
cout<<"  After merging: axis size "<<m_newAxisVCol.size()<<", Check the overlap"<<endl;

    //convert to constant object and save in HalfCluster:
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_axisCol; tmp_axisCol.clear();
    for(int ic=0; ic<m_newAxisVCol.size(); ic++) tmp_axisCol.push_back(m_newAxisVCol[ic]);
    p_HalfClustersV->at(ih)->setHalfClusters(settings.map_stringPars["OutputAxisName"], tmp_axisCol);
  }

  return StatusCode::SUCCESS;
};

StatusCode AxisMergingAlg::ClearAlgorithm(){
  p_HalfClustersU = NULL;
  p_HalfClustersV = NULL;
  m_axisUCol.clear();
  m_axisVCol.clear();
  m_newAxisUCol.clear();
  m_newAxisVCol.clear();

  return StatusCode::SUCCESS;
};



#endif
