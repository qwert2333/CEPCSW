#ifndef _AXISMERGING_ALG_C
#define _AXISMERGING_ALG_C

#include "Algorithm/AxisMergingAlg.h"

StatusCode AxisMergingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_stringPars.find("OutputAxisName")==settings.map_stringPars.end()) settings.map_stringPars["OutputAxisName"] = "MergedAxis";
  if(settings.map_floatPars.find("th_overlap")==settings.map_floatPars.end()) settings.map_floatPars["th_overlap"] = 0.5;
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

    if(m_newAxisUCol.size()<2){ //No need to merge, save into OutputAxis. 
      std::vector<const PandoraPlus::CaloHalfCluster*> tmp_axisCol; tmp_axisCol.clear();
      for(int ic=0; ic<m_newAxisUCol.size(); ic++) tmp_axisCol.push_back(m_newAxisUCol[ic]);
      p_HalfClustersU->at(ih)->setHalfClusters(settings.map_stringPars["OutputAxisName"], tmp_axisCol);
      continue; 
    }


    std::sort( m_newAxisUCol.begin(), m_newAxisUCol.end(), compLayer );
/*
printf("  In HalfClusterU #%d: readin axis size %d \n", ih, m_newAxisUCol.size());
std::map<std::string, std::vector<const PandoraPlus::CaloHalfCluster*> > tmp_HClusMap =  p_HalfClustersU->at(ih)->getHalfClusterMap();
cout<<"Print Readin AxisU: "<<endl;
for(auto iter : tmp_HClusMap){
  cout<<"  Axis name: "<<iter.first<<endl;
  for(int ia=0; ia<iter.second.size(); ia++){
    cout<<"    No. #"<<ia<<endl;
    for(int il=0 ;il<iter.second[ia]->getCluster().size(); il++)
    printf("    LocalMax %d: (%.3f, %.3f, %.3f, %.3f), %p \n", il, iter.second[ia]->getCluster()[il]->getPos().x(), iter.second[ia]->getCluster()[il]->getPos().y(), iter.second[ia]->getCluster()[il]->getPos().z(), iter.second[ia]->getCluster()[il]->getEnergy(), iter.second[ia]->getCluster()[il]);
  }
cout<<endl;
}
*/
//for(int ia=0; ia<m_newAxisUCol.size(); ia++){
//cout<<"  Axis #"<<ia<<endl;
//for(int il=0 ;il<m_newAxisUCol[ia]->getCluster().size(); il++)
//  printf("    LocalMax %d: (%.3f, %.3f, %.3f, %.3f), %p \n", il, m_newAxisUCol[ia]->getCluster()[il]->getPos().x(), m_newAxisUCol[ia]->getCluster()[il]->getPos().y(), m_newAxisUCol[ia]->getCluster()[il]->getPos().z(), m_newAxisUCol[ia]->getCluster()[il]->getEnergy(), m_newAxisUCol[ia]->getCluster()[il]);
//cout<<endl;
//}
    //Case1: Merge axes associated to the same track. 
    TrkMatchedMerging(m_newAxisUCol);    

    //Case1: Merge axes that share same localMax.
    OverlapMerging(m_newAxisUCol);

    //Case2: Merge nearby axes. 
    ConeMerging(m_newAxisUCol);

/*
cout<<"  After merging: axis size "<<m_newAxisUCol.size()<<", Check the overlap"<<endl;
cout<<"Print Merged AxisU: "<<endl;
for(int ia=0; ia<m_newAxisUCol.size(); ia++){
cout<<"  Axis #"<<ia<<endl;
for(int il=0 ;il<m_newAxisUCol[ia]->getCluster().size(); il++)
  printf("    LocalMax %d: (%.3f, %.3f, %.3f, %.3f), towerID [%d, %d, %d], %p \n", il, 
    m_newAxisUCol[ia]->getCluster()[il]->getPos().x(), 
    m_newAxisUCol[ia]->getCluster()[il]->getPos().y(), 
    m_newAxisUCol[ia]->getCluster()[il]->getPos().z(), 
    m_newAxisUCol[ia]->getCluster()[il]->getEnergy(), 
    m_newAxisUCol[ia]->getCluster()[il]->getTowerID()[0][0], 
    m_newAxisUCol[ia]->getCluster()[il]->getTowerID()[0][1], 
    m_newAxisUCol[ia]->getCluster()[il]->getTowerID()[0][2], 
    m_newAxisUCol[ia]->getCluster()[il]);
cout<<endl;
}
*/
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
    std::sort( m_newAxisVCol.begin(), m_newAxisVCol.end(), compLayer );

printf("  In HalfClusterV #%d: readin axis size %d \n", ih, m_newAxisVCol.size());
std::map<std::string, std::vector<const PandoraPlus::CaloHalfCluster*> > tmp_HClusMap =  p_HalfClustersV->at(ih)->getHalfClusterMap();
cout<<"Print Readin AxisV: "<<endl;
for(auto iter : tmp_HClusMap){
  cout<<"  Axis name: "<<iter.first<<endl;
  for(int ia=0; ia<iter.second.size(); ia++){
    cout<<"    No. #"<<ia<<": track size "<<iter.second[ia]->getAssociatedTracks().size()<<endl;
    for(int il=0 ;il<iter.second[ia]->getCluster().size(); il++)
    printf("    LocalMax %d: (%.3f, %.3f, %.3f, %.3f), %p \n", il, iter.second[ia]->getCluster()[il]->getPos().x(), iter.second[ia]->getCluster()[il]->getPos().y(), iter.second[ia]->getCluster()[il]->getPos().z(), iter.second[ia]->getCluster()[il]->getEnergy(), iter.second[ia]->getCluster()[il]);
  }
cout<<endl;
}

    //Case1: Merge axes associated to the same track. 
    TrkMatchedMerging(m_newAxisVCol);
printf("  In HalfClusterV #%d: After Step1: axis size %d \n", ih, m_newAxisVCol.size());

    //Case1: Merge axes that share same localMax.
    OverlapMerging(m_newAxisVCol);
printf("  In HalfClusterV #%d: After Step2: axis size %d \n", ih, m_newAxisVCol.size());

    //Case2: Merge nearby axes.
    ConeMerging(m_newAxisVCol);

/*
cout<<"  After merging: axis size "<<m_newAxisVCol.size()<<", Check the overlap"<<endl;
cout<<"Print Merged AxisV: "<<endl;
for(int ia=0; ia<m_newAxisVCol.size(); ia++){
cout<<"  Axis #"<<ia<<endl;
for(int il=0 ;il<m_newAxisVCol[ia]->getCluster().size(); il++)
  printf("    LocalMax %d: (%.3f, %.3f, %.3f, %.3f), %p \n", il, m_newAxisVCol[ia]->getCluster()[il]->getPos().x(), m_newAxisVCol[ia]->getCluster()[il]->getPos().y(), m_newAxisVCol[ia]->getCluster()[il]->getPos().z(), m_newAxisVCol[ia]->getCluster()[il]->getEnergy(), m_newAxisVCol[ia]->getCluster()[il]);
cout<<endl;
}
*/
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


StatusCode AxisMergingAlg::TrkMatchedMerging( std::vector<PandoraPlus::CaloHalfCluster*>& m_axisCol ){
  if(m_axisCol.size()<2) return StatusCode::SUCCESS;

  for(int iax=0; iax<m_axisCol.size(); iax++){
    std::vector<const PandoraPlus::Track*> m_trkCol = m_axisCol[iax]->getAssociatedTracks();
    for(int jax=iax+1; jax<m_axisCol.size(); jax++){
      std::vector<const PandoraPlus::Track*> p_trkCol = m_axisCol[jax]->getAssociatedTracks();

printf("  In pair (%d, %d): associated trk size (%d, %d) \n", iax, jax, m_trkCol.size(), p_trkCol.size());

      bool fl_match = false;
      for(int itrk=0; itrk<m_trkCol.size(); itrk++){
        if( find(p_trkCol.begin(), p_trkCol.end(), m_trkCol[itrk])!=p_trkCol.end() ){
          fl_match = true; break;
        }
      }
printf("  In pair (%d, %d): match = %d \n", iax, jax, fl_match);

      if(fl_match){
        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );
        delete m_axisCol[jax]; m_axisCol[jax]=NULL;
        m_axisCol.erase(m_axisCol.begin()+jax);
        jax--;
        if(iax>jax+1) iax--;
      }
    }
  }
  return StatusCode::SUCCESS;
};


StatusCode AxisMergingAlg::OverlapMerging( std::vector<PandoraPlus::CaloHalfCluster*>& m_axisCol ){
  if(m_axisCol.size()<2) return StatusCode::SUCCESS;

  for(int iax=0; iax<m_axisCol.size(); iax++){
    const PandoraPlus::CaloHalfCluster* m_axis = m_axisCol[iax];
    for(int jax=iax+1; jax<m_axisCol.size(); jax++){
      const PandoraPlus::CaloHalfCluster* p_axis = m_axisCol[jax];
      std::vector<const Calo1DCluster*> tmp_localMax = p_axis->getCluster();

      int nsharedHits = 0;
      for(int ihit=0; ihit<m_axis->getCluster().size(); ihit++)
        if( find(tmp_localMax.begin(), tmp_localMax.end(), m_axis->getCluster()[ihit])!=tmp_localMax.end() ) nsharedHits++;

printf("  In pair (%d, %d): hit size (%d, %d), shared hit size %d \n",iax, jax, m_axis->getCluster().size(), p_axis->getCluster().size(), nsharedHits);
      
      if( (m_axis->getCluster().size()<=p_axis->getCluster().size() && (float)nsharedHits/m_axis->getCluster().size()>settings.map_floatPars["th_overlap"] ) || 
          (p_axis->getCluster().size()<m_axis->getCluster().size() && (float)nsharedHits/p_axis->getCluster().size()>settings.map_floatPars["th_overlap"] ) ){

cout<<"  Merge: Yes. "<<endl;
        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );

        //if(m_axisCol[iax]->getType()==0 || m_axisCol[jax]->getType()==0)  m_axisCol[iax]->setType(0);
        //else if( m_axisCol[iax]->getType()==1 || m_axisCol[jax]->getType()==1 ) m_axisCol[iax]->setType(1);
        //else m_axisCol[iax]->setType(2);
        m_axisCol[iax]->setType(3);

        delete m_axisCol[jax]; m_axisCol[jax]=NULL;
        m_axisCol.erase(m_axisCol.begin()+jax);
        jax--;
        if(iax>jax+1) iax--;

      }

      p_axis=nullptr;
    }
    m_axis=nullptr;
  }


  //iterate
  for(int iax=0; iax<m_axisCol.size(); iax++){
    const PandoraPlus::CaloHalfCluster* m_axis = m_axisCol[iax];
    for(int jax=iax+1; jax<m_axisCol.size(); jax++){
      const PandoraPlus::CaloHalfCluster* p_axis = m_axisCol[jax];
      std::vector<const Calo1DCluster*> tmp_localMax = p_axis->getCluster();

      int nsharedHits = 0;
      for(int ihit=0; ihit<m_axis->getCluster().size(); ihit++)
        if( find(tmp_localMax.begin(), tmp_localMax.end(), m_axis->getCluster()[ihit])!=tmp_localMax.end() ) nsharedHits++;

printf("  In pair (%d, %d): hit size (%d, %d), shared hit size %d \n",iax, jax, m_axis->getCluster().size(), p_axis->getCluster().size(), nsharedHits);

      if( (m_axis->getCluster().size()<=p_axis->getCluster().size() && (float)nsharedHits/m_axis->getCluster().size()>settings.map_floatPars["th_overlap"] ) ||
          (p_axis->getCluster().size()<m_axis->getCluster().size() && (float)nsharedHits/p_axis->getCluster().size()>settings.map_floatPars["th_overlap"] ) ){

cout<<"  Merge: Yes. "<<endl;
        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );  
        //if(m_axisCol[iax]->getType()==0 || m_axisCol[jax]->getType()==0)  m_axisCol[iax]->setType(0);
        //else if( m_axisCol[iax]->getType()==1 || m_axisCol[jax]->getType()==1 ) m_axisCol[iax]->setType(1);
        //else m_axisCol[iax]->setType(2);
        m_axisCol[iax]->setType(3);

        delete m_axisCol[jax]; m_axisCol[jax]=NULL;
        m_axisCol.erase(m_axisCol.begin()+jax);
        jax--;
        if(iax>jax+1) iax--;

      }

      p_axis=nullptr;
    }
    m_axis=nullptr;
  }


  return StatusCode::SUCCESS;
};


StatusCode AxisMergingAlg::ConeMerging( std::vector<PandoraPlus::CaloHalfCluster*>& m_axisCol ){
  if(m_axisCol.size()<2) return StatusCode::SUCCESS;

  for(int iax=0; iax<m_axisCol.size(); iax++){
    const PandoraPlus::CaloHalfCluster* m_axis = m_axisCol[iax];
    for(int jax=iax+1; jax<m_axisCol.size(); jax++){
      const PandoraPlus::CaloHalfCluster* p_axis = m_axisCol[jax];

      double minRelAngle = min( sin( (m_axis->getPos()-p_axis->getPos()).Angle(m_axis->getAxis())),
                                sin( (m_axis->getPos()-p_axis->getPos()).Angle(p_axis->getAxis())) );
      int skipLayer = m_axis->getBeginningDlayer()<p_axis->getBeginningDlayer() ?
                      (p_axis->getBeginningDlayer()-m_axis->getEndDlayer()) : (m_axis->getBeginningDlayer()-p_axis->getEndDlayer());

      if( sin( m_axis->getAxis().Angle(p_axis->getAxis()) ) < sin(settings.map_floatPars["axis_Angle"]) &&
          sin(minRelAngle) < sin(settings.map_floatPars["relP_Angle"]) &&
          skipLayer <= settings.map_floatPars["skipLayer"] ){

        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );
        m_axisCol[iax]->setType(4);
        delete m_axisCol[jax]; m_axisCol[jax]=NULL;
        m_axisCol.erase(m_axisCol.begin()+jax);
        jax--;
        if(iax>jax+1) iax--;
      }
      p_axis = nullptr;
    }
    m_axis = nullptr;
  }


  //iterate
  for(int iax=0; iax<m_axisCol.size(); iax++){
    const PandoraPlus::CaloHalfCluster* m_axis = m_axisCol[iax];
    for(int jax=iax+1; jax<m_axisCol.size(); jax++){
      const PandoraPlus::CaloHalfCluster* p_axis = m_axisCol[jax];

      double relDis = (m_axis->getEnergyCenter()-p_axis->getEnergyCenter()).Mag();
      double minRelAngle = min( sin( (m_axis->getEnergyCenter()-p_axis->getEnergyCenter()).Angle(m_axis->getAxis())),
                                sin( (m_axis->getEnergyCenter()-p_axis->getEnergyCenter()).Angle(p_axis->getAxis())) );
      int skipLayer = m_axis->getBeginningDlayer()<p_axis->getBeginningDlayer() ?
                      (p_axis->getBeginningDlayer()-m_axis->getEndDlayer()) : (m_axis->getBeginningDlayer()-p_axis->getEndDlayer());

      if( sin( m_axis->getAxis().Angle(p_axis->getAxis()) ) < sin(settings.map_floatPars["axis_Angle"]) &&
          ( sin(minRelAngle) < sin(settings.map_floatPars["relP_Angle"]) )&&
          skipLayer <= settings.map_floatPars["skipLayer"] ){

        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );
        m_axisCol[iax]->setType(2);
        delete m_axisCol[jax]; m_axisCol[jax]=NULL;
        m_axisCol.erase(m_axisCol.begin()+jax);
        jax--;
        if(iax>jax+1) iax--;
      }
      p_axis = nullptr;
    }
    m_axis = nullptr;
  }

  return StatusCode::SUCCESS;
};

#endif
