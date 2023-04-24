#ifndef _AXISMERGING_ALG_C
#define _AXISMERGING_ALG_C

#include "Algorithm/AxisMergingAlg.h"

StatusCode AxisMergingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_stringPars.find("OutputAxisName")==settings.map_stringPars.end()) settings.map_stringPars["OutputAxisName"] = "MergedAxis";
  if(settings.map_floatPars.find("th_overlap")==settings.map_floatPars.end()) settings.map_floatPars["th_overlap"] = 0.5;
  if(settings.map_intPars.find("th_CoreNhit")==settings.map_intPars.end()) settings.map_intPars["th_CoreNhit"] = 7;
  if(settings.map_floatPars.find("axis_Angle")==settings.map_floatPars.end()) settings.map_floatPars["axis_Angle"] = TMath::Pi()/4.;
  if(settings.map_floatPars.find("relP_Angle")==settings.map_floatPars.end()) settings.map_floatPars["relP_Angle"] = TMath::Pi()/4.;
  if(settings.map_floatPars.find("relP_Dis")==settings.map_floatPars.end()) settings.map_floatPars["relP_Dis"] = 5*PandoraPlus::CaloUnit::barsize;
  if(settings.map_intPars.find("th_Nhit")==settings.map_intPars.end()) settings.map_intPars["th_Nhit"] = 5;

  return StatusCode::SUCCESS;
};

StatusCode AxisMergingAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  m_axisUCol.clear();
  m_axisVCol.clear();
  m_newAxisUCol.clear();
  m_newAxisVCol.clear();

  p_HalfClusterU.clear();
  p_HalfClusterV.clear();

  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColU"].size(); ih++)
    p_HalfClusterU.push_back( m_datacol.map_HalfCluster["HalfClusterColU"][ih].get() );
  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColV"].size(); ih++)
    p_HalfClusterV.push_back( m_datacol.map_HalfCluster["HalfClusterColV"][ih].get() );

  return StatusCode::SUCCESS;
};

StatusCode AxisMergingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  //if( !p_HalfClusterU || !p_HalfClusterV ){
  //  std::cout<<"AxisMergingAlg: Input HalfCluster is NULL! "<<std::endl;
  //  return StatusCode::SUCCESS;
  //}
  if( p_HalfClusterU.size() + p_HalfClusterV.size()==0 ) {
    std::cout<<"AxisMergingAlg: No HalfCluster input"<<std::endl;
    return StatusCode::SUCCESS;
  }

cout<<"AxisMergingAlg: Readin halfcluster size: "<<p_HalfClusterU.size()<<", "<<p_HalfClusterV.size()<<endl;


  //Merge axis in HalfClusterU: 
  for(int ih=0; ih<p_HalfClusterU.size(); ih++){
    m_axisUCol.clear(); m_newAxisUCol.clear();
    m_axisUCol = p_HalfClusterU[ih]->getAllHalfClusterCol();

    for(int ic=0; ic<m_axisUCol.size(); ic++){ 
      std::shared_ptr<CaloHalfCluster> ptr_cloned = m_axisUCol[ic]->Clone();
      m_datacol.map_HalfCluster["bkHalfCluster"].push_back(ptr_cloned);
      m_newAxisUCol.push_back( ptr_cloned.get() );
    }

    if(m_newAxisUCol.size()<2){ //No need to merge, save into OutputAxis. 
      std::vector<const PandoraPlus::CaloHalfCluster*> tmp_axisCol; tmp_axisCol.clear();
      for(int ic=0; ic<m_newAxisUCol.size(); ic++) tmp_axisCol.push_back(m_newAxisUCol[ic]);
      p_HalfClusterU[ih]->setHalfClusters(settings.map_stringPars["OutputAxisName"], tmp_axisCol);
      continue; 
    }


    std::sort( m_newAxisUCol.begin(), m_newAxisUCol.end(), compLayer );
/*
printf("  In HalfClusterU #%d: readin axis size %d \n", ih, m_newAxisUCol.size());
std::map<std::string, std::vector<const PandoraPlus::CaloHalfCluster*> > tmp_HClusMap =  p_HalfClusterU->at(ih)->getHalfClusterMap();
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
//printf("  In HalfClusterU #%d: After Step1: axis size %d \n", ih, m_newAxisUCol.size());

    //Case2: Merge axes that share same localMax.
    OverlapMerging(m_newAxisUCol);
//printf("  In HalfClusterU #%d: After Step2: axis size %d \n", ih, m_newAxisUCol.size());

    //Case3: Merge fragments to core axes. 
    FragmentsMerging(m_newAxisUCol);
//printf("  In HalfClusterU #%d: After Step3: axis size %d \n", ih, m_newAxisUCol.size());

    //Case4: Merge nearby axes. 
    //ConeMerging(m_newAxisUCol);
//printf("  In HalfClusterU #%d: After Step4: axis size %d \n", ih, m_newAxisUCol.size());

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
    //Check axis quality
    std::vector<PandoraPlus::CaloHalfCluster*> tmp_goodAxis; tmp_goodAxis.clear(); 
    std::vector<PandoraPlus::CaloHalfCluster*> tmp_badAxis; tmp_badAxis.clear();
    for(int ic=0; ic<m_newAxisUCol.size(); ic++){
      if( (m_newAxisUCol[ic]->getType()==1 && m_newAxisUCol[ic]->getCluster().size()>=settings.map_intPars["th_CoreNhit"] ) ||
          m_newAxisUCol[ic]->getType()==3 || m_newAxisUCol[ic]->getType()==4 )
        tmp_goodAxis.push_back( m_newAxisUCol[ic] );

      else tmp_badAxis.push_back( m_newAxisUCol[ic] );
    }  
printf("  In HalfClusterU: Good axis %d, bad axis %d \n",tmp_goodAxis.size(), tmp_badAxis.size() );

    for(int ic=0; ic<tmp_badAxis.size(); ic++)
      MergeToClosestCluster( tmp_badAxis[ic], tmp_goodAxis );

    //convert to constant object and save in HalfCluster: 
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_axisCol; tmp_axisCol.clear();
    for(int ic=0; ic<tmp_goodAxis.size(); ic++) tmp_axisCol.push_back(tmp_goodAxis[ic]);
    p_HalfClusterU[ih]->setHalfClusters(settings.map_stringPars["OutputAxisName"], tmp_axisCol);

  }


  //Merge axis in HalfClusterV:
  for(int ih=0; ih<p_HalfClusterV.size(); ih++){
    m_axisVCol.clear(); m_newAxisVCol.clear();
    m_axisVCol = p_HalfClusterV[ih]->getAllHalfClusterCol();

    for(int ic=0; ic<m_axisVCol.size(); ic++){
      std::shared_ptr<CaloHalfCluster> ptr_cloned = m_axisVCol[ic]->Clone();
      m_datacol.map_HalfCluster["bkHalfCluster"].push_back(ptr_cloned);
      m_newAxisVCol.push_back( ptr_cloned.get() );
    }


    if(m_newAxisVCol.size()<2){
      std::vector<const PandoraPlus::CaloHalfCluster*> tmp_axisCol; tmp_axisCol.clear();
      for(int ic=0; ic<m_newAxisVCol.size(); ic++) tmp_axisCol.push_back(m_newAxisVCol[ic]);
      p_HalfClusterV[ih]->setHalfClusters(settings.map_stringPars["OutputAxisName"], tmp_axisCol);
      continue;
    }
    std::sort( m_newAxisVCol.begin(), m_newAxisVCol.end(), compLayer );
/*
printf("  In HalfClusterV #%d: readin axis size %d \n", ih, m_newAxisVCol.size());
std::map<std::string, std::vector<const PandoraPlus::CaloHalfCluster*> > tmp_HClusMap =  p_HalfClusterV->at(ih)->getHalfClusterMap();
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
*/
    //Case1: Merge axes associated to the same track. 
    TrkMatchedMerging(m_newAxisVCol);
//printf("  In HalfClusterV #%d: After Step1: axis size %d \n", ih, m_newAxisVCol.size());

    //Case2: Merge axes that share same localMax.
    OverlapMerging(m_newAxisVCol);
//printf("  In HalfClusterV #%d: After Step2: axis size %d \n", ih, m_newAxisVCol.size());

    //Case3: Merge fragments to core axes. 
    FragmentsMerging(m_newAxisVCol);
//printf("  In HalfClusterV #%d: After Step3: axis size %d \n", ih, m_newAxisVCol.size());

    //Case4: Merge nearby axes.
    //ConeMerging(m_newAxisVCol);
//printf("  In HalfClusterV #%d: After Step4: axis size %d \n", ih, m_newAxisVCol.size());

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

    //Check axis quality
    std::vector<PandoraPlus::CaloHalfCluster*> tmp_goodAxis; tmp_goodAxis.clear();
    std::vector<PandoraPlus::CaloHalfCluster*> tmp_badAxis; tmp_badAxis.clear();
    for(int ic=0; ic<m_newAxisVCol.size(); ic++){
      if( (m_newAxisVCol[ic]->getType()==1 && m_newAxisVCol[ic]->getCluster().size()>=settings.map_intPars["th_CoreNhit"] ) ||
          m_newAxisVCol[ic]->getType()==3 || m_newAxisVCol[ic]->getType()==4 )
        tmp_goodAxis.push_back( m_newAxisVCol[ic] );

      else tmp_badAxis.push_back( m_newAxisVCol[ic] );
    }
    for(int ic=0; ic<tmp_badAxis.size(); ic++)
      MergeToClosestCluster( tmp_badAxis[ic], tmp_goodAxis );

    //convert to constant object and save in HalfCluster:
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_axisCol; tmp_axisCol.clear();
    for(int ic=0; ic<tmp_goodAxis.size(); ic++) tmp_axisCol.push_back(tmp_goodAxis[ic]);
    p_HalfClusterV[ih]->setHalfClusters(settings.map_stringPars["OutputAxisName"], tmp_axisCol);

  }

  return StatusCode::SUCCESS;
};

StatusCode AxisMergingAlg::ClearAlgorithm(){
  p_HalfClusterU.clear();
  p_HalfClusterV.clear();
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

//printf("  In pair (%d, %d): associated trk size (%d, %d) \n", iax, jax, m_trkCol.size(), p_trkCol.size());

      bool fl_match = false;
      for(int itrk=0; itrk<m_trkCol.size(); itrk++){
        if( find(p_trkCol.begin(), p_trkCol.end(), m_trkCol[itrk])!=p_trkCol.end() ){
          fl_match = true; break;
        }
      }
//printf("  In pair (%d, %d): match = %d \n", iax, jax, fl_match);

      if(fl_match){
        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );
        //delete m_axisCol[jax]; m_axisCol[jax]=nullptr;
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

//printf("  In pair (%d, %d): hit size (%d, %d), shared hit size %d \n",iax, jax, m_axis->getCluster().size(), p_axis->getCluster().size(), nsharedHits);
      
      if( (m_axis->getCluster().size()<=p_axis->getCluster().size() && (float)nsharedHits/m_axis->getCluster().size()>settings.map_floatPars["th_overlap"] ) || 
          (p_axis->getCluster().size()<m_axis->getCluster().size() && (float)nsharedHits/p_axis->getCluster().size()>settings.map_floatPars["th_overlap"] ) ){

//cout<<"  Merge: Yes. "<<endl;
        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );

        if(m_axisCol[iax]->getType()==1 || m_axisCol[jax]->getType()==1)  m_axisCol[iax]->setType(3);
        else m_axisCol[iax]->setType(4);

        //delete m_axisCol[jax]; m_axisCol[jax]=nullptr;
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

//printf("  In pair (%d, %d): hit size (%d, %d), shared hit size %d \n",iax, jax, m_axis->getCluster().size(), p_axis->getCluster().size(), nsharedHits);

      if( (m_axis->getCluster().size()<=p_axis->getCluster().size() && (float)nsharedHits/m_axis->getCluster().size()>settings.map_floatPars["th_overlap"] ) ||
          (p_axis->getCluster().size()<m_axis->getCluster().size() && (float)nsharedHits/p_axis->getCluster().size()>settings.map_floatPars["th_overlap"] ) ){

//cout<<"  Merge: Yes. "<<endl;
        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );  
        if(m_axisCol[iax]->getType()==1 || m_axisCol[jax]->getType()==1 || m_axisCol[iax]->getType()==3 || m_axisCol[jax]->getType()==3)  m_axisCol[iax]->setType(3);
        else m_axisCol[iax]->setType(4);

        //delete m_axisCol[jax]; m_axisCol[jax]=nullptr;
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


StatusCode AxisMergingAlg::FragmentsMerging( std::vector<PandoraPlus::CaloHalfCluster*>& m_axisCol ){
  if(m_axisCol.size()<2) return StatusCode::SUCCESS;

/*
cout<<"FragmentsMerging: print readin axes "<<endl;
for(int ia=0; ia<m_axisCol.size(); ia++){
cout<<"  Axis #"<<ia<<", Energy center ";
printf("(%.3f, %.3f, %.3f, %.3f) \n", m_axisCol[ia]->getEnergyCenter().x(), m_axisCol[ia]->getEnergyCenter().y(), m_axisCol[ia]->getEnergyCenter().z());
for(int il=0 ;il<m_axisCol[ia]->getCluster().size(); il++)
  printf("    LocalMax %d: (%.3f, %.3f, %.3f, %.3f) \n", il,
    m_axisCol[ia]->getCluster()[il]->getPos().x(),
    m_axisCol[ia]->getCluster()[il]->getPos().y(),
    m_axisCol[ia]->getCluster()[il]->getPos().z(),
    m_axisCol[ia]->getCluster()[il]->getEnergy() );
cout<<endl;
}
*/

  //Merge fragments to core.
  std::sort( m_axisCol.begin(), m_axisCol.end(), compLayer );
  for(int iax=0; iax<m_axisCol.size(); iax++){
    const PandoraPlus::CaloHalfCluster* m_axis = m_axisCol[iax];
    //Define the Core: Hough+Nhit || Trk+Hough || Trk+Cone.
//cout<<"  Readin core axis #"<<iax<<" type: "<<m_axis->getType()<<", Nclus: "<<m_axis->getCluster().size();
    if( !( (m_axis->getType()==1 && m_axis->getCluster().size()>=settings.map_intPars["th_CoreNhit"] ) ||
        m_axis->getType()==3 || m_axis->getType()==4 ) ) 
      { m_axis=nullptr; continue; }

    for(int jax=0; jax<m_axisCol.size(); jax++){
      if(jax==iax) continue;
      const PandoraPlus::CaloHalfCluster* p_axis = m_axisCol[jax];
      //Do not merge 2 cores. 
      if( (p_axis->getType()==1 && p_axis->getCluster().size()>=settings.map_intPars["th_CoreNhit"] ) ||
          p_axis->getType()==3 || p_axis->getType()==4 )
        {p_axis=nullptr; continue; }

//printf("    Check Core #%d and Fragment #%d: frag axis: (%.3f, %.3f, %.3f) \n", iax, jax, p_axis->getAxis().x(), p_axis->getAxis().y(), p_axis->getAxis().z());

      TVector3 relP1 = p_axis->getClusterInLayer( p_axis->getBeginningDlayer() )[0]->getPos() - m_axis->getClusterInLayer( m_axis->getEndDlayer() )[0]->getPos();
      TVector3 relP2 = p_axis->getClusterInLayer( p_axis->getBeginningDlayer() )[0]->getPos() - m_axis->getEnergyCenter();
      TVector3 relP3 = m_axis->getEnergyCenter() - p_axis->getClusterInLayer( p_axis->getEndDlayer() )[0]->getPos();
      TVector3 relP4 = p_axis->getPos() - m_axis->getPos();

//printf("      R[head-end]: (%.3f, %.3f, %.3f), R[head-cent]: (%.3f, %.3f, %.3f), R[tail-cent]: (%.3f, %.3f, %.3f) \n", 
//relP1.x(), relP1.y(), relP1.z(),
//relP2.x(), relP2.y(), relP2.z(),
//relP3.x(), relP3.y(), relP3.z() );

      if( (relP1.Angle(p_axis->getAxis()) < settings.map_floatPars["axis_Angle"] && relP1.Mag()<=settings.map_floatPars["relP_Dis"]) ||
          (relP2.Angle(p_axis->getAxis()) < settings.map_floatPars["axis_Angle"] && relP2.Mag()<=settings.map_floatPars["relP_Dis"]) ||
          (relP3.Angle(p_axis->getAxis()) < settings.map_floatPars["axis_Angle"] && relP3.Mag()<=settings.map_floatPars["relP_Dis"]) || 
          ( p_axis->getType()==2 && relP4.Angle(m_axis->getAxis())<settings.map_floatPars["axis_Angle"] && relP4.Mag()<=settings.map_floatPars["relP_Dis"]) )
        {
//cout<<"      Merge!"<<endl;
 
        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );
        //delete m_axisCol[jax]; m_axisCol[jax]=nullptr;
        m_axisCol.erase(m_axisCol.begin()+jax);
        jax--;
        if(iax>jax+1) iax--;
      }

      p_axis=nullptr;
    }
    m_axis = nullptr;
  }


  return StatusCode::SUCCESS;
};


StatusCode AxisMergingAlg::ConeMerging( std::vector<PandoraPlus::CaloHalfCluster*>& m_axisCol ){
  if(m_axisCol.size()<2) return StatusCode::SUCCESS;

  std::sort( m_axisCol.begin(), m_axisCol.end(), compLayer );
  for(int iax=0; iax<m_axisCol.size(); iax++){
    const PandoraPlus::CaloHalfCluster* m_axis = m_axisCol[iax];
    for(int jax=iax+1; jax<m_axisCol.size(); jax++){
      const PandoraPlus::CaloHalfCluster* p_axis = m_axisCol[jax];

//printf("  Layer range: axis %d: [%d, %d], axis %d: [%d, %d] \n", iax, m_axis->getBeginningDlayer(), m_axis->getEndDlayer(), jax, p_axis->getBeginningDlayer(), p_axis->getEndDlayer());

      TVector3 m_beginPoint = m_axis->getClusterInLayer( m_axis->getBeginningDlayer() )[0]->getPos();
      TVector3 m_endPoint   = m_axis->getClusterInLayer( m_axis->getEndDlayer() )[0]->getPos();
      TVector3 p_beginPoint = p_axis->getClusterInLayer( p_axis->getBeginningDlayer() )[0]->getPos();
      TVector3 p_endPoint   = p_axis->getClusterInLayer( p_axis->getEndDlayer() )[0]->getPos();
      double relDis = -1.; 
      if( m_endPoint.Mag()<p_beginPoint.Mag() )      relDis = p_beginPoint.Mag() - m_endPoint.Mag(); 
      else if( p_endPoint.Mag()<m_beginPoint.Mag() ) relDis = m_beginPoint.Mag() - p_endPoint.Mag();
      else if( m_beginPoint.Mag()>p_beginPoint.Mag() && m_endPoint.Mag()>p_endPoint.Mag() )  relDis = p_endPoint.Mag() - m_beginPoint.Mag();
      else if( p_beginPoint.Mag()>m_beginPoint.Mag() && p_endPoint.Mag()>m_endPoint.Mag() )  relDis = m_endPoint.Mag() - p_beginPoint.Mag();
      else if( (m_beginPoint.Mag()>p_beginPoint.Mag() && m_endPoint.Mag()<p_endPoint.Mag()) ||
               (m_beginPoint.Mag()<p_beginPoint.Mag() && m_endPoint.Mag()>p_endPoint.Mag()) ) relDis = 9999.;

      double minRelAngle = min( sin( (m_axis->getPos()-p_axis->getPos()).Angle(m_axis->getAxis())),
                                sin( (m_axis->getPos()-p_axis->getPos()).Angle(p_axis->getAxis())) );
//      int skipLayer = m_axis->getBeginningDlayer()<p_axis->getBeginningDlayer() ?
//                      (p_axis->getBeginningDlayer()-m_axis->getEndDlayer()) : (m_axis->getBeginningDlayer()-p_axis->getEndDlayer());

//printf("  AxisMerging: axes pair (%d, %d): axisAngel = %.3f, RelPAngle = %.3f, dis = %d \n", iax, jax, m_axis->getAxis().Angle(p_axis->getAxis()), minRelAngle, relDis);


      if( ( sin( m_axis->getAxis().Angle(p_axis->getAxis()) ) < sin(settings.map_floatPars["axis_Angle"]) ||
            sin(minRelAngle) < sin(settings.map_floatPars["relP_Angle"]) )  &&
          relDis <= settings.map_floatPars["relP_Dis"] && relDis>=0 ){
//cout<<"    Merge! "<<endl;

        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );
        if(m_axisCol[iax]->getType()==1 || m_axisCol[jax]->getType()==1 || m_axisCol[iax]->getType()==5 || m_axisCol[jax]->getType()==5) m_axisCol[iax]->setType(5);
        else if(m_axisCol[iax]->getType()==3 || m_axisCol[jax]->getType()==3 || m_axisCol[iax]->getType()==6 || m_axisCol[jax]->getType()==6) m_axisCol[iax]->setType(6);
        else if(m_axisCol[iax]->getType()==4 || m_axisCol[jax]->getType()==4) m_axisCol[iax]->setType(4);
        else m_axisCol[iax]->setType(7);

        //delete m_axisCol[jax]; m_axisCol[jax]=nullptr;
        m_axisCol.erase(m_axisCol.begin()+jax);
        jax--;
        if(iax>jax+1) iax--;
      }
      p_axis = nullptr;
    }
    m_axis = nullptr;
  }

//cout<<"Iteration merge"<<endl;

  //iterate
  for(int iax=0; iax<m_axisCol.size(); iax++){
    const PandoraPlus::CaloHalfCluster* m_axis = m_axisCol[iax];
    for(int jax=iax+1; jax<m_axisCol.size(); jax++){
      const PandoraPlus::CaloHalfCluster* p_axis = m_axisCol[jax];

      TVector3 m_beginPoint = m_axis->getClusterInLayer( m_axis->getBeginningDlayer() )[0]->getPos();
      TVector3 m_endPoint   = m_axis->getClusterInLayer( m_axis->getEndDlayer() )[0]->getPos();
      TVector3 p_beginPoint = p_axis->getClusterInLayer( p_axis->getBeginningDlayer() )[0]->getPos();
      TVector3 p_endPoint   = p_axis->getClusterInLayer( p_axis->getEndDlayer() )[0]->getPos();
      double relDis = -1.;
      if( m_endPoint.Mag()<p_beginPoint.Mag() )      relDis = p_beginPoint.Mag() - m_endPoint.Mag();
      else if( p_endPoint.Mag()<m_beginPoint.Mag() ) relDis = m_beginPoint.Mag() - p_endPoint.Mag();
      else if( m_beginPoint.Mag()>p_beginPoint.Mag() && m_endPoint.Mag()>p_endPoint.Mag() )  relDis = p_endPoint.Mag() - m_beginPoint.Mag();
      else if( p_beginPoint.Mag()>m_beginPoint.Mag() && p_endPoint.Mag()>m_endPoint.Mag() )  relDis = m_endPoint.Mag() - p_beginPoint.Mag();
      else if( (m_beginPoint.Mag()>p_beginPoint.Mag() && m_endPoint.Mag()<p_endPoint.Mag()) ||
               (m_beginPoint.Mag()<p_beginPoint.Mag() && m_endPoint.Mag()>p_endPoint.Mag()) ) relDis = 9999.;
      double minRelAngle = min( sin( (m_axis->getEnergyCenter()-p_axis->getEnergyCenter()).Angle(m_axis->getAxis())),
                                sin( (m_axis->getEnergyCenter()-p_axis->getEnergyCenter()).Angle(p_axis->getAxis())) );
      //int skipLayer = m_axis->getBeginningDlayer()<p_axis->getBeginningDlayer() ?
      //                (p_axis->getBeginningDlayer()-m_axis->getEndDlayer()) : (m_axis->getBeginningDlayer()-p_axis->getEndDlayer());

//printf("  AxisMerging: axes pair (%d, %d): axisAngel = %.3f, RelPAngle = %.3f, dis = %d \n", iax, jax, m_axis->getAxis().Angle(p_axis->getAxis()), minRelAngle, relDis);
//printf("  Layer range: axis %d: [%d, %d], axis %d: [%d, %d] \n", iax, m_axis->getBeginningDlayer(), m_axis->getEndDlayer(), jax, p_axis->getBeginningDlayer(), p_axis->getEndDlayer());

      if( ( sin( m_axis->getAxis().Angle(p_axis->getAxis()) ) < sin(settings.map_floatPars["axis_Angle"]) ||
            sin(minRelAngle) < sin(settings.map_floatPars["relP_Angle"]) ) &&
          relDis <= settings.map_floatPars["relP_Dis"] && relDis>=0 ){

        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );
        if(m_axisCol[iax]->getType()==1 || m_axisCol[jax]->getType()==1 || m_axisCol[iax]->getType()==5 || m_axisCol[jax]->getType()==5) m_axisCol[iax]->setType(5);
        else if(m_axisCol[iax]->getType()==3 || m_axisCol[jax]->getType()==3 || m_axisCol[iax]->getType()==6 || m_axisCol[jax]->getType()==6) m_axisCol[iax]->setType(6);
        else if(m_axisCol[iax]->getType()==4 || m_axisCol[jax]->getType()==4) m_axisCol[iax]->setType(4);
        else m_axisCol[iax]->setType(7);

        //delete m_axisCol[jax]; m_axisCol[jax]=nullptr;
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


StatusCode AxisMergingAlg::MergeToClosestCluster( PandoraPlus::CaloHalfCluster* m_badaxis, std::vector<PandoraPlus::CaloHalfCluster*>& m_axisCol ){
  if(m_axisCol.size()==0){ 
    m_axisCol.push_back(m_badaxis);
    return StatusCode::SUCCESS;
  }

  float minR = 9999;
  int index = -1;
  for(int iax=0; iax<m_axisCol.size(); iax++){
    float tmp_R = (m_badaxis->getEnergyCenter()-m_axisCol[iax]->getEnergyCenter()).Mag();
    if( tmp_R<minR ){
      minR = tmp_R;
      index = iax;
    }
  }
  if(index>=0) m_axisCol[index]->mergeHalfCluster(m_badaxis);

  //delete m_badaxis; m_badaxis=nullptr;
  return StatusCode::SUCCESS;
};

#endif
