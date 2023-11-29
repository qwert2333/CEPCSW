#ifndef _AXISMERGING_ALG_C
#define _AXISMERGING_ALG_C

#include "Algorithm/AxisMergingAlg.h"

StatusCode AxisMergingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_stringPars.find("OutputAxisName")==settings.map_stringPars.end()) settings.map_stringPars["OutputAxisName"] = "MergedAxis";
  if(settings.map_floatPars.find("th_overlap")==settings.map_floatPars.end()) settings.map_floatPars["th_overlap"] = 0.5;
  if(settings.map_intPars.find("th_CoreNhit")==settings.map_intPars.end()) settings.map_intPars["th_CoreNhit"] = 3;
  if(settings.map_floatPars.find("axis_Angle")==settings.map_floatPars.end()) settings.map_floatPars["axis_Angle"] = TMath::Pi()/4.;
  if(settings.map_floatPars.find("relP_Angle")==settings.map_floatPars.end()) settings.map_floatPars["relP_Angle"] = TMath::Pi()/4.;
  if(settings.map_floatPars.find("relP_Dis")==settings.map_floatPars.end()) settings.map_floatPars["relP_Dis"] = 5*PandoraPlus::CaloUnit::barsize;
  if(settings.map_intPars.find("th_Nhit")==settings.map_intPars.end()) settings.map_intPars["th_Nhit"] = 5;
  if(settings.map_floatPars.find("th_branch_distance")==settings.map_floatPars.end()) settings.map_floatPars["th_branch_distance"] = 30;

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
  //cout << "yyy: ----------------------------  Running Axis MergingAlg  -----------------------------------" << endl;

  if( p_HalfClusterU.size() + p_HalfClusterV.size()==0 ) {
    std::cout<<"AxisMergingAlg: No HalfCluster input"<<std::endl;
    return StatusCode::SUCCESS;
  }

  // cout << "yyy: Readin halfcluster size in U, V: "<<p_HalfClusterU.size()<<", "<<p_HalfClusterV.size()<<endl;

  // cout << "yyy: Merge axis in U pHalfClusterU" << endl;
  //Merge axis in HalfClusterU: 
  for(int ih=0; ih<p_HalfClusterU.size(); ih++){
    m_axisUCol.clear(); m_newAxisUCol.clear();
    m_axisUCol = p_HalfClusterU[ih]->getAllHalfClusterCol();

    for(int ic=0; ic<m_axisUCol.size(); ic++){ 
      std::shared_ptr<CaloHalfCluster> ptr_cloned = m_axisUCol[ic]->Clone();
      m_datacol.map_HalfCluster["bkHalfCluster"].push_back(ptr_cloned);
      m_newAxisUCol.push_back( ptr_cloned.get() );
    }

    // cout << "  yyy: For p_HalfClusterU[" << ih << "], m_newAxisUCol.size() = " << m_newAxisUCol.size() 
    //      << ", if size<2, no need to merge" << endl;

    if(m_newAxisUCol.size()<2){ //No need to merge, save into OutputAxis. 
      std::vector<const PandoraPlus::CaloHalfCluster*> tmp_axisCol; tmp_axisCol.clear();
      for(int ic=0; ic<m_newAxisUCol.size(); ic++) tmp_axisCol.push_back(m_newAxisUCol[ic]);
      p_HalfClusterU[ih]->setHalfClusters(settings.map_stringPars["OutputAxisName"], tmp_axisCol);
      continue; 
    }


    std::sort( m_newAxisUCol.begin(), m_newAxisUCol.end(), compLayer );


//printf("  In HalfClusterU #%d: readin axis size %d \n", ih, m_newAxisUCol.size());
/*
std::map<std::string, std::vector<const PandoraPlus::CaloHalfCluster*> > tmp_HClusMap =  p_HalfClusterU.at(ih)->getHalfClusterMap();
cout<<"Print Readin AxisU: "<<endl;
for(auto iter : tmp_HClusMap){
  cout<<"  Axis name: "<<iter.first<<endl;
  for(int ia=0; ia<iter.second.size(); ia++){
    cout<<"    No. #"<<ia<<endl;
    for(int il=0 ;il<iter.second[ia]->getCluster().size(); il++)
    printf("    LocalMax %d: (%.3f, %.3f, %.3f, %.3f), towerID [%d, %d, %d], %p \n", il, 
        iter.second[ia]->getCluster()[il]->getPos().x(), 
        iter.second[ia]->getCluster()[il]->getPos().y(), 
        iter.second[ia]->getCluster()[il]->getPos().z(), 
        iter.second[ia]->getCluster()[il]->getEnergy(), 
        iter.second[ia]->getCluster()[il]);
  }
cout<<endl;
}
*/
/*
cout<<"Print Readin AxisU: "<<endl;
for(int ia=0; ia<m_newAxisUCol.size(); ia++){
cout<<"  Axis #"<<ia<<": type "<<m_newAxisUCol[ia]->getType()<<endl;
for(int il=0 ;il<m_newAxisUCol[ia]->getCluster().size(); il++)
  printf("    LocalMax %d: (%.3f, %.3f, %.3f, %.3f), towerID [%d, %d, %d], %p \n", il, 
    m_newAxisUCol[ia]->getCluster()[il]->getPos().x(), 
    m_newAxisUCol[ia]->getCluster()[il]->getPos().y(), 
    m_newAxisUCol[ia]->getCluster()[il]->getPos().z(), 
    m_newAxisUCol[ia]->getCluster()[il]->getTowerID()[0][0],
    m_newAxisUCol[ia]->getCluster()[il]->getTowerID()[0][1],
    m_newAxisUCol[ia]->getCluster()[il]->getTowerID()[0][2],
    m_newAxisUCol[ia]->getCluster()[il]->getEnergy(), 
    m_newAxisUCol[ia]->getCluster()[il]);
cout<<endl;
}
*/
    //Case1: Merge axes associated to the same track. 
    TrkMatchedMerging(m_newAxisUCol);    
    // cout << "  yyy: after TrkMatchedMerging(), axis size = " <<  m_newAxisUCol.size() << endl;

    //Case2: Merge axes that share same localMax.
    OverlapMerging(m_newAxisUCol);
    // cout << "  yyy: after OverlapMerging(), axis size = " <<  m_newAxisUCol.size() << endl;

    // Case3: Merge fake photon to track axis.
    BranchMerging(m_newAxisUCol);
    // cout << "  yyy: after BranchMerging(), axis size = " <<  m_newAxisUCol.size() << endl;

    //Case4: Merge fragments to core axes. 
    FragmentsMerging(m_newAxisUCol);
    // cout << "  yyy: after FragmentsMerging(), axis size = " <<  m_newAxisUCol.size() << endl;

    //Case4: Merge nearby axes. 
    //ConeMerging(m_newAxisUCol);
    //printf("  In HalfClusterU #%d: After Step4: axis size %d \n", ih, m_newAxisUCol.size());


    // cout<<"  yyy: after all merging: axis size "<<m_newAxisUCol.size()<<endl;
    // for(int ia=0; ia<m_newAxisUCol.size(); ia++){
    //   cout<<"  yyy: axis #"<<ia<<endl;
    //   for(int il=0 ;il<m_newAxisUCol[ia]->getCluster().size(); il++)
    //     printf("         LocalMax %d: (%.2f, %.2f, %.2f, %.4f), towerID [%d, %d, %d], %p \n", il, 
    //       m_newAxisUCol[ia]->getCluster()[il]->getPos().x(), 
    //       m_newAxisUCol[ia]->getCluster()[il]->getPos().y(), 
    //       m_newAxisUCol[ia]->getCluster()[il]->getPos().z(), 
    //       m_newAxisUCol[ia]->getCluster()[il]->getEnergy(), 
    //       m_newAxisUCol[ia]->getCluster()[il]->getTowerID()[0][0], 
    //       m_newAxisUCol[ia]->getCluster()[il]->getTowerID()[0][1], 
    //       m_newAxisUCol[ia]->getCluster()[il]->getTowerID()[0][2], 
    //       m_newAxisUCol[ia]->getCluster()[il]);
    // }


    // cout << "  yyy: check axis quality" << endl;
//printf("  In HalfClusterU #%d: After Step4: axis size %d \n", ih, m_newAxisUCol.size());

/*
cout<<"  After merging: axis size "<<m_newAxisUCol.size()<<", Check the overlap"<<endl;
cout<<"Print Merged AxisU: "<<endl;
for(int ia=0; ia<m_newAxisUCol.size(); ia++){
cout<<"  Axis #"<<ia<<": type "<<m_newAxisUCol[ia]->getType()<<endl;
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
      if( (m_newAxisUCol[ic]->getType()==100 && m_newAxisUCol[ic]->getCluster().size()>=settings.map_intPars["th_CoreNhit"] ) ||
          m_newAxisUCol[ic]->getType()>100 ){
            tmp_goodAxis.push_back( m_newAxisUCol[ic] );
          }
        

      else {
        tmp_badAxis.push_back( m_newAxisUCol[ic] );
      }
    }  

    for(int ic=0; ic<tmp_badAxis.size(); ic++){
      MergeToClosestCluster( tmp_badAxis[ic], tmp_goodAxis );
    }
      

    //convert to constant object and save in HalfCluster: 
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_axisCol; tmp_axisCol.clear();
    for(int ic=0; ic<tmp_goodAxis.size(); ic++) tmp_axisCol.push_back(tmp_goodAxis[ic]);
    // for(int ic=0; ic<m_newAxisUCol.size(); ic++) tmp_axisCol.push_back(m_newAxisUCol[ic]);
    p_HalfClusterU[ih]->setHalfClusters(settings.map_stringPars["OutputAxisName"], tmp_axisCol);

  }

  // cout << "yyy: Merge axis in V pHalfClusterV" << endl;
  //Merge axis in HalfClusterV:
  for(int ih=0; ih<p_HalfClusterV.size(); ih++){
    m_axisVCol.clear(); m_newAxisVCol.clear();
    m_axisVCol = p_HalfClusterV[ih]->getAllHalfClusterCol();

    for(int ic=0; ic<m_axisVCol.size(); ic++){
      std::shared_ptr<CaloHalfCluster> ptr_cloned = m_axisVCol[ic]->Clone();
      m_datacol.map_HalfCluster["bkHalfCluster"].push_back(ptr_cloned);
      m_newAxisVCol.push_back( ptr_cloned.get() );
    }

    // cout << "  yyy: For p_HalfClusterV[" << ih << "], m_newAxisVCol.size() = " << m_newAxisVCol.size() 
    //      << ", if size<2, no need to merge" << endl;

    if(m_newAxisVCol.size()<2){
      std::vector<const PandoraPlus::CaloHalfCluster*> tmp_axisCol; tmp_axisCol.clear();
      for(int ic=0; ic<m_newAxisVCol.size(); ic++) tmp_axisCol.push_back(m_newAxisVCol[ic]);
      p_HalfClusterV[ih]->setHalfClusters(settings.map_stringPars["OutputAxisName"], tmp_axisCol);
      continue;
    }
    std::sort( m_newAxisVCol.begin(), m_newAxisVCol.end(), compLayer );
    //Case1: Merge axes associated to the same track. 
    TrkMatchedMerging(m_newAxisVCol);
    // cout << "  yyy: after TrkMatchedMerging(), axis size = " <<  m_newAxisVCol.size() << endl;

    //Case2: Merge axes that share same localMax.
    OverlapMerging(m_newAxisVCol);
    // cout << "  yyy: after OverlapMerging(), axis size = " <<  m_newAxisVCol.size() << endl;

    //Case3: Merge fake photon to track axis.
    BranchMerging(m_newAxisVCol);
    // cout << "  yyy: after BranchMerging(), axis size = " <<  m_newAxisVCol.size() << endl;

    //Case4: Merge fragments to core axes. 
    FragmentsMerging(m_newAxisVCol);
    // cout << "  yyy: after FragmentsMerging(), axis size = " <<  m_newAxisVCol.size() << endl;

    //Case4: Merge nearby axes.
    //ConeMerging(m_newAxisVCol);
    //printf("  In HalfClusterV #%d: After Step4: axis size %d \n", ih, m_newAxisVCol.size());

    // cout<<"  yyy: after all merging: axis size "<<m_newAxisVCol.size()<<endl;
    // for(int ia=0; ia<m_newAxisVCol.size(); ia++){
    //   cout<<"  yyy: axis #"<<ia<<endl;
    //   for(int il=0 ;il<m_newAxisVCol[ia]->getCluster().size(); il++)
    //     printf("         LocalMax %d: (%.2f, %.2f, %.2f, %.4f), towerID [%d, %d, %d], %p \n", il, 
    //       m_newAxisVCol[ia]->getCluster()[il]->getPos().x(), 
    //       m_newAxisVCol[ia]->getCluster()[il]->getPos().y(), 
    //       m_newAxisVCol[ia]->getCluster()[il]->getPos().z(), 
    //       m_newAxisVCol[ia]->getCluster()[il]->getEnergy(), 
    //       m_newAxisVCol[ia]->getCluster()[il]->getTowerID()[0][0], 
    //       m_newAxisVCol[ia]->getCluster()[il]->getTowerID()[0][1], 
    //       m_newAxisVCol[ia]->getCluster()[il]->getTowerID()[0][2], 
    //       m_newAxisVCol[ia]->getCluster()[il]);
    // }

    // cout << "  yyy: check axis quality" << endl;
    //Check axis quality
    std::vector<PandoraPlus::CaloHalfCluster*> tmp_goodAxis; tmp_goodAxis.clear();
    std::vector<PandoraPlus::CaloHalfCluster*> tmp_badAxis; tmp_badAxis.clear();
    for(int ic=0; ic<m_newAxisVCol.size(); ic++){
      if( (m_newAxisVCol[ic]->getType()==100 && m_newAxisVCol[ic]->getCluster().size()>=settings.map_intPars["th_CoreNhit"] ) ||
          m_newAxisVCol[ic]->getType()>100 ){
            tmp_goodAxis.push_back( m_newAxisVCol[ic] );
          }
        

      else {
        tmp_badAxis.push_back( m_newAxisVCol[ic] );
      }
    }

    for(int ic=0; ic<tmp_badAxis.size(); ic++){
      MergeToClosestCluster( tmp_badAxis[ic], tmp_goodAxis );
    }
      

    //convert to constant object and save in HalfCluster:
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_axisCol; tmp_axisCol.clear();
    for(int ic=0; ic<tmp_goodAxis.size(); ic++) tmp_axisCol.push_back(tmp_goodAxis[ic]);
    // for(int ic=0; ic<m_newAxisVCol.size(); ic++) tmp_axisCol.push_back(m_newAxisVCol[ic]);
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
  // cout << "  yyy: calling TrkMatchedMerging(), m_axisCol.size()=" << m_axisCol.size() << endl;
  if(m_axisCol.size()<2) return StatusCode::SUCCESS;

  for(int iax=0; iax<m_axisCol.size(); iax++){
    std::vector<const PandoraPlus::Track*> m_trkCol = m_axisCol[iax]->getAssociatedTracks();
    for(int jax=iax+1; jax<m_axisCol.size(); jax++){
      std::vector<const PandoraPlus::Track*> p_trkCol = m_axisCol[jax]->getAssociatedTracks();


      bool fl_match = false;
      for(int itrk=0; itrk<m_trkCol.size(); itrk++){
        if( find(p_trkCol.begin(), p_trkCol.end(), m_trkCol[itrk])!=p_trkCol.end() ){
          fl_match = true; break;
        }
      }

      if(fl_match){
        // cout << "    yyy: m_axisCol[" << iax << "] with type=" << m_axisCol[iax]->getType() 
        //      << " and [" << jax << "] with type=" << m_axisCol[jax]->getType() << "share the same associated track, merge them" << endl;
        // cout << "         hit in m_axisCol["<<iax<<"]:"<<endl;
        // for(int yii=0; yii<m_axisCol[iax]->getCluster().size(); yii++){
        //   cout << "           (" << std::fixed << std::setprecision(2) << m_axisCol[iax]->getCluster()[yii]->getPos().x()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[iax]->getCluster()[yii]->getPos().y()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[iax]->getCluster()[yii]->getPos().z()
        //        << ")" << endl;
        // }
        // cout << "         hit in m_axisCol["<<jax<<"]:"<<endl;
        // for(int yii=0; yii<m_axisCol[jax]->getCluster().size(); yii++){
        //   cout << "           (" << std::fixed << std::setprecision(2) << m_axisCol[jax]->getCluster()[yii]->getPos().x()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[jax]->getCluster()[yii]->getPos().y()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[jax]->getCluster()[yii]->getPos().z()
        //        << ")" << endl;
        // }

        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );
        // cout << "    yyy: after merge them, new type=" << m_axisCol[iax]->getType() << endl;
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
  // cout << "  yyy: calling OverlapMerging(), m_axisCol.size()=" << m_axisCol.size() << endl;
  if(m_axisCol.size()<2) return StatusCode::SUCCESS;

  // cout << "  yyy: first interate" << endl;
  for(int iax=0; iax<m_axisCol.size(); iax++){
    const PandoraPlus::CaloHalfCluster* m_axis = m_axisCol[iax];
    for(int jax=iax+1; jax<m_axisCol.size(); jax++){
      const PandoraPlus::CaloHalfCluster* p_axis = m_axisCol[jax];

      // Do not merge two different track axes
      if ( m_axis->getAssociatedTracks().size()>0 && p_axis->getAssociatedTracks().size()>0 ) continue;

      std::vector<const Calo1DCluster*> tmp_localMax = p_axis->getCluster();

      int nsharedHits = 0;
      for(int ihit=0; ihit<m_axis->getCluster().size(); ihit++)
        if( find(tmp_localMax.begin(), tmp_localMax.end(), m_axis->getCluster()[ihit])!=tmp_localMax.end() ) nsharedHits++;

      
      if( (m_axis->getCluster().size()<=p_axis->getCluster().size() && (float)nsharedHits/m_axis->getCluster().size()>settings.map_floatPars["th_overlap"] ) || 
          (p_axis->getCluster().size()<m_axis->getCluster().size() && (float)nsharedHits/p_axis->getCluster().size()>settings.map_floatPars["th_overlap"] ) ){

        // cout << "    yyy: m_axisCol[" << iax << "] with type=" << m_axisCol[iax]->getType() 
        //      << " and [" << jax << "] with type=" << m_axisCol[jax]->getType() 
        //      << " share " << nsharedHits << " hits, merge them" << endl;
        // cout << "         hit in m_axisCol["<<iax<<"]:"<<endl;
        // for(int yii=0; yii<m_axisCol[iax]->getCluster().size(); yii++){
        //   cout << "           (" << std::fixed << std::setprecision(2) << m_axisCol[iax]->getCluster()[yii]->getPos().x()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[iax]->getCluster()[yii]->getPos().y()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[iax]->getCluster()[yii]->getPos().z()
        //        << ")" << endl;
        // }
        // cout << "         hit in m_axisCol["<<jax<<"]:"<<endl;
        // for(int yii=0; yii<m_axisCol[jax]->getCluster().size(); yii++){
        //   cout << "           (" << std::fixed << std::setprecision(2) << m_axisCol[jax]->getCluster()[yii]->getPos().x()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[jax]->getCluster()[yii]->getPos().y()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[jax]->getCluster()[yii]->getPos().z()
        //        << ")" << endl;
        // }

        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );

        // track axis + neutral axis = track axis
        int axis_type = m_axisCol[iax]->getType() + m_axisCol[jax]->getType();
        m_axisCol[iax]->setType(axis_type);

        // cout << "    yyy: after merge them, new type=" << m_axisCol[iax]->getType() << endl;
        //delete m_axisCol[jax]; m_axisCol[jax]=nullptr;
        m_axisCol.erase(m_axisCol.begin()+jax);
        jax--;
        if(iax>jax+1) iax--;

      }

      p_axis=nullptr;
    }
    m_axis=nullptr;
  }


  // cout << "  yyy: second interate" << endl;
  //iterate
  for(int iax=0; iax<m_axisCol.size(); iax++){
    const PandoraPlus::CaloHalfCluster* m_axis = m_axisCol[iax];
    for(int jax=iax+1; jax<m_axisCol.size(); jax++){
      const PandoraPlus::CaloHalfCluster* p_axis = m_axisCol[jax];

      // Do not merge two different track axes
      if ( m_axis->getAssociatedTracks().size()>0 && p_axis->getAssociatedTracks().size()>0 ) continue;
      
      std::vector<const Calo1DCluster*> tmp_localMax = p_axis->getCluster();

      int nsharedHits = 0;
      for(int ihit=0; ihit<m_axis->getCluster().size(); ihit++)
        if( find(tmp_localMax.begin(), tmp_localMax.end(), m_axis->getCluster()[ihit])!=tmp_localMax.end() ) nsharedHits++;

  //printf("  In pair (%d, %d): hit size (%d, %d), shared hit size %d \n",iax, jax, m_axis->getCluster().size(), p_axis->getCluster().size(), nsharedHits);

      if( (m_axis->getCluster().size()<=p_axis->getCluster().size() && (float)nsharedHits/m_axis->getCluster().size()>settings.map_floatPars["th_overlap"] ) ||
          (p_axis->getCluster().size()<m_axis->getCluster().size() && (float)nsharedHits/p_axis->getCluster().size()>settings.map_floatPars["th_overlap"] ) ){

  //cout<<"  Merge: Yes. "<<endl;
        // cout << "    yyy: m_axisCol[" << iax << "] with type=" << m_axisCol[iax]->getType() 
        //      << " and [" << jax << "] with type=" << m_axisCol[jax]->getType() 
        //      << " share " << nsharedHits << " hits, merge them" << endl;
        // cout << "         hit in m_axisCol["<<iax<<"]:"<<endl;
        // for(int yii=0; yii<m_axisCol[iax]->getCluster().size(); yii++){
        //   cout << "           (" << std::fixed << std::setprecision(2) << m_axisCol[iax]->getCluster()[yii]->getPos().x()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[iax]->getCluster()[yii]->getPos().y()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[iax]->getCluster()[yii]->getPos().z()
        //        << ")" << endl;
        // }
        // cout << "         hit in m_axisCol["<<jax<<"]:"<<endl;
        // for(int yii=0; yii<m_axisCol[jax]->getCluster().size(); yii++){
        //   cout << "           (" << std::fixed << std::setprecision(2) << m_axisCol[jax]->getCluster()[yii]->getPos().x()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[jax]->getCluster()[yii]->getPos().y()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[jax]->getCluster()[yii]->getPos().z()
        //        << ")" << endl;
        // }

        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );  

        int axis_type = m_axisCol[iax]->getType() + m_axisCol[jax]->getType();
        m_axisCol[iax]->setType(axis_type);

        // cout << "    yyy: after merge them, new type=" << m_axisCol[iax]->getType() << endl;

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


StatusCode AxisMergingAlg::BranchMerging( std::vector<PandoraPlus::CaloHalfCluster*>& m_axisCol ){
  // cout << "  yyy: calling BranchMerging(), m_axisCol.size()=" << m_axisCol.size() << endl;
  if(m_axisCol.size()<2) return StatusCode::SUCCESS;

  // Merge fake Hough axis to track axis
  std::sort( m_axisCol.begin(), m_axisCol.end(), compLayer );
  for(int iax=0; iax<m_axisCol.size(); iax++){  // select track axis
    const PandoraPlus::CaloHalfCluster* m_axis = m_axisCol[iax];

    if (m_axis->getType()<10000){
      // cout << "yyy:  m_axisCol[" << iax << "] is not a track axis. skip" << endl;
      continue;
    }

    for(int jax=0; jax<m_axisCol.size(); jax++){ // select Hough axis
      if(jax==iax) continue;
      const PandoraPlus::CaloHalfCluster* p_axis = m_axisCol[jax];
      
      if (p_axis->getType()>=10000) continue; // Do not merge two track axis
      if (p_axis->getType()/100%100<1){
        // cout << "yyy:  m_axisCol[" << jax << "] is not a Hough axis. skip" << endl;
        continue;
      }

      // // Determine if two axes are close to each other
      // bool is_close = false;
      // for(int ibar=0; ibar<m_axis->getCluster().size(); ibar++){
      //   for(int jbar=0; jbar<p_axis->getCluster().size(); jbar++){
      //     double distance = ( m_axis->getCluster()[ibar]->getPos() - p_axis->getCluster()[jbar]->getPos() ).Mag();
      //     if(distance<100){ // yyy: hard coding here. min distance of the two axes must be < 50 mm
      //       is_close = true;
      //       break;
      //     }
      //   }
      //   if(is_close) break;
      // }


      // if(!is_close) continue;  // if the two axes are not close to each other, no need to merge

      double hough_rho = p_axis->getHoughRho();
      double hough_alpha = p_axis->getHoughAlpha();

      // V plane
      if(m_axis->getSlayer()==1){
        double x0 = m_axis->getEnergyCenter().x();
        double y0 = m_axis->getEnergyCenter().y();
        double distance = TMath::Abs( x0*TMath::Cos(hough_alpha) + y0*TMath::Sin(hough_alpha) - hough_rho );

        // cout << "    yyy:V rho = " << hough_rho << ", alpha = " << hough_alpha << ", x0 = " << x0 << ", y0 = " << y0 << endl;
        // cout << "         distanceV = " << distance << endl;

        if (distance<settings.map_floatPars["th_branch_distance"]){
          m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );
          int axis_type = m_axisCol[iax]->getType() + m_axisCol[jax]->getType();
          m_axisCol[iax]->setType(axis_type);
          m_axisCol.erase(m_axisCol.begin()+jax);
          jax--;
          if(iax>jax+1) iax--;
          // cout << "  yyy: axis " << jax << " is merged into axis " << iax << endl;
        }
      }

      // U plane
      else{
        int m_module = m_axis->getEnergyCenterTower()[0];
        // cout << "  yyy:branceU: m_module = " << m_module << endl;
        int p_module = p_axis->getTowerID()[0][0];
        // cout << "  yyy:branceU: p_module = " << p_module << endl;
        // Do not merge the two axis in two different modules
        if (m_module != p_module) continue;

        TVector3 t_pos = m_axis->getEnergyCenter();
        // cout << "  yyy:branceU: t_pos = " << t_pos.x() << ", " << t_pos.y() << ", " << t_pos.z() << endl;
        t_pos.RotateZ( TMath::Pi()/4.*(6-m_module) );
        // cout << "  yyy:branceU: t_pos after rotate to module 6 = " << t_pos.x() << ", " << t_pos.y() << ", " << t_pos.z() << endl;
        double x0 = t_pos.x();
        double y0 = t_pos.z();
        double distance = TMath::Abs( x0*TMath::Cos(hough_alpha) + y0*TMath::Sin(hough_alpha) - hough_rho );

        // cout << "    yyy:U rho = " << hough_rho << ", alpha = " << hough_alpha << ", x0 = " << x0 << ", z0 = " << y0 << endl;
        // cout << "         distanceU = " << distance << endl;

        if (distance<settings.map_floatPars["th_branch_distance"]){
          m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );
          int axis_type = m_axisCol[iax]->getType() + m_axisCol[jax]->getType();
          m_axisCol[iax]->setType(axis_type);
          m_axisCol.erase(m_axisCol.begin()+jax);
          jax--;
          if(iax>jax+1) iax--;
          // cout << "  yyy: axis " << jax << " is merged into axis " << iax << endl;
        }


      }

      p_axis=nullptr;
    }
    m_axis=nullptr;
  }
}


StatusCode AxisMergingAlg::FragmentsMerging( std::vector<PandoraPlus::CaloHalfCluster*>& m_axisCol ){
  // cout << "  yyy: calling FragmentsMerging(), m_axisCol.size()=" << m_axisCol.size() << endl;
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
    if( !( (m_axis->getType()==100 && m_axis->getCluster().size()>=settings.map_intPars["th_CoreNhit"] ) ||
        m_axis->getType()>100 ) ) 
      { m_axis=nullptr; continue; }

    // cout << "    yyy: m_axisCol[" << iax << "] is a core" << endl;

    for(int jax=0; jax<m_axisCol.size(); jax++){
      if(jax==iax) continue;
      const PandoraPlus::CaloHalfCluster* p_axis = m_axisCol[jax];
      //Do not merge 2 cores. 
      if( (p_axis->getType()==100 && p_axis->getCluster().size()>=settings.map_intPars["th_CoreNhit"] ) ||
          p_axis->getType()>100 ){
            // cout << "    yyy: m_axisCol[" << jax << "] is also a core, do not merge 2 cores" << endl;
            // cout << "         hit in m_axisCol["<<iax<<"]:"<<endl;
            //for(int yii=0; yii<m_axisCol[iax]->getCluster().size(); yii++){
            //  cout << "           (" << std::fixed << std::setprecision(2) << m_axisCol[iax]->getCluster()[yii]->getPos().x()
            //      << ", " << std::fixed << std::setprecision(2) << m_axisCol[iax]->getCluster()[yii]->getPos().y()
            //      << ", " << std::fixed << std::setprecision(2) << m_axisCol[iax]->getCluster()[yii]->getPos().z()
            //      << ")" << endl;
            //}
            // cout << "         hit in m_axisCol["<<jax<<"]:"<<endl;
            //for(int yii=0; yii<m_axisCol[jax]->getCluster().size(); yii++){
            //  cout << "           (" << std::fixed << std::setprecision(2) << m_axisCol[jax]->getCluster()[yii]->getPos().x()
            //      << ", " << std::fixed << std::setprecision(2) << m_axisCol[jax]->getCluster()[yii]->getPos().y()
            //      << ", " << std::fixed << std::setprecision(2) << m_axisCol[jax]->getCluster()[yii]->getPos().z()
            //      << ")" << endl;
            //}
            
            p_axis=nullptr; 
            continue; 
          }

      TVector3 relP1 = p_axis->getClusterInLayer( p_axis->getBeginningDlayer() )[0]->getPos() - m_axis->getClusterInLayer( m_axis->getEndDlayer() )[0]->getPos();
      TVector3 relP2 = p_axis->getClusterInLayer( p_axis->getBeginningDlayer() )[0]->getPos() - m_axis->getEnergyCenter();
      TVector3 relP3 = m_axis->getEnergyCenter() - p_axis->getClusterInLayer( p_axis->getEndDlayer() )[0]->getPos();
      TVector3 relP4 = p_axis->getPos() - m_axis->getPos();

      // cout << "    yyy: m_axisCol[" << jax << "] is NOT a core, comparing them" << endl;
      // cout << "    yyy: relP1 = (" << std::fixed << std::setprecision(2) << relP1.x() 
      //      <<                 ", " << std::fixed << std::setprecision(2) << relP1.y() 
      //      <<                 ", " << std::fixed << std::setprecision(2) << relP1.z()
      //      << endl
      //      << "         relP2 = (" << std::fixed << std::setprecision(2) << relP2.x() 
      //      <<                 ", " << std::fixed << std::setprecision(2) << relP2.y() 
      //      <<                 ", " << std::fixed << std::setprecision(2) << relP2.z() 
      //      << endl
      //      << "         relP3 = (" << std::fixed << std::setprecision(2) << relP3.x() 
      //      <<                 ", " << std::fixed << std::setprecision(2) << relP3.y() 
      //      <<                 ", " << std::fixed << std::setprecision(2) << relP3.z() 
      //      << endl
      //      << "         relP4 = (" << std::fixed << std::setprecision(2) << relP4.x() 
      //      <<                 ", " << std::fixed << std::setprecision(2) << relP4.y() 
      //      <<                 ", " << std::fixed << std::setprecision(2) << relP4.z() 
      //      << endl
      //      << "         p_axis->getAxis() = (" << std::fixed << std::setprecision(2) << p_axis->getAxis().x() 
      //      <<                             ", " << std::fixed << std::setprecision(2) << p_axis->getAxis().y() 
      //      <<                             ", " << std::fixed << std::setprecision(2) << p_axis->getAxis().z() 
      //      << endl
      //      << "         m_axis->getAxis() = (" << std::fixed << std::setprecision(2) << m_axis->getAxis().x() 
      //      <<                             ", " << std::fixed << std::setprecision(2) << m_axis->getAxis().y() 
      //      <<                             ", " << std::fixed << std::setprecision(2) << m_axis->getAxis().z() 
      //      << endl;

      if( (relP1.Angle(p_axis->getAxis()) < settings.map_floatPars["axis_Angle"] && relP1.Mag()<=settings.map_floatPars["relP_Dis"]) ||
          (relP2.Angle(p_axis->getAxis()) < settings.map_floatPars["axis_Angle"] && relP2.Mag()<=settings.map_floatPars["relP_Dis"]) ||
          (relP3.Angle(p_axis->getAxis()) < settings.map_floatPars["axis_Angle"] && relP3.Mag()<=settings.map_floatPars["relP_Dis"]) || 
          ( p_axis->getType()<100 && relP4.Angle(m_axis->getAxis())<settings.map_floatPars["axis_Angle"] && relP4.Mag()<=settings.map_floatPars["relP_Dis"]) )
      {
        // cout << "    yyy: m_axisCol[" << iax << "] (core) with type=" << m_axisCol[iax]->getType() 
        //      << " and [" << jax << "] (branch) with type=" << m_axisCol[jax]->getType() 
        //      << " are merged" << endl;
        // cout << "         hit in m_axisCol["<<iax<<"]:"<<endl;
        // for(int yii=0; yii<m_axisCol[iax]->getCluster().size(); yii++){
        //   cout << "           (" << std::fixed << std::setprecision(2) << m_axisCol[iax]->getCluster()[yii]->getPos().x()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[iax]->getCluster()[yii]->getPos().y()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[iax]->getCluster()[yii]->getPos().z()
        //        << ")" << endl;
        // }
        // cout << "         hit in m_axisCol["<<jax<<"]:"<<endl;
        // for(int yii=0; yii<m_axisCol[jax]->getCluster().size(); yii++){
        //   cout << "           (" << std::fixed << std::setprecision(2) << m_axisCol[jax]->getCluster()[yii]->getPos().x()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[jax]->getCluster()[yii]->getPos().y()
        //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[jax]->getCluster()[yii]->getPos().z()
        //        << ")" << endl;
        // }
 
        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax] );

        int axis_type = m_axisCol[iax]->getType() + m_axisCol[jax]->getType();
        m_axisCol[iax]->setType(axis_type);

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
  // cout << "  yyy: calling MergeToClosestCluster()" << endl;
  if(m_axisCol.size()==0){ 
    m_axisCol.push_back(m_badaxis);
    // cout << "  yyy: no good axis, treat the bad axis as good axis" << endl;
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
  if(index>=0){
    m_axisCol[index]->mergeHalfCluster(m_badaxis);

    int axis_type = m_axisCol[index]->getType() + m_badaxis->getType();
    m_axisCol[index]->setType(axis_type);

    // cout << "  yyy: merge bad axis to good axis " << index << endl;
    // cout << "       hit in bad axis: " << endl;
    // for(int yii=0; yii<m_badaxis->getCluster().size(); yii++){
    //   cout << "           (" << std::fixed << std::setprecision(2) << m_badaxis->getCluster()[yii]->getPos().x()
    //        << ", " << std::fixed << std::setprecision(2) << m_badaxis->getCluster()[yii]->getPos().y()
    //        << ", " << std::fixed << std::setprecision(2) << m_badaxis->getCluster()[yii]->getPos().z()
    //        << ")" << endl;
    // }
    // cout << "       hit in good axis " << index << ":" << endl;
    // for(int yii=0; yii<m_axisCol[index]->getCluster().size(); yii++){
    //   cout << "           (" << std::fixed << std::setprecision(2) << m_axisCol[index]->getCluster()[yii]->getPos().x()
    //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[index]->getCluster()[yii]->getPos().y()
    //        << ", " << std::fixed << std::setprecision(2) << m_axisCol[index]->getCluster()[yii]->getPos().z()
    //        << ")" << endl;
    // }
  } 

  //delete m_badaxis; m_badaxis=nullptr;
  return StatusCode::SUCCESS;
};

#endif
