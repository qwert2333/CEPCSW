#ifndef _CONECLUSTERING2D_ALG_C
#define _CONECLUSTERING2D_ALG_C

#include "Algorithm/ConeClustering2DAlg.h"

StatusCode ConeClustering2DAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_floatPars.find("th_beginLayer")==settings.map_floatPars.end())   settings.map_floatPars["th_beginLayer"] = 1;
  if(settings.map_floatPars.find("th_stopLayer")==settings.map_floatPars.end())    settings.map_floatPars["th_stopLayer"] = 15;
  if(settings.map_floatPars.find("th_ConeTheta")==settings.map_floatPars.end())    settings.map_floatPars["th_ConeTheta"] = TMath::Pi()/4.;
  if(settings.map_floatPars.find("th_ConeR")==settings.map_floatPars.end())        settings.map_floatPars["th_ConeR"] = 50;
  if(settings.map_intPars.find("th_Nshowers")==settings.map_intPars.end())       settings.map_intPars["th_Nshowers"] = 4;
  if(settings.map_stringPars.find("ReadinLocalMaxName")==settings.map_stringPars.end()) settings.map_stringPars["ReadinLocalMaxName"] = "LeftLocalMax";
  if(settings.map_stringPars.find("OutputLongiClusName")==settings.map_stringPars.end()) settings.map_stringPars["OutputLongiClusName"] = "ConeAxis";

  return StatusCode::SUCCESS;
};


StatusCode ConeClustering2DAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  p_HalfClusterU.clear();
  p_HalfClusterV.clear();


  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColU"].size(); ih++)
    p_HalfClusterU.push_back( m_datacol.map_HalfCluster["HalfClusterColU"][ih].get() );
  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColV"].size(); ih++)
    p_HalfClusterV.push_back( m_datacol.map_HalfCluster["HalfClusterColV"][ih].get() );

  return StatusCode::SUCCESS;
}

StatusCode ConeClustering2DAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  if( (p_HalfClusterU.size()+p_HalfClusterV.size())<1 ){
    std::cout << "ConeClustering2DAlg: No HalfCluster input"<<std::endl;
    return StatusCode::SUCCESS;
  }   
 
  std::vector<PandoraPlus::CaloHalfCluster*> tmp_longiClusCol; tmp_longiClusCol.clear(); 
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > m_orderedLocalMax; 

//cout<<"  ConeClustering: Input HalfCluster size "<<p_HalfClusterU.size()<<", "<<p_HalfClusterV.size()<<endl;

  //Processing U plane:
  for(int ic=0; ic<p_HalfClusterU.size(); ic++){
    //Get LocalMax
    m_localMaxUCol.clear(); 
    m_localMaxUCol = p_HalfClusterU[ic]->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

//cout<<"    In ClusU #"<<ic<<": localMax size "<<m_localMaxUCol.size()<<endl;

    //Get ordered LocalMax
    m_orderedLocalMax.clear();
    for(int is=0; is<m_localMaxUCol.size(); is++)
      m_orderedLocalMax[m_localMaxUCol[is]->getDlayer()].push_back(m_localMaxUCol[is]);

    tmp_longiClusCol.clear();
    LongiConeLinking( m_orderedLocalMax, tmp_longiClusCol, m_datacol.map_HalfCluster["bkHalfCluster"] );

//cout<<"    Cluster size after ConeLinking: "<<tmp_longiClusCol.size()<<endl;

    //Convert LongiClusters to const object.    
    const_longiClusUCol.clear();
    for(int ic=0; ic<tmp_longiClusCol.size(); ic++){
      tmp_longiClusCol[ic]->Check();
      if(tmp_longiClusCol[ic]->getCluster().size()>=settings.map_intPars["th_Nshowers"]) const_longiClusUCol.push_back(tmp_longiClusCol[ic]);
    }
//cout<<"    Cluster size after requiring Nhit: "<<const_longiClusUCol.size()<<endl;

    p_HalfClusterU[ic]->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], const_longiClusUCol);    
  }

  //Processing V plane:
  for(int ic=0; ic<p_HalfClusterV.size(); ic++){
    m_localMaxVCol.clear();
    m_localMaxVCol = p_HalfClusterV[ic]->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

    m_orderedLocalMax.clear();
    for(int is=0; is<m_localMaxVCol.size(); is++)
      m_orderedLocalMax[m_localMaxVCol[is]->getDlayer()].push_back(m_localMaxVCol[is]);

    tmp_longiClusCol.clear();
    LongiConeLinking(m_orderedLocalMax, tmp_longiClusCol, m_datacol.map_HalfCluster["bkHalfCluster"]);


    //Convert LongiClusters to const object.
    const_longiClusVCol.clear();
    for(int ic=0; ic<tmp_longiClusCol.size(); ic++){
      tmp_longiClusCol[ic]->Check();
      if(tmp_longiClusCol[ic]->getCluster().size()>=settings.map_intPars["th_Nshowers"]) const_longiClusVCol.push_back(tmp_longiClusCol[ic]);
    }

    p_HalfClusterV[ic]->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], const_longiClusVCol);
  }

  tmp_longiClusCol.clear();
  m_orderedLocalMax.clear();
  return StatusCode::SUCCESS;
}

StatusCode ConeClustering2DAlg::ClearAlgorithm(){
  p_HalfClusterV.clear(); 
  p_HalfClusterU.clear();
  m_localMaxVCol.clear();
  m_localMaxUCol.clear();
  const_longiClusVCol.clear();
  const_longiClusUCol.clear();

  return StatusCode::SUCCESS;
}


StatusCode ConeClustering2DAlg::LongiConeLinking(  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> >& orderedShower,
                                                   std::vector<PandoraPlus::CaloHalfCluster*>& ClusterCol, 
                                                   std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& bk_HFclus )
{
  if(orderedShower.size()==0) return StatusCode::SUCCESS;
 
  vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_clusCol; m_clusCol.clear();
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*>>::iterator iter = orderedShower.begin();
  //In first layer: initial clusters. All showers in the first layer are regarded as cluster seed.
  //cluster initial direction = R.
  std::vector<const PandoraPlus::Calo1DCluster*> ShowersinFirstLayer;  ShowersinFirstLayer.clear();
  ShowersinFirstLayer = iter->second;
  for(int i=0;i<ShowersinFirstLayer.size(); i++){
    if(iter->first < settings.map_floatPars["th_beginLayer"] || iter->first > settings.map_floatPars["th_stopLayer"] ) continue; 
    std::shared_ptr<PandoraPlus::CaloHalfCluster> m_clus = std::make_shared<PandoraPlus::CaloHalfCluster>();
    m_clus->addUnit(ShowersinFirstLayer[i]);
    m_clus->setType(1);
    m_clusCol.push_back(m_clus);
  }
  iter++;


  for(iter; iter!=orderedShower.end(); iter++){
    if(iter->first < settings.map_floatPars["th_beginLayer"] || iter->first > settings.map_floatPars["th_stopLayer"] ) continue; 
    std::vector<const PandoraPlus::Calo1DCluster*> ShowersinLayer = iter->second;
//cout<<"  In Layer "<<iter->first<<": hit size "<<ShowersinLayer.size()<<endl;

    for(int ic=0; ic<m_clusCol.size(); ic++){
      const PandoraPlus::Calo1DCluster* shower_in_clus = m_clusCol[ic].get()->getCluster().back();
      if(!shower_in_clus) continue;
      for(int is=0; is<ShowersinLayer.size(); is++){
        TVector2 relR = GetProjectedRelR(shower_in_clus, ShowersinLayer[is]);  //Return vec: 1->2.
        TVector2 clusaxis = GetProjectedAxis( m_clusCol[ic].get() );
//printf("    Cluster axis (%.3f, %.3f), last hit (%.3f, %.3f, %.3f), coming hit (%.3f, %.3f, %.3f), projected relR (%.3f, %.3f) \n",
//clusaxis.Px(), clusaxis.Py(), shower_in_clus->getPos().x(),  shower_in_clus->getPos().y(),  shower_in_clus->getPos().z(),
//ShowersinLayer[is]->getPos().x(), ShowersinLayer[is]->getPos().y(), ShowersinLayer[is]->getPos().z(), relR.Px(), relR.Py() );

        if( relR.DeltaPhi(clusaxis)<settings.map_floatPars["th_ConeTheta"] && relR.Mod()<settings.map_floatPars["th_ConeR"] ){
          m_clusCol[ic].get()->addUnit(ShowersinLayer[is]);
          ShowersinLayer.erase(ShowersinLayer.begin()+is);
          is--;
        }
      }
    }//end loop showers in layer.

    if(ShowersinLayer.size()>0){
      for(int i=0;i<ShowersinLayer.size(); i++){
        std::shared_ptr<PandoraPlus::CaloHalfCluster> m_clus = std::make_shared<PandoraPlus::CaloHalfCluster>();
        m_clus->addUnit(ShowersinLayer[i]);
        m_clus->setType(1);
        m_clusCol.push_back(m_clus);
    }}//end new cluster
  }//end loop layers.

  bk_HFclus.insert(bk_HFclus.end(), m_clusCol.begin(), m_clusCol.end());
  for(int icl=0; icl<m_clusCol.size(); icl++)
    ClusterCol.push_back( m_clusCol[icl].get() );  

  return StatusCode::SUCCESS;
}



TVector2 ConeClustering2DAlg::GetProjectedAxis( const PandoraPlus::CaloHalfCluster* m_shower ){
  TVector2 axis(0., 0.);
  TVector3 m_axis = m_shower->getAxis();
  if( m_shower->getSlayer()==1 )  axis.Set(m_axis.x(), m_axis.y()); //For V-bars: (x, y) 
  else                            axis.Set( sqrt(m_axis.x()*m_axis.x() + m_axis.y()*m_axis.y()), m_axis.z()); //For U-bars: (R, z)
  axis.Unit();

  return axis;
}

TVector2 ConeClustering2DAlg::GetProjectedRelR( const PandoraPlus::Calo1DCluster* m_shower1, const PandoraPlus::Calo1DCluster* m_shower2 ){
  TVector2 paxis1, paxis2;
  if(m_shower1->getSlayer()==1){ //For V-bars
    paxis1.Set(m_shower1->getPos().x(), m_shower1->getPos().y());
    paxis2.Set(m_shower2->getPos().x(), m_shower2->getPos().y());
  }
  else{
    float rotAngle = -m_shower1->getTowerID()[0][0]*PI/4.;
    TVector3 raw_pos1 = m_shower1->getPos(); 
    TVector3 raw_pos2 = m_shower2->getPos(); 
    raw_pos1.RotateZ(rotAngle); raw_pos2.RotateZ(rotAngle);
    paxis1.Set(raw_pos1.y(), raw_pos1.z());
    paxis2.Set(raw_pos2.y(), raw_pos2.z());

    //paxis1.Set( sqrt(m_shower1->getPos().x()*m_shower1->getPos().x()+m_shower1->getPos().y()*m_shower1->getPos().y()), m_shower1->getPos().z());
    //paxis2.Set( sqrt(m_shower2->getPos().x()*m_shower2->getPos().x()+m_shower2->getPos().y()*m_shower2->getPos().y()), m_shower2->getPos().z());
  }

  return paxis2 - paxis1;
}

#endif
