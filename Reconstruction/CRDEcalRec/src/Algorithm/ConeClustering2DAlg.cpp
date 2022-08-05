#ifndef _CONECLUSTERING2D_ALG_C
#define _CONECLUSTERING2D_ALG_C

#include "Algorithm/ConeClustering2DAlg.h"

StatusCode ConeClustering2DAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_floatPars.find("th_beginLayer")==settings.map_floatPars.end())   settings.map_floatPars["th_beginLayer"] = 1;
  if(settings.map_floatPars.find("th_stopLayer")==settings.map_floatPars.end())    settings.map_floatPars["th_stopLayer"] = 15;
  if(settings.map_floatPars.find("th_ConeTheta")==settings.map_floatPars.end())    settings.map_floatPars["th_ConeTheta"] = TMath::Pi()/2.;
  if(settings.map_floatPars.find("th_ConeR")==settings.map_floatPars.end())        settings.map_floatPars["th_ConeR"] = 60;
  if(settings.map_floatPars.find("th_Nshowers")==settings.map_floatPars.end())     settings.map_floatPars["th_Nshowers"] = 4;
  if(settings.map_stringPars.find("ReadinLocalMaxName")==settings.map_stringPars.end()) settings.map_stringPars["ReadinLocalMaxName"] = "LeftLocalMax";
  if(settings.map_stringPars.find("OutputLongiClusName")==settings.map_stringPars.end()) settings.map_stringPars["OutputLongiClusName"] = "ConeLongiCluster";

  return StatusCode::SUCCESS;
};


StatusCode ConeClustering2DAlg::Initialize(){

  return StatusCode::SUCCESS;
}

StatusCode ConeClustering2DAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  std::vector<PandoraPlus::Calo3DCluster*>* p_3DClusters = &(m_datacol.Cluster3DCol);

  for(int ic=0; ic<p_3DClusters->size(); ic++){
    std::vector<const PandoraPlus::CaloBarShower*> m_localMaxUCol = p_3DClusters->at(ic)->getLocalMaxUCol(settings.map_stringPars["ReadinLocalMaxName"]);
    std::vector<const PandoraPlus::CaloBarShower*> m_localMaxVCol = p_3DClusters->at(ic)->getLocalMaxVCol(settings.map_stringPars["ReadinLocalMaxName"]);

    if(m_localMaxUCol.size()==0 && m_localMaxVCol.size()==0) continue;

//cout<<"ConeClustering2DAlg: Print input localMax"<<endl;
//cout<<"LocalMaxU: "<<endl;
//for(int i=0; i<m_localMaxUCol.size(); i++)
//  printf("  #%d: (%.3f, %.3f, %.3f) \n ", i, m_localMaxUCol[i]->getPos().x(), m_localMaxUCol[i]->getPos().y(), m_localMaxUCol[i]->getPos().z());
//cout<<"LocalMaxV: "<<endl;
//for(int i=0; i<m_localMaxVCol.size(); i++)
//  printf("  #%d: (%.3f, %.3f, %.3f) \n ", i, m_localMaxVCol[i]->getPos().x(), m_localMaxVCol[i]->getPos().y(), m_localMaxVCol[i]->getPos().z());


    std::map<int, std::vector<const PandoraPlus::CaloBarShower*> > m_orderedShowerU; m_orderedShowerU.clear(); 
    std::map<int, std::vector<const PandoraPlus::CaloBarShower*> > m_orderedShowerV; m_orderedShowerV.clear(); 
    for(int is=0; is<m_localMaxUCol.size(); is++)
      m_orderedShowerU[m_localMaxUCol[is]->getDlayer()].push_back(m_localMaxUCol[is]);
    for(int is=0; is<m_localMaxVCol.size(); is++)
      m_orderedShowerV[m_localMaxVCol[is]->getDlayer()].push_back(m_localMaxVCol[is]);
   

    std::vector<PandoraPlus::LongiCluster*> m_longiClusUCol; m_longiClusUCol.clear(); 
    std::vector<PandoraPlus::LongiCluster*> m_longiClusVCol; m_longiClusVCol.clear(); 

    LongiConeLinking( m_orderedShowerU, m_longiClusUCol );
    LongiConeLinking( m_orderedShowerV, m_longiClusVCol );
   
    //Convert LongiClusters to const object.
    std::vector<const PandoraPlus::LongiCluster*> const_longiClusUCol; const_longiClusUCol.clear();
    std::vector<const PandoraPlus::LongiCluster*> const_longiClusVCol; const_longiClusVCol.clear();
    for(int ic=0; ic<m_longiClusUCol.size(); ic++)
      if(m_longiClusUCol[ic]->getBarShowers().size()>=settings.map_floatPars["th_Nshowers"]) const_longiClusUCol.push_back(m_longiClusUCol[ic]);
     
    for(int ic=0; ic<m_longiClusVCol.size(); ic++) 
      if(m_longiClusVCol[ic]->getBarShowers().size()>=settings.map_floatPars["th_Nshowers"]) const_longiClusVCol.push_back(m_longiClusVCol[ic]);

    for(int ic=0; ic<m_longiClusUCol.size(); ic++) m_datacol.bk_LongiClusCol.push_back( m_longiClusUCol[ic] );
    for(int ic=0; ic<m_longiClusVCol.size(); ic++) m_datacol.bk_LongiClusCol.push_back( m_longiClusVCol[ic] );

    p_3DClusters->at(ic)->setLongiClusters( settings.map_stringPars["OutputLongiClusName"], const_longiClusUCol,
                                            settings.map_stringPars["OutputLongiClusName"], const_longiClusVCol);
  }
  p_3DClusters = nullptr;
  return StatusCode::SUCCESS;
}

StatusCode ConeClustering2DAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}


StatusCode ConeClustering2DAlg::LongiConeLinking(  std::map<int, std::vector<const PandoraPlus::CaloBarShower*> >& orderedShower,
                                                   std::vector<PandoraPlus::LongiCluster*>& ClusterCol)
{
  if(orderedShower.size()==0) return StatusCode::SUCCESS;

  std::map<int, std::vector<const PandoraPlus::CaloBarShower*>>::iterator iter = orderedShower.begin();
  //In first layer: initial clusters. All showers in the first layer are regarded as cluster seed.
  //cluster initial direction = R.
  std::vector<const PandoraPlus::CaloBarShower*> ShowersinFirstLayer;  ShowersinFirstLayer.clear();
  ShowersinFirstLayer = iter->second;
  for(int i=0;i<ShowersinFirstLayer.size(); i++){
    if(iter->first < settings.map_floatPars["th_beginLayer"] || iter->first > settings.map_floatPars["th_stopLayer"] ) continue; 
    PandoraPlus::LongiCluster* m_clus = new  PandoraPlus::LongiCluster();
    m_clus->addBarShower(ShowersinFirstLayer[i]);
    ClusterCol.push_back(m_clus);
  }
  iter++;


  for(iter; iter!=orderedShower.end(); iter++){
    if(iter->first < settings.map_floatPars["th_beginLayer"] || iter->first > settings.map_floatPars["th_stopLayer"] ) continue; 
    std::vector<const PandoraPlus::CaloBarShower*> ShowersinLayer = iter->second;

    for(int is=0; is<ShowersinLayer.size(); is++){
      for(int ic=0; ic<ClusterCol.size(); ic++ ){
        const PandoraPlus::CaloBarShower* shower_in_clus = ClusterCol[ic]->getBarShowers().back();
        if(!shower_in_clus) continue; 

        TVector3 relR = ShowersinLayer[is]->getPos() - shower_in_clus->getPos();
        if( relR.Angle(ClusterCol[ic]->getAxis())<settings.map_floatPars["th_ConeTheta"] && relR.Mag()<settings.map_floatPars["th_ConeR"] ){
          ClusterCol[ic]->addBarShower(ShowersinLayer[is]);
          ShowersinLayer.erase(ShowersinLayer.begin()+is);
          is--;
          break;  
        } 
      }
    }//end loop showers in layer.
    if(ShowersinLayer.size()>0){
      for(int i=0;i<ShowersinLayer.size(); i++){
        PandoraPlus::LongiCluster* m_clus = new PandoraPlus::LongiCluster();
        m_clus->addBarShower(ShowersinLayer[i]);
        ClusterCol.push_back(m_clus);
    }}//end new cluster
  }//end loop layers.

  return StatusCode::SUCCESS;
}

#endif
