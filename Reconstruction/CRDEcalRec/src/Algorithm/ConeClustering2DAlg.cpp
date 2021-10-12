#ifndef _CONECLUSTERING2D_ALG_C
#define _CONECLUSTERING2D_ALG_C

#include "Algorithm/ConeClustering2DAlg.h"

void ConeClustering2DAlg::Settings::SetInitialValue(){
  th_stopLayer = 15; 
  th_ConeTheta = PI/2.;
  th_ConeR = 50.; 

}


StatusCode ConeClustering2DAlg::Initialize(){

  return StatusCode::SUCCESS;
}

StatusCode ConeClustering2DAlg::RunAlgorithm( ConeClustering2DAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol ){
  settings = m_settings; 

  std::vector<CRDEcalEDM::CRDCaloLayer> m_layerCol = m_datacol.LayerCol; 
  std::map<int, std::vector<CRDEcalEDM::CRDCaloBarShower> > m_orderedShowerX; m_orderedShowerX.clear(); 
  std::map<int, std::vector<CRDEcalEDM::CRDCaloBarShower> > m_orderedShowerY; m_orderedShowerY.clear(); 
  for(int il=0; il<m_layerCol.size(); il++){
      m_orderedShowerX[m_layerCol[il].getDlayer()].insert( m_orderedShowerX[m_layerCol[il].getDlayer()].end(), m_layerCol[il].barShowerXCol.begin(), m_layerCol[il].barShowerXCol.end() );
      m_orderedShowerY[m_layerCol[il].getDlayer()].insert( m_orderedShowerY[m_layerCol[il].getDlayer()].end(), m_layerCol[il].barShowerYCol.begin(), m_layerCol[il].barShowerYCol.end() );
  }

  std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_ClusterColX; m_ClusterColX.clear(); 
  std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_ClusterColY; m_ClusterColY.clear(); 
  LongiConeLinking( m_orderedShowerX, m_ClusterColX );
  LongiConeLinking( m_orderedShowerY, m_ClusterColY );

  m_datacol.LongiClusXCol = m_ClusterColX;
  m_datacol.LongiClusYCol = m_ClusterColY; 

  return StatusCode::SUCCESS;
}

StatusCode ConeClustering2DAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}


StatusCode ConeClustering2DAlg::LongiConeLinking(  
  std::map<int, std::vector<CRDEcalEDM::CRDCaloBarShower> >& orderedShower,
  std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& ClusterCol)
{
  if(orderedShower.size()==0) return StatusCode::SUCCESS;

  std::map<int, std::vector<CRDEcalEDM::CRDCaloBarShower>>::iterator iter = orderedShower.begin();
  //In first layer: initial clusters. All showers in the first layer are regarded as cluster seed.
  //cluster initial direction = R.
  std::vector<CRDEcalEDM::CRDCaloBarShower> ShowersinFirstLayer;  ShowersinFirstLayer.clear();
  ShowersinFirstLayer = iter->second;
  for(int i=0;i<ShowersinFirstLayer.size(); i++){
    if(iter->first > settings.th_stopLayer) continue; 
    CRDEcalEDM::CRDCaloHitLongiCluster m_clus;
    m_clus.AddBarShower(ShowersinFirstLayer[i]);
    ClusterCol.push_back(m_clus);
  }
  iter++;


  for(iter;iter!=orderedShower.end();iter++){
    if(iter->first > settings.th_stopLayer) continue; 
    std::vector<CRDEcalEDM::CRDCaloBarShower> ShowersinLayer = iter->second;

    for(int is=0; is<ShowersinLayer.size(); is++){
      CRDEcalEDM::CRDCaloBarShower m_shower = ShowersinLayer[is];
      for(int ic=0; ic<ClusterCol.size(); ic++ ){
        CRDEcalEDM::CRDCaloBarShower shower_in_clus = ClusterCol[ic].getBarShowers().back();
        dd4hep::Position relR = m_shower.getPos()-shower_in_clus.getPos();

        TVector3 relR_vec(relR.x(), relR.y(), relR.z());
        if( relR_vec.Angle(ClusterCol[ic].getAxis())<settings.th_ConeTheta && relR_vec.Mag()<settings.th_ConeR ){
          ClusterCol[ic].AddBarShower(m_shower);
          ShowersinLayer.erase(ShowersinLayer.begin()+is);
          is--;
          break;  
        } 
      }
    }//end loop showers in layer.
    if(ShowersinLayer.size()>0){
      for(int i=0;i<ShowersinLayer.size(); i++){
        CRDEcalEDM::CRDCaloHitLongiCluster m_clus;
        m_clus.AddBarShower(ShowersinLayer[i]);
        ClusterCol.push_back(m_clus);
    }}//end new cluster
  }//end loop layers.

  return StatusCode::SUCCESS;
}

#endif
