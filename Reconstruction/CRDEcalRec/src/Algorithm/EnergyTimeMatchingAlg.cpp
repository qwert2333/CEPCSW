#ifndef _ENERGYTIMEMATCHING_ALG_C
#define _ENERGYTIMEMATCHING_ALG_C

#include "Algorithm/EnergyTimeMatchingAlg.h"

void EnergyTimeMatchingAlg::Settings::SetInitialValue(){
  chi2Wi_E = 1;
  chi2Wi_T = 10;
  th_chi2  = -1;
  sigmaE   = 0.10;
  fl_UseChi2 = true;
  Debug = 0;
}


EnergyTimeMatchingAlg::EnergyTimeMatchingAlg(){

}

StatusCode Initialize(){
	return StatusCode::SUCCESS;
}

StatusCode EnergyTimeMatchingAlg::RunAlgorithm( EnergyTimeMatchingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol ){
  settings = m_settings;

  std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_2DshowerCol; m_2DshowerCol.clear();
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_ClusterCol;  m_ClusterCol.clear();  

  //Basic unit: tower. 
  std::vector<CRDEcalEDM::CRDCaloTower> m_towerCol = m_datacol.TowerCol;
  if(m_towerCol.size()==0){ std::cout<<"Warning: Empty input in EnergyTimeMatchingAlg! Please check previous algorithm!"<<std::endl; return StatusCode::SUCCESS; }

  for(int it=0; it<m_towerCol.size(); it++){
    std::vector<CRDEcalEDM::CRDCaloBlock> m_blocks = m_towerCol[it].getBlocks(); 

    //Case1: Don't use shadow cluster or chi2 matching, Save all combinations for 3D clustering.
    if(!settings.fl_UseChi2){
      for(int ib=0;ib<m_blocks.size();ib++){
        std::vector<CRDEcalEDM::CRDCaloBarShower> showerXCol = m_blocks[ib].barShowerXCol;
        std::vector<CRDEcalEDM::CRDCaloBarShower> showerYCol = m_blocks[ib].barShowerYCol;
   
        std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_showerinlayer; m_showerinlayer.clear();
        GetFullMatchedShowers( showerXCol, showerYCol, m_showerinlayer );
   
        m_2DshowerCol.insert( m_2DshowerCol.end(), m_showerinlayer.begin(), m_showerinlayer.end() );
   
      }
      m_datacol.Shower2DCol = m_2DshowerCol;
    }

    //Case2: chi2 matching (combined chi2 in all layers)
    else{
      std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_longiClXCol = m_towerCol[it].getLongiClusterXCol();
      std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_longiClYCol = m_towerCol[it].getLongiClusterYCol();

      const int NclusX = m_longiClXCol.size(); 
      const int NclusY = m_longiClYCol.size(); 
      if(NclusX==0 || NclusY==0) continue; 

      double chi2[14][NclusX][NclusY] = {0}; 
      double sumchi2[NclusX][NclusY] = {0};

      std::vector<CRDEcalEDM::CRDCaloHit3DCluster> tmp_clusters; tmp_clusters.clear(); 
      CRDEcalEDM::CRDCaloHit3DCluster tmp_clus; tmp_clus.Clear(); 
      
      //Case 2.1
      if(NclusX==1 && NclusY==1){
        //chi2=0; sumchi2=0; 
        XYShowerMatchingL0(m_longiClXCol[0], m_longiClYCol[0], tmp_clus); 
        tmp_clusters.push_back(m_3Dclus);
        continue;
      }
      //Case 2.2
      else if(NclusX==1){ XYShowerMatchingL1(m_longiClXCol[0], m_longiClYCol, tmp_clusters); continue; }
      else if(NclusY==1){ XYShowerMatchingL1(m_longiClYCol[0], m_longiClXCol, tmp_clusters); continue; }

      //Case 2.3
      else if( NclusX==NclusY ){ 
        XYShowerChi2Matching(m_longiClXCol, m_longiClYCol, chi2);


      }
      //Case 2.4
      else{ 
        XYShowerChi2MatchingL1(m_longiClXCol, m_longiClYCol, chi2);

      
      }

    }


  }



  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetFullMatchedShowers( 
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, 
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol, 
           std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol )
{
  outshCol.clear(); 

  for(int is=0;is<barShowerXCol.size();is++){
  for(int js=0;js<barShowerYCol.size();js++){
    if(barShowerXCol[is].getBars().size()==0 || barShowerYCol[js].getBars().size()==0) continue;
    CRDEcalEDM::CRDCaloHitTransShower tmp_shower; tmp_shower.Clear();
    XYShowerMatchingL0(barShowerXCol[is], barShowerYCol[js], tmp_shower);
    outshCol.push_back(tmp_shower);
  }}

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::XYShowerMatchingL0( CRDEcalEDM::CRDCaloHitLongiCluster& m_longiClX, 
                                                      CRDEcalEDM::CRDCaloHitLongiCluster& m_longiClY, 
                                                      CRDEcalEDM::CRDCaloHit3DCluster& m_clus )
{
  m_clus.Clear(); 

  std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_showers;  m_showers.clear(); 

  std::vector<int> layerindex; layerindex.clear(); 
  for(int is=0; is<m_longiClX.getBarShowers(); is++){
    if( find( layerindex.begin(), layerindex.end(), m_longiClX.getBarShowers()[is].getDlayer() )==layerindex.end() ) layerindex.push_back(m_longiClX.getBarShowers()[is].getDlayer());
  }
  for(int is=0; is<m_longiClY.getBarShowers(); is++){
    if( find( layerindex.begin(), layerindex.end(), m_longiClY.getBarShowers()[is].getDlayer() )==layerindex.end() ) layerindex.push_back(m_longiClY.getBarShowers()[is].getDlayer());
  }

  
  for(int il=0; il<layerindex.size(); il++){
    




  }


  return StatusCode::SUCCESS;
}


StatusCode XYShowerMatchingL1( CRDEcalEDM::CRDCaloHitLongiCluster& m_longiCl1, 
                               std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClN, 
                               std::vector<CRDEcalEDM::CRDCaloHit3DCluster>& m_clusters )
{


  return StatusCode::SUCCESS;
}

StatusCode XYShowerChi2Matching( std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClXCol, 
                                 std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClYCol, 
                                 double** chi2 )
{


  return StatusCode::SUCCESS;
}
StatusCode XYShowerChi2MatchingL1( std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClXCol, 
                                   std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClYCol, 
                                   double** chi2 )
{



  return StatusCode::SUCCESS;
}


#endif
