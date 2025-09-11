#ifndef ETMATCHING_ALG_C
#define ETMATCHING_ALG_C

#include "Algorithm/EnergyTimeMatchingAlg.h"
StatusCode EnergyTimeMatchingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Set initial value
  if(settings.map_floatPars.find("chi2Wi_E")==settings.map_floatPars.end())          settings.map_floatPars["chi2Wi_E"] = 1.;
  if(settings.map_floatPars.find("chi2Wi_T")==settings.map_floatPars.end())          settings.map_floatPars["chi2Wi_T"] = 10.;
  if(settings.map_floatPars.find("sigmaE")==settings.map_floatPars.end())            settings.map_floatPars["sigmaE"] = 0.10;
  if(settings.map_floatPars.find("sigmaPos")==settings.map_floatPars.end())          settings.map_floatPars["sigmaPos"] = 34.89;
  if(settings.map_floatPars.find("nMat")==settings.map_floatPars.end())              settings.map_floatPars["nMat"] = 2.15;
  if(settings.map_floatPars.find("Eth_HFClus")==settings.map_floatPars.end())        settings.map_floatPars["Eth_HFClus"] = 0.05;
  if(settings.map_floatPars.find("th_overlapE")==settings.map_floatPars.end())       settings.map_floatPars["th_overlapE"] = 0.3;
  if(settings.map_floatPars.find("th_UVdeltaE")==settings.map_floatPars.end())       settings.map_floatPars["th_UVdeltaE"] = 0.3;
  if(settings.map_intPars.find("compressLayer")==settings.map_intPars.end())       settings.map_intPars["compressLayer"] = 3;
  if(settings.map_floatPars.find("th_UVdeltaEinLayer")==settings.map_floatPars.end()) settings.map_floatPars["th_UVdeltaEinLayer"] = 0.5;
  if(settings.map_floatPars.find("th_ConeTheta")==settings.map_floatPars.end())      settings.map_floatPars["th_ConeTheta"] = TMath::Pi()/6.;
  if(settings.map_floatPars.find("th_ConeR")==settings.map_floatPars.end())          settings.map_floatPars["th_ConeR"] = 30.;

  if(settings.map_stringPars.find("ReadinHFClusterName")==settings.map_stringPars.end()) settings.map_stringPars["ReadinHFClusterName"] = "ESHalfCluster";
  if(settings.map_stringPars.find("ReadinTowerName")==settings.map_stringPars.end()) settings.map_stringPars["ReadinTowerName"] = "ESTower";
  if(settings.map_stringPars.find("OutputClusterName")==settings.map_stringPars.end()) settings.map_stringPars["OutputClusterName"] = "EcalCluster";

  return StatusCode::SUCCESS;
};


StatusCode EnergyTimeMatchingAlg::Initialize( CyberDataCol& m_datacol ){
  m_HFClusUCol.clear();
  m_HFClusVCol.clear();
  m_clusterCol.clear();
  m_towerCol.clear();
  m_bkCol.Clear();

  int ntower = m_datacol.map_CaloCluster[settings.map_stringPars["ReadinTowerName"]].size(); 
  for(int it=0; it<ntower; it++)
    m_towerCol.push_back( m_datacol.map_CaloCluster[settings.map_stringPars["ReadinTowerName"]][it].get() );

//cout<<"  EnergyTimeMatchingAlg: Readin tower size: "<<m_towerCol.size()<<endl;
	return StatusCode::SUCCESS;
};


StatusCode EnergyTimeMatchingAlg::RunAlgorithm( CyberDataCol& m_datacol ){


  m_clusterCol.clear();
  //std::vector<const Cyber::CaloHalfCluster*> m_leftHFClusterUCol; 
  //std::vector<const Cyber::CaloHalfCluster*> m_leftHFClusterVCol; 
  //Loop for towers:  
  for(int it=0; it<m_towerCol.size(); it++){
    m_HFClusUCol.clear(); m_HFClusVCol.clear();
    m_emptyHFClusUCol.clear(); m_emptyHFClusVCol.clear();

//cout<<"Check tower ID: ";
//for(int i=0; i<m_towerCol[it]->getTowerID().size(); i++) printf("[%d, %d, %d], ", m_towerCol[it]->getTowerID()[i][0], m_towerCol[it]->getTowerID()[i][1], m_towerCol[it]->getTowerID()[i][2]);
//cout<<endl;

    m_HFClusUCol = m_towerCol.at(it)->getHalfClusterUCol(settings.map_stringPars["ReadinHFClusterName"]+"U");
    m_HFClusVCol = m_towerCol.at(it)->getHalfClusterVCol(settings.map_stringPars["ReadinHFClusterName"]+"V");
    //cout<<"  Axis HFCluster: ["<<m_HFClusUCol.size()<<", "<<m_HFClusVCol.size()<<"] "<<endl;

    double tmp_towerEn = m_towerCol.at(it)->getEnergy();
    if(tmp_towerEn<1.){
      for(int ih=0; ih<m_HFClusUCol.size(); ih++){
        if(m_HFClusUCol[ih]->getEnergy()<0.05*tmp_towerEn){
          m_HFClusUCol.erase(m_HFClusUCol.begin()+ih);
        }
      }
      for(int ih=0; ih<m_HFClusVCol.size(); ih++){
        if(m_HFClusVCol[ih]->getEnergy()<0.05*tmp_towerEn){
          m_HFClusVCol.erase(m_HFClusVCol.begin()+ih);
        }
      }
    }
  
    //m_emptyHFClusUCol = m_towerCol.at(it)->getHalfClusterUCol("emptyHalfClusterU");
    //m_emptyHFClusVCol = m_towerCol.at(it)->getHalfClusterVCol("emptyHalfClusterV");

//cout<<"  Axis HFCluster after filter: ["<<m_HFClusUCol.size()<<", "<<m_HFClusVCol.size()<<"] "<<endl;
//cout<<"  Empty HFCluster: ["<<m_emptyHFClusUCol.size()<<", "<<m_emptyHFClusVCol.size()<<"] "<<endl;
//cout<<"  3D cluster size (step 0): "<<m_clusterCol.size()<<endl;
    if(m_HFClusUCol.size()>0 && m_HFClusVCol.size()>0)
      Matching(m_HFClusUCol, m_HFClusVCol, m_clusterCol);
    //else if(m_HFClusUCol.size()>0 && m_emptyHFClusVCol.size()>0)
    //  Matching(m_HFClusUCol, m_emptyHFClusVCol, m_clusterCol);
    //else if(m_HFClusVCol.size()>0 && m_emptyHFClusUCol.size()>0)
    //  Matching(m_emptyHFClusUCol, m_HFClusVCol, m_clusterCol);
    //else      
    //  Matching(m_emptyHFClusUCol, m_emptyHFClusVCol, m_clusterCol);

  }//End loop towers

  ClusterReconnecting( m_clusterCol );
  
//cout<<"  3D cluster size (after reconnection): "<<m_clusterCol.size()<<endl;
  m_datacol.map_CaloCluster[settings.map_stringPars["OutputClusterName"]] = m_clusterCol;

//cout<<" Save backup collections in to main datacol. "<<endl;
  //Save backup collections in to main datacol. 
  m_datacol.map_CaloHit["bkHit"].insert( m_datacol.map_CaloHit["bkHit"].end(), m_bkCol.map_CaloHit["bkHit"].begin(), m_bkCol.map_CaloHit["bkHit"].end() );
  m_datacol.map_BarCol["bkBar"].insert( m_datacol.map_BarCol["bkBar"].end(), m_bkCol.map_BarCol["bkBar"].begin(), m_bkCol.map_BarCol["bkBar"].end() );
  m_datacol.map_1DCluster["bk1DCluster"].insert( m_datacol.map_1DCluster["bk1DCluster"].end(), m_bkCol.map_1DCluster["bk1DCluster"].begin(), m_bkCol.map_1DCluster["bk1DCluster"].end() );
  m_datacol.map_2DCluster["bk2DCluster"].insert( m_datacol.map_2DCluster["bk2DCluster"].end(), m_bkCol.map_2DCluster["bk2DCluster"].begin(), m_bkCol.map_2DCluster["bk2DCluster"].end() );
  m_datacol.map_HalfCluster["bkHalfCluster"].insert( m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_bkCol.map_HalfCluster["bkHalfCluster"].begin(), m_bkCol.map_HalfCluster["bkHalfCluster"].end() );

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::ClearAlgorithm(){
  m_HFClusUCol.clear();
  m_HFClusVCol.clear();
  m_clusterCol.clear();
  m_towerCol.clear();
  m_bkCol.Clear();
  

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::Matching( std::vector<const Cyber::CaloHalfCluster*> m_HFClusUCol,
                                            std::vector<const Cyber::CaloHalfCluster*> m_HFClusVCol,
                                            std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_clusterCol )
{

  const int NclusU = m_HFClusUCol.size();
  const int NclusV = m_HFClusVCol.size();

cout<<" cluster size ["<<NclusU<<", "<<NclusV<<"] "<<endl;
cout<<"ClusterU energy: ";
for(int i=0; i<NclusU; i++) cout<<m_HFClusUCol[i]->getEnergy()<<'\t';
cout<<endl;
cout<<"ClusterV energy: ";
for(int i=0; i<NclusV; i++) cout<<m_HFClusVCol[i]->getEnergy()<<'\t';
cout<<endl;


  std::vector<std::shared_ptr<Cyber::Calo3DCluster>> tmp_clusters; tmp_clusters.clear();
  if(NclusU==0 || NclusV==0) return StatusCode::SUCCESS;

  else if(NclusU==1 && NclusV==1){
      std::shared_ptr<Cyber::Calo3DCluster> tmp_clus = std::make_shared<Cyber::Calo3DCluster>();
      XYClusterMatchingL0(m_HFClusUCol[0], m_HFClusVCol[0], tmp_clus);
      tmp_clusters.push_back(tmp_clus);
  }  
  else if(NclusU==1){
    double sumE_V = 0.;
    for(int icl=0; icl<NclusV; icl++) sumE_V += m_HFClusVCol[icl]->getEnergy();

    //Split the HalfClusterU into NclusV part
    for(int icl=0; icl<NclusV; icl++){
      std::shared_ptr<Cyber::CaloHalfCluster> newClus = std::make_shared<Cyber::CaloHalfCluster>();
      double frac = m_HFClusVCol[icl]->getEnergy()/sumE_V;
      for(int ish=0; ish<m_HFClusUCol[0]->getCluster().size(); ish++){
        const Cyber::Calo1DCluster* p_shower = m_HFClusUCol[0]->getCluster()[ish];

        std::vector<const Cyber::CaloUnit*> Bars; Bars.clear();
        for(int ibar=0; ibar<p_shower->getCluster().size(); ibar++){
          auto bar = p_shower->getCluster()[ibar]->Clone();
          bar->setQ(bar->getQ1()*frac, bar->getQ2()*frac );
          Bars.push_back(bar.get());
          m_bkCol.map_BarCol["bkBar"].push_back(bar);
        }

        std::shared_ptr<Cyber::Calo1DCluster> shower = std::make_shared<Cyber::Calo1DCluster>();
        shower->setBars(Bars);
        shower->setSeed();
        shower->setIDInfo();

        newClus->addUnit(shower.get());
        m_bkCol.map_1DCluster["bk1DCluster"].push_back( shower );
      }
      newClus->setHoughPars(m_HFClusUCol[0]->getHoughAlpha(), m_HFClusUCol[0]->getHoughRho());
      newClus->setIntercept(m_HFClusUCol[0]->getHoughIntercept());
      for(auto iter: m_HFClusUCol[0]->getHalfClusterMap()) newClus->setHalfClusters( iter.first, iter.second );
      for(int itrk=0; itrk<m_HFClusUCol[0]->getAssociatedTracks().size(); itrk++) newClus->addAssociatedTrack(m_HFClusUCol[0]->getAssociatedTracks()[itrk]);
      newClus->setLinkedMCP(m_HFClusUCol[0]->getLinkedMCP());
      newClus->setType(m_HFClusUCol[0]->getType());
      m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(newClus);

      //Match the splitted HalfClusterU and HalfClusterV[icl]
      std::shared_ptr<Cyber::Calo3DCluster> tmp_clus = std::make_shared<Cyber::Calo3DCluster>();
      XYClusterMatchingL0(newClus.get(), m_HFClusVCol[icl], tmp_clus);
      tmp_clusters.push_back(tmp_clus);
    }
  }

  else if(NclusV==1){
    double sumE_U = 0.;
    for(int icl=0; icl<NclusU; icl++) sumE_U += m_HFClusUCol[icl]->getEnergy();

    //Split the HalfClusterV into NclusU part
    for(int icl=0; icl<NclusU; icl++){
      std::shared_ptr<Cyber::CaloHalfCluster> newClus = std::make_shared<Cyber::CaloHalfCluster>();
      double frac = m_HFClusUCol[icl]->getEnergy()/sumE_U;
      for(int ish=0; ish<m_HFClusVCol[0]->getCluster().size(); ish++){
        const Cyber::Calo1DCluster* p_shower = m_HFClusVCol[0]->getCluster()[ish];

        std::vector<const Cyber::CaloUnit*> Bars; Bars.clear();
        for(int ibar=0; ibar<p_shower->getCluster().size(); ibar++){
          auto bar = p_shower->getCluster()[ibar]->Clone();
          bar->setQ(bar->getQ1()*frac, bar->getQ2()*frac );
          Bars.push_back(bar.get());
          m_bkCol.map_BarCol["bkBar"].push_back(bar);
        }

        std::shared_ptr<Cyber::Calo1DCluster> shower = std::make_shared<Cyber::Calo1DCluster>();
        shower->setBars(Bars);
        shower->setSeed();
        shower->setIDInfo();

        newClus->addUnit(shower.get());
        m_bkCol.map_1DCluster["bk1DCluster"].push_back( shower );
      }
      newClus->setHoughPars(m_HFClusVCol[0]->getHoughAlpha(), m_HFClusVCol[0]->getHoughRho());
      newClus->setIntercept(m_HFClusVCol[0]->getHoughIntercept());
      for(auto iter: m_HFClusVCol[0]->getHalfClusterMap()) newClus->setHalfClusters( iter.first, iter.second );
      for(int itrk=0; itrk<m_HFClusVCol[0]->getAssociatedTracks().size(); itrk++) newClus->addAssociatedTrack(m_HFClusVCol[0]->getAssociatedTracks()[itrk]);
      newClus->setLinkedMCP(m_HFClusVCol[0]->getLinkedMCP());
      newClus->setType(m_HFClusVCol[0]->getType());
      m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(newClus);

      //Match the splitted HalfClusterU and HalfClusterV[icl]
      std::shared_ptr<Cyber::Calo3DCluster> tmp_clus = std::make_shared<Cyber::Calo3DCluster>();
      XYClusterMatchingL0(m_HFClusUCol[icl], newClus.get(), tmp_clus);
      tmp_clusters.push_back(tmp_clus);

    }
  }

  else{
    vector<vector<double>> Chi2Matrix(NclusV, vector<double>(NclusU, 0.));
    std::vector< std::pair<int, int> > Chi2MatrixOrder; Chi2MatrixOrder.clear();
    //Calculate the chi2 matrix with energy
    Chi2MatrixCalculation(m_HFClusUCol, m_HFClusVCol, Chi2Matrix, Chi2MatrixOrder);

//cout<<"Print chi2 matrix "<<endl;
//for (const auto& row : Chi2Matrix) {
//    for (const double& val : row) {
//        std::cout << val << " ";
//    }
//    std::cout << std::endl;
//}

    //Calculate the initial pattern matrix with track + neighbor
    vector<vector<int>> PatternMatrix(NclusV, vector<int>(NclusU, 0));
    PatternMatrixCalculation(m_HFClusUCol, m_HFClusVCol, PatternMatrix);
//cout<<"Print initial pattern matrix"<<endl;
//for (const auto& row : PatternMatrix) {
//    for (const int& val : row) {
//        std::cout << val << " ";
//    }
//    std::cout << std::endl;
//}

    //If the pattern is empty: use chi2 map to determine the first element.
    bool isEmpty = true;
    for(int i=0; i<NclusV && isEmpty; i++){
      for(int j=0; j<NclusU; j++){
        if(PatternMatrix[i][j]==1) isEmpty = false;
        break;
    }}
    if(isEmpty){
      for(int ipoint=0; ipoint<Chi2MatrixOrder.size(); ipoint++){
        if(PatternMatrix[Chi2MatrixOrder[ipoint].first][Chi2MatrixOrder[ipoint].second]==0 ){
          PatternMatrix[Chi2MatrixOrder[ipoint].first][Chi2MatrixOrder[ipoint].second] = 1;
          break;
        }
      }

//cout<<"Print initial pattern matrix with chi2"<<endl;
//for (const auto& row : PatternMatrix) {
//    for (const int& val : row) {
//        std::cout << val << " ";
//    }
//    std::cout << std::endl;
//}
    }

    //Initialize
    double min_residual = 9999;
    int index = -1;
    std::vector<double> Eij(NclusU*NclusV, 0.);
    std::vector<std::vector<int>> ParameterMatrix;
    std::vector<double> En_clusters; En_clusters.clear();
    std::vector<double> map_residual; map_residual.clear();
    std::vector<vector<vector<int>> > map_pattern; map_pattern.clear();
    std::vector<std::vector<double>> map_Eij; map_Eij.clear();
    for(int icl=0; icl<NclusV; icl++) En_clusters.push_back( m_HFClusVCol[icl]->getEnergy() );
    for(int icl=0; icl<NclusU; icl++) En_clusters.push_back( m_HFClusUCol[icl]->getEnergy() );

    //Loop up to NclusV+NclusU
    for(int iter=0; iter<NclusV+NclusU-1; iter++){
      //Modify the chi2 map with pattern matrix
      bool modified = false;
      bool hasEmptyRow = false;
      bool hasEmptyColumn = false; 
      for (size_t i = 0; i < PatternMatrix.size(); ++i) {
          int rowSum = 0;
          for (size_t j = 0; j < PatternMatrix[i].size(); ++j) {
              rowSum += PatternMatrix[i][j];
          }
          if(rowSum==0){
            hasEmptyRow = true;
            break;
          }
      }
      for (size_t j = 0; j < PatternMatrix[0].size(); ++j) {
          int colSum = 0;
          for (size_t i = 0; i < PatternMatrix.size(); ++i) {
              colSum += PatternMatrix[i][j];
          }
          if(colSum==0){
            hasEmptyColumn = true;
            break;
          }
      }

      vector<vector<double>> Chi2Matrix_copy; 
      if(hasEmptyRow && hasEmptyColumn){

        Chi2Matrix_copy = Chi2Matrix;

        for (size_t i = 0; i < PatternMatrix.size(); ++i) {
            int rowSum = 0;
            for (size_t j = 0; j < PatternMatrix[i].size(); ++j) {
                rowSum += PatternMatrix[i][j];
            }
            if (rowSum != 0) {
                Chi2Matrix_copy[i] = vector<double>(PatternMatrix[0].size(), 9999.);
                modified = true;
            }
        }
        for (size_t j = 0; j < PatternMatrix[0].size(); ++j) {
            int colSum = 0;
            for (size_t i = 0; i < PatternMatrix.size(); ++i) {
                colSum += PatternMatrix[i][j];
            }
            if (colSum != 0) {
                for (size_t i = 0; i < PatternMatrix.size(); ++i) {
                    Chi2Matrix_copy[i][j] = 9999.;
                }
                modified = true;
            }
        }
      }
      else if(hasEmptyRow || hasEmptyColumn){
        Chi2Matrix_copy = vector<vector<double>>(NclusV, vector<double>(NclusU, 9999.));
        for (size_t i = 0; i < PatternMatrix.size(); ++i) {
            int rowSum = 0;
            for (size_t j = 0; j < PatternMatrix[i].size(); ++j) {
                rowSum += PatternMatrix[i][j];
            }
            if (rowSum == 0) {
                Chi2Matrix_copy[i] = Chi2Matrix[i];
                modified = true;
            }
        }

        for (size_t j = 0; j < PatternMatrix[0].size(); ++j) {
            int colSum = 0;
            for (size_t i = 0; i < PatternMatrix.size(); ++i) {
                colSum += PatternMatrix[i][j];
            }
            if (colSum == 0) {
                for (size_t i = 0; i < PatternMatrix.size(); ++i) {
                    Chi2Matrix_copy[i][j] = Chi2Matrix[i][j];
                }
                modified = true;
            }
        }
      }
      else Chi2Matrix_copy = Chi2Matrix;
    

//cout<<"In iteration #"<<iter<<": Print updated chi2 matrix"<<endl;
//for (const auto& row : Chi2Matrix_copy) {
//    for (const double& val : row) {
//        std::cout << val << " ";
//    }
//    std::cout << std::endl;
//}

      std::vector< std::pair<int, int> > Chi2MatrixOrder_copy;
      if(modified){
        std::vector<std::pair<double, int>> indices;
        for(size_t i = 0; i < Chi2Matrix_copy.size(); ++i) {
          for (size_t j = 0; j < Chi2Matrix_copy[i].size(); ++j) {
            indices.push_back(std::make_pair(Chi2Matrix_copy[i][j], i * Chi2Matrix_copy[i].size() + j));
          }
        }
        std::sort(indices.begin(), indices.end());
        for(auto &iter1: indices)
          Chi2MatrixOrder_copy.push_back( std::make_pair( iter1.second/Chi2Matrix_copy[0].size(), iter1.second%Chi2Matrix_copy[0].size() ) );
      }
      else
        Chi2MatrixOrder_copy = Chi2MatrixOrder;


      //Update the pattern matrix with chi2 map
      bool updateFlag = false;
      if(iter>0){
        for(int ipoint=0; ipoint<Chi2MatrixOrder_copy.size(); ipoint++){
          if(PatternMatrix[Chi2MatrixOrder_copy[ipoint].first][Chi2MatrixOrder_copy[ipoint].second]==0 ){
            PatternMatrix[Chi2MatrixOrder_copy[ipoint].first][Chi2MatrixOrder_copy[ipoint].second] = 1;
            updateFlag = true;
            break;
          }
        }
      }
      if(iter>0 && !updateFlag){
        //std::cout<<"ERROR: Already loop for all pattern element. Break. "<<std::endl;
        break;
      }

//cout<<"  Print new pattern matrix "<<endl;
//for (const auto& row : PatternMatrix) {
//    for (const int& val : row) {
//        std::cout << val << " ";
//    }
//    std::cout << std::endl;
//}

      //Re-calculate the parameter matrix, simplify it with pattern matrix
      ParameterMatrix.clear(); Eij.clear();
      Eij = std::vector<double>(NclusU*NclusV, 0.);
      ParameterMatrixCalculation(NclusV, NclusU, ParameterMatrix);
      SimplityMatrix(ParameterMatrix, Eij, PatternMatrix);

      //Recalculate the residual
      double residual = SolveMatrix(NclusU, NclusV, ParameterMatrix, Eij, En_clusters);
//cout<<"  In iteration "<<iter<<": residual = "<<residual<<endl;

      if(residual<min_residual){
        min_residual = residual;
        index = iter;
      }
      map_residual.push_back(residual);
      map_Eij.push_back(Eij);
      map_pattern.push_back(PatternMatrix);
    }

    if(index<0) index=0;


//cout<<"Minimal residual = "<<min_residual<<", index "<<index<<endl;
//cout<<"Pattern map size "<<map_pattern.size()<<", energy map size "<<map_Eij.size()<<endl;
//cout<<"Chosen pattern matrix "<<endl;
//for (const auto& row : map_pattern[index]) {
//    for (const int& val : row) {
//        std::cout << val << " ";
//    }
//    std::cout << std::endl;
//}

    //Derive final energy matrix
    vector<vector<double>> EnergyMatrix(NclusV, vector<double>(NclusU, 0));
    int icount = 0;
    for(int icl=0; icl<NclusV; icl++){
      for(int jcl=0; jcl<NclusU; jcl++){
        if(map_pattern[index][icl][jcl]==1){
          EnergyMatrix[icl][jcl] = map_Eij[index][icount];
          icount++;
        }
      }
    }
//cout<<"Calculated energy matrix"<<endl;
//for (const auto& row : EnergyMatrix) {
//    for (const double& val : row) {
//        std::cout << val << " ";
//    }
//    std::cout << std::endl;
//}

    //Create Calo3DClusters with the result
    ClusterBuilding(EnergyMatrix, m_HFClusUCol, m_HFClusVCol, tmp_clusters);
  }

//cout<<"End matching, clean empty clusters"<<endl;
  //Clean empty Calo3DClusters
  for(int ic=0; ic<tmp_clusters.size(); ic++){
    if( !tmp_clusters[ic].get() || tmp_clusters[ic]->getEnergy()==0 || isnan(tmp_clusters[ic]->getShowerCenter().x()) ){
      tmp_clusters.erase(tmp_clusters.begin()+ic);
      ic--;
    }
  }

  m_clusterCol.insert( m_clusterCol.end(), tmp_clusters.begin(), tmp_clusters.end() );
  tmp_clusters.clear();

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::PatternMatrixCalculation( std::vector<const Cyber::CaloHalfCluster*> m_HFClusUCol,
                                                            std::vector<const Cyber::CaloHalfCluster*> m_HFClusVCol,
                                                            vector<vector<int>>& matrix )
{

  const int NclusU = m_HFClusUCol.size();
  const int NclusV = m_HFClusVCol.size();  

  if(NclusU==0 || NclusV==0){
    std::cout<<"ERROR: empty HalfCluster input. Check! "<<std::endl;
    matrix.clear();
    return StatusCode::SUCCESS;
  }
  else if(NclusU<=1 || NclusV<=1){
    for(int icl=0; icl<NclusV; icl++){
      for(int jcl=0; jcl<NclusU; jcl++) matrix[icl][jcl] = 1;
    }
    return StatusCode::SUCCESS;
  }
  else{
    //Fill the matrix with track info
    std::vector<const Cyber::Track*> m_linkedTrk; m_linkedTrk.clear();

    for(int icl=0; icl<m_HFClusUCol.size(); ++icl){
      std::vector<const Cyber::Track*> tmp_linkedTrk = m_HFClusUCol[icl]->getAssociatedTracks();
      if(tmp_linkedTrk.size()==0) continue;
      m_linkedTrk.insert(m_linkedTrk.end(), tmp_linkedTrk.begin(), tmp_linkedTrk.end() );
    }
    for(int icl=0; icl<m_HFClusVCol.size(); ++icl){
      std::vector<const Cyber::Track*> tmp_linkedTrk = m_HFClusVCol[icl]->getAssociatedTracks();
      if(tmp_linkedTrk.size()==0) continue;
      m_linkedTrk.insert(m_linkedTrk.end(), tmp_linkedTrk.begin(), tmp_linkedTrk.end());
    }

    for(int itrk=0; itrk<m_linkedTrk.size(); itrk++){
      int indexU = -1;
      int indexV = -1;
      for(int icl=0; icl<m_HFClusUCol.size(); ++icl){
        std::vector<const Cyber::Track*> m_linkedTrkU = m_HFClusUCol[icl]->getAssociatedTracks();
        if(m_linkedTrkU.size()==0) continue;
        if(find(m_linkedTrkU.begin(), m_linkedTrkU.end(), m_linkedTrk[itrk])!=m_linkedTrkU.end()){
          indexU = icl;
          break;
        }
      }

      for(int icl=0; icl<m_HFClusVCol.size(); ++icl){
        std::vector<const Cyber::Track*> m_linkedTrkV = m_HFClusVCol[icl]->getAssociatedTracks();
        if(m_linkedTrkV.size()==0) continue;
        if(find(m_linkedTrkV.begin(), m_linkedTrkV.end(), m_linkedTrk[itrk])!=m_linkedTrkV.end()){
          indexV = icl;
          break;
        }
      }

      if(indexU>=0 && indexV>=0) matrix[indexV][indexU] = 1;
    }

    //Fill the matrix with neighbor info
    for(int icl=0; icl<m_HFClusUCol.size(); ++icl){
      if(m_HFClusUCol[icl]->getHalfClusterCol("CousinCluster").size()==0) continue;
      for(int jcl=0; jcl<m_HFClusVCol.size(); ++jcl){
        if(m_HFClusVCol[jcl]->getHalfClusterCol("CousinCluster").size()==0) continue;

        std::vector<const Cyber::CaloHalfCluster*> m_cousinU = m_HFClusUCol[icl]->getHalfClusterCol("CousinCluster");
        std::vector<const Cyber::CaloHalfCluster*> m_cousinV = m_HFClusVCol[jcl]->getHalfClusterCol("CousinCluster");

        //Make map: <towerID, (E_clU - E_clV)/mean(E_clU, E_clV)  >
        std::map<std::vector<int>, pair<float, float>> map_pairE; map_pairE.clear();
        map_pairE[m_HFClusUCol[icl]->getTowerID()[0]] = make_pair(m_HFClusUCol[icl]->getEnergy(), m_HFClusVCol[jcl]->getEnergy());
        for(int ics=0; ics<m_cousinU.size(); ++ics){
        for(int jcs=0; jcs<m_cousinV.size(); ++jcs){
          if(m_cousinU[ics]->getTowerID()[0] == m_cousinV[jcs]->getTowerID()[0])
            map_pairE[ m_cousinU[ics]->getTowerID()[0] ] = make_pair(m_cousinU[ics]->getEnergy(), m_cousinV[jcs]->getEnergy());

        }}
        if(map_pairE.size()<=1) continue; //No common tower cousin clusters.

        bool isLink = true;
        float totE_U = 0.;
        float totE_V = 0.;
        for(auto iter: map_pairE){
          float deltaE = 2*fabs(iter.second.first - iter.second.second)/(iter.second.first + iter.second.second);
          totE_U += iter.second.first;
          totE_V += iter.second.second;
          //In tower: if too large difference between U/V then do not link.
          if(deltaE > settings.map_floatPars["th_UVdeltaE"]) { isLink = false; break; }
          //if(iter.second.first>1 && iter.second.second>1 && deltaE > settings.map_floatPars["th_UVdeltaE1"]){ isLink = false; break; } //High energy: deltaE
          //else if(iter.second.first/iter.second.second > 10 || iter.second.first/iter.second.second<0.1) { isLink = false; break; }  //Low energy: 0.1<E_u/E_v<10.
        }
        if( 2*fabs(totE_U-totE_V)/(totE_U+totE_V)>settings.map_floatPars["th_UVdeltaE"] ) isLink=false;

        if(isLink) matrix[jcl][icl] = 1;
      }
    }
  }

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::Chi2MatrixCalculation( std::vector<const Cyber::CaloHalfCluster*> m_HFClusUCol,
                                                         std::vector<const Cyber::CaloHalfCluster*> m_HFClusVCol,
                                                         vector<vector<double>>& matrix,
                                                         std::vector< std::pair<int, int> >& chi2order)
{

  const int NclusU = m_HFClusUCol.size();
  const int NclusV = m_HFClusVCol.size();
//cout<<"  Chi2MatrixCalculation: input cluster size "<<NclusU<<", "<<NclusV<<endl;

  if(NclusU==0 || NclusV==0){
    std::cout<<"ERROR: empty HalfCluster input. Check! "<<std::endl;
    matrix.clear();
    chi2order.clear();
    return StatusCode::SUCCESS;
  }
  else if(NclusU<=1 || NclusV<=1){
    for(int icl=0; icl<NclusV; icl++){
      for(int jcl=0; jcl<NclusU; jcl++) matrix[icl][jcl] = 0.;
    }
    chi2order.clear();
    return StatusCode::SUCCESS;
  }
  else{ 
//cout<<" Print sum chi2 matrix"<<endl;
//for (const auto& row : matrix) {
//    for (const double& val : row) {
//        std::cout << val << " ";
//    }
//    std::cout << std::endl;
//}

    int Nlayer = Cyber::CaloUnit::Nlayer; 
    for(int il=0; il<Nlayer; il++){
      std::vector<std::vector<const Cyber::Calo1DCluster*>> tmp_showerU(NclusU);
      std::vector<std::vector<const Cyber::Calo1DCluster*>> tmp_showerV(NclusV);
      for(int icl=0; icl<NclusU; icl++){
        for(int jcl=0; jcl<NclusV; jcl++){
          tmp_showerU[icl] = m_HFClusUCol[icl]->getClusterInLayer(il+1);
          tmp_showerV[jcl] = m_HFClusVCol[jcl]->getClusterInLayer(il+1);
        }
      }
//cout<<"  In layer #"<<il<<": cluster size "<<tmp_showerU.size()<<", "<<tmp_showerV.size()<<endl;
      vector<vector<double>> chi2Map;
      chi2Map = GetClusterChi2Map(tmp_showerU, tmp_showerV);

//cout<<" Print chi2 matrix"<<endl;
//for (const auto& row : chi2Map) {
//    for (const double& val : row) {
//        std::cout << val << " ";
//    }
//    std::cout << std::endl;
//}

      for(int icl=0; icl<NclusU; icl++){
        for(int jcl=0; jcl<NclusV; jcl++){
          matrix[jcl][icl] += chi2Map[jcl][icl];
          //Use 1/E as weight: 
          //matrix[jcl][icl] += chi2Map[jcl][icl]*( 1./m_HFClusUCol[icl]->getEnergy() + 1./m_HFClusVCol[jcl]->getEnergy() );
        }
      }
//cout<<" Print sum chi2 matrix"<<endl;
//for (const auto& row : matrix) {
//    for (const double& val : row) {
//        std::cout << val << " ";
//    }
//    std::cout << std::endl;
//}
    }

    
    //Get the order of chi2 matrix map
    std::vector<std::pair<double, int>> indices;
    for(size_t i = 0; i < matrix.size(); ++i) {
      for (size_t j = 0; j < matrix[i].size(); ++j) {
        indices.push_back(std::make_pair(matrix[i][j], i * matrix[i].size() + j));
      }
    }
    std::sort(indices.begin(), indices.end());
//cout<<"  matrix size: "<<matrix.size()<<", "<<matrix[0].size()<<endl;
    for(auto &iter: indices)
      chi2order.push_back( std::make_pair( iter.second/matrix[0].size(), iter.second%matrix[0].size() ) );

  }

  return StatusCode::SUCCESS;
}

vector<vector<double>> EnergyTimeMatchingAlg::GetClusterChi2Map( std::vector<std::vector<const Cyber::Calo1DCluster*>>& barShowerUCol,
                                                                 std::vector<std::vector<const Cyber::Calo1DCluster*>>& barShowerVCol )
{


  const int NclusU = barShowerUCol.size();
  const int NclusV = barShowerVCol.size();
  vector<vector<double>> chi2map_E(NclusV, vector<double>(NclusU, 0.));
  vector<vector<double>> chi2map_tx(NclusV, vector<double>(NclusU, 0.));
  vector<vector<double>> chi2map_ty(NclusV, vector<double>(NclusU, 0.));
  vector<vector<double>> chi2map(NclusV, vector<double>(NclusU, 0.));

  if(NclusU==0 || NclusV==0) return chi2map;

  double wi_E = settings.map_floatPars["chi2Wi_E"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);
  double wi_T = settings.map_floatPars["chi2Wi_T"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);

  //Rotate angle and tower center position for ECAL Barrel
  TVector3 m_vec(0,0,0);
  double rotAngle = -999;
  TVector3 Ctower(0,0,0);
  int system = -1;
  for(int ish=0; ish<barShowerUCol.size(); ish++){
    if(barShowerUCol[ish].size()==0) continue;
    rotAngle = -(barShowerUCol[ish][0]->getBars())[0]->getModule()*TMath::TwoPi()/Cyber::CaloUnit::Nmodule;
    Ctower.SetX( (barShowerUCol[ish][0]->getBars())[0]->getPosition().x() );
    Ctower.SetY( (barShowerUCol[ish][0]->getBars())[0]->getPosition().y() );
    system = (barShowerUCol[ish][0]->getBars())[0]->getSystem();
    break;
  }
  for(int ish=0; ish<barShowerVCol.size(); ish++){
    if(barShowerVCol[ish].size()==0) continue;
    Ctower.SetZ( (barShowerVCol[ish][0]->getBars())[0]->getPosition().z() );
    system = (barShowerVCol[ish][0]->getBars())[0]->getSystem();
    break;
  }
  Ctower.RotateZ(rotAngle);

  for(int ix=0;ix<NclusU;ix++){
  for(int iy=0;iy<NclusV;iy++){
    std::vector<const Cyber::Calo1DCluster*> clusterU = barShowerUCol[ix];
    std::vector<const Cyber::Calo1DCluster*> clusterV = barShowerVCol[iy];

    if(clusterU.size()==0){
      double totE_V = 0;
      for(int icy=0; icy<clusterV.size(); icy++) totE_V += clusterV[icy]->getEnergy();
      chi2map[iy][ix] = pow(totE_V/settings.map_floatPars["sigmaE"], 2);
      continue;
    }
    if(clusterV.size()==0){
      double totE_U = 0;
      for(int icy=0; icy<clusterU.size(); icy++) totE_U += clusterU[icy]->getEnergy();
      chi2map[iy][ix] = pow(totE_U/settings.map_floatPars["sigmaE"], 2);
      continue;
    }

    double min_chi2E = 999;
    double min_chi2tx = 999;
    double min_chi2ty = 999;

    for(int icx=0; icx<clusterU.size(); icx++){
    for(int icy=0; icy<clusterV.size(); icy++){
      const Cyber::Calo1DCluster* showerX = clusterU[icx];
      const Cyber::Calo1DCluster* showerY = clusterV[icy];

      double Ex = showerX->getEnergy();
      double Ey = showerY->getEnergy();
      double chi2_E = pow(fabs(Ex-Ey)/settings.map_floatPars["sigmaE"], 2);

      double PosTx, chi2_tx, PosTy, chi2_ty;
      if(system == Cyber::CaloUnit::System_Barrel){
        PosTx = C*(showerY->getT1()-showerY->getT2())/(2*settings.map_floatPars["nMat"]) + showerY->getPos().z();
        chi2_tx = pow( fabs(PosTx-showerX->getPos().z())/settings.map_floatPars["sigmaPos"], 2 );

        PosTy = C*(showerX->getT1()-showerX->getT2())/(2*settings.map_floatPars["nMat"]);
        m_vec = showerY->getPos();
        m_vec.RotateZ(rotAngle);
        chi2_ty = pow( fabs(PosTy - (m_vec-Ctower).x() )/settings.map_floatPars["sigmaPos"], 2);
      }
      else if(system == Cyber::CaloUnit::System_Endcap){
        PosTx = C*(showerY->getT1()-showerY->getT2())/(2*settings.map_floatPars["nMat"]) + showerY->getPos().x();
        chi2_tx = pow( fabs(PosTx-showerX->getPos().x())/settings.map_floatPars["sigmaPos"], 2 );

        PosTy = C*(showerX->getT1()-showerX->getT2())/(2*settings.map_floatPars["nMat"]) + showerX->getPos().y();
        chi2_ty = pow( fabs(PosTy - showerY->getPos().y() )/settings.map_floatPars["sigmaPos"], 2);
      }
      else{
        std::cout<<"ERROR in EnergyTimeMatchingAlg: Unknown system ID: "<<system<<endl;
      }

      if(chi2_E<min_chi2E) min_chi2E=chi2_E;
      if(chi2_tx<min_chi2tx) min_chi2tx=chi2_tx;
      if(chi2_ty<min_chi2ty) min_chi2ty=chi2_ty;
    }}

    chi2map_E[iy][ix] = min_chi2E;
    chi2map_tx[iy][ix] = min_chi2tx;
    chi2map_ty[iy][ix] = min_chi2ty;
    chi2map[iy][ix] = chi2map_E[iy][ix]*wi_E + (chi2map_tx[iy][ix]+chi2map_ty[iy][ix])*wi_T ;
//cout<<"  chi2 value "<<min_chi2E<<endl;

  }}

  return chi2map;
}


StatusCode EnergyTimeMatchingAlg::ParameterMatrixCalculation(const int& M, const int& N, vector<vector<int>>& parMatrix){

  parMatrix = vector<vector<int>>(M + N, std::vector<int>(M * N, 0));
  //Fill first M raws
  for (int i = 0; i < M; i++) {
      for (int j = i * N; j < (i + 1) * N; j++) {
          parMatrix[i][j] = 1;
      }
  }
  //Fill second N raws
  for(int i=M; i<N+M; i++){
    for(int j=0; j<M; j++)
      parMatrix[i][i-M+j*N] = 1;
  }

  return StatusCode::SUCCESS; 
}


StatusCode EnergyTimeMatchingAlg::SimplityMatrix(vector<vector<int>>& parMatrix, vector<double>& Eij, vector<vector<int>>& pattern){
  int M = pattern.size();
  int N = pattern[0].size();

  int count = 0;
  for(int i=0; i<M; i++){
    for(int j=0; j<N; j++){
      int k = i*N+j-count;
      if(pattern[i][j]==0){
          for (int a = 0; a < M + N; a++) {
            parMatrix[a].erase(parMatrix[a].begin() + k);
          }
          Eij.erase(Eij.begin()+k);
          count++;
      }
    }
  }

  return StatusCode::SUCCESS;
}


double EnergyTimeMatchingAlg::SolveMatrix(const int& M, const int& N, vector<vector<int>>& parMatrix, vector<double>& Eij, vector<double>& En_clusters){


  leastSquares(parMatrix, En_clusters, Eij);
  vector<double> residual; residual.clear();

  for (int i = 0; i < parMatrix.size(); i++) {
      double r = 0;
      for (int j = 0; j < parMatrix[i].size(); j++) {
          r += parMatrix[i][j] * Eij[j];
      }
      r -= En_clusters[i];
      residual.push_back(r);
  }

  double residualSum = 0;
  for (int i = 0; i < residual.size(); i++) {
      residualSum += residual[i] * residual[i];
  }

  return residualSum;
}


StatusCode EnergyTimeMatchingAlg::leastSquares(const vector<vector<int>>& A, const vector<double>& b, vector<double>& x){
  int m = A.size();
  int n = A[0].size();

  vector<vector<double>> Q, R;
  vector<vector<double>> doubleA(m, vector<double>(n));
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++)
      doubleA[i][j] = (double)A[i][j];
  }

  //qrDecomposition(doubleA, Q, R);
  Q = vector<vector<double>>(m, vector<double>(n, 0));
  R = vector<vector<double>>(n, vector<double>(n, 0));
  vector<vector<double>> A_copy = doubleA;
  for (int j = 0; j < n; j++) {
      double colNorm = 0.0;
      for (int i = 0; i < m; i++) {
          colNorm += A_copy[i][j] * A_copy[i][j];
      }
      colNorm = sqrt(colNorm);

      R[j][j] = colNorm;

      for (int i = 0; i < m; i++) {
          Q[i][j] = A_copy[i][j] / colNorm;
      }

      for (int k = j + 1; k < n; k++) {
          double dotProduct = 0.0;
          for (int i = 0; i < m; i++) {
              dotProduct += Q[i][j] * A_copy[i][k];
          }
          R[j][k] = dotProduct;
          for (int i = 0; i < m; i++) {
              A_copy[i][k] -= Q[i][j] * dotProduct;
          }
      }
  }


  vector<double> Qt_b(n, 0);
  for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
          Qt_b[j] += Q[i][j] * b[i];
      }
  }

  x = vector<double>(n, 0);
  for (int i = n - 1; i >= 0; i--) {
      double sum = 0.0;
      for (int j = i + 1; j < n; j++) {
          sum += R[i][j] * x[j];
      }
      if(R[i][i]!=0) x[i] = (Qt_b[i] - sum) / R[i][i];
      if(x[i]<0) x[i] = 0.;
  }

  return StatusCode::SUCCESS;
}

StatusCode EnergyTimeMatchingAlg::ClusterBuilding( std::vector<std::vector<double>> Ematrix,
                                                   std::vector<const Cyber::CaloHalfCluster*>& m_HFClusUCol,
                                                   std::vector<const Cyber::CaloHalfCluster*>& m_HFClusVCol,
                                                   std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_clusterCol)
{

  int NclusU = m_HFClusUCol.size();
  int NclusV = m_HFClusVCol.size();
  if(NclusU==0 || NclusV==0){
    std::cout<<"ERROR: Empty readin HalfCluster! "<<NclusU<<", "<<NclusV<<endl;
  }
  else if(NclusU==1 && NclusV==1){
    std::shared_ptr<Cyber::Calo3DCluster> tmp_clus = std::make_shared<Cyber::Calo3DCluster>();
    XYClusterMatchingL0(m_HFClusUCol[0], m_HFClusVCol[0], tmp_clus);
    m_clusterCol.push_back(tmp_clus);
  }
  else{
    std::vector<double> sumE_U(NclusU, 0.); //column
    std::vector<double> sumE_V(NclusV, 0.); //raw
    for(int icl=0; icl<NclusV; icl++){
      for(int jcl=0; jcl<NclusU; jcl++){
        sumE_V[icl] += Ematrix[icl][jcl];
        sumE_U[jcl] += Ematrix[icl][jcl];
      }
    }

    for(int icl=0; icl<NclusV; icl++){
      for(int jcl=0; jcl<NclusU; jcl++){

        double fracU = Ematrix[icl][jcl]/sumE_U[jcl];
        if(fracU==0) continue;
        std::shared_ptr<Cyber::CaloHalfCluster> newClusU = std::make_shared<Cyber::CaloHalfCluster>();
        if(fracU==1) newClusU = m_HFClusUCol[jcl]->Clone();
        else{
          for(int ish=0; ish<m_HFClusUCol[jcl]->getCluster().size(); ish++){
            const Cyber::Calo1DCluster* p_shower = m_HFClusUCol[jcl]->getCluster()[ish];
   
            std::vector<const Cyber::CaloUnit*> Bars; Bars.clear();
            for(int ibar=0; ibar<p_shower->getCluster().size(); ibar++){
              auto bar = p_shower->getCluster()[ibar]->Clone();
              bar->setQ(bar->getQ1()*fracU, bar->getQ2()*fracU );
              Bars.push_back(bar.get());
              m_bkCol.map_BarCol["bkBar"].push_back(bar);
            }
   
            std::shared_ptr<Cyber::Calo1DCluster> shower = std::make_shared<Cyber::Calo1DCluster>();
            shower->setBars(Bars);
            shower->setSeed();
            shower->setIDInfo();
   
            newClusU->addUnit(shower.get());
            m_bkCol.map_1DCluster["bk1DCluster"].push_back( shower );
          }
          newClusU->setHoughPars(m_HFClusUCol[jcl]->getHoughAlpha(), m_HFClusUCol[jcl]->getHoughRho());
          newClusU->setIntercept(m_HFClusUCol[jcl]->getHoughIntercept());
          for(auto iter: m_HFClusUCol[jcl]->getHalfClusterMap()) newClusU->setHalfClusters( iter.first, iter.second );
          for(int itrk=0; itrk<m_HFClusUCol[jcl]->getAssociatedTracks().size(); itrk++) newClusU->addAssociatedTrack(m_HFClusUCol[jcl]->getAssociatedTracks()[itrk]);
          newClusU->setLinkedMCP(m_HFClusUCol[jcl]->getLinkedMCP());
          newClusU->setType(m_HFClusUCol[jcl]->getType());
        }

        std::shared_ptr<Cyber::CaloHalfCluster> newClusV = std::make_shared<Cyber::CaloHalfCluster>();
        double fracV = Ematrix[icl][jcl]/sumE_V[icl];
        if(fracV==0) continue;
        if(fracV==1) newClusV = m_HFClusVCol[icl]->Clone();
        else{
          for(int ish=0; ish<m_HFClusVCol[icl]->getCluster().size(); ish++){
            const Cyber::Calo1DCluster* p_shower = m_HFClusVCol[icl]->getCluster()[ish];
   
            std::vector<const Cyber::CaloUnit*> Bars; Bars.clear();
            for(int ibar=0; ibar<p_shower->getCluster().size(); ibar++){
              auto bar = p_shower->getCluster()[ibar]->Clone();
              bar->setQ(bar->getQ1()*fracV, bar->getQ2()*fracV );
              Bars.push_back(bar.get());
              m_bkCol.map_BarCol["bkBar"].push_back(bar);
            }
   
            std::shared_ptr<Cyber::Calo1DCluster> shower = std::make_shared<Cyber::Calo1DCluster>();
            shower->setBars(Bars);
            shower->setSeed();
            shower->setIDInfo();
   
            newClusV->addUnit(shower.get());
            m_bkCol.map_1DCluster["bk1DCluster"].push_back( shower );
          }
          newClusV->setHoughPars(m_HFClusVCol[icl]->getHoughAlpha(), m_HFClusVCol[icl]->getHoughRho());
          newClusV->setIntercept(m_HFClusVCol[icl]->getHoughIntercept());
          for(auto iter: m_HFClusVCol[icl]->getHalfClusterMap()) newClusV->setHalfClusters( iter.first, iter.second );
          for(int itrk=0; itrk<m_HFClusVCol[icl]->getAssociatedTracks().size(); itrk++) newClusV->addAssociatedTrack(m_HFClusVCol[icl]->getAssociatedTracks()[itrk]);
          newClusV->setLinkedMCP(m_HFClusVCol[icl]->getLinkedMCP());
          newClusV->setType(m_HFClusVCol[icl]->getType());
        }
        m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(newClusU);
        m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(newClusV);

        std::shared_ptr<Cyber::Calo3DCluster> tmp_clus = std::make_shared<Cyber::Calo3DCluster>();
        XYClusterMatchingL0(newClusU.get(), newClusV.get(), tmp_clus);
        m_clusterCol.push_back(tmp_clus);
    }}
  }

  return StatusCode::SUCCESS;
}


//Longitudinal cluster: 1*1
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL0( const Cyber::CaloHalfCluster* m_longiClU,
                                                       const Cyber::CaloHalfCluster* m_longiClV,
                                                       std::shared_ptr<Cyber::Calo3DCluster>& m_clus )
{
//cout<<"  Cluster matching for case: 1 * 1. Input HalfCluster En: "<<m_longiClU->getEnergy()<<", "<<m_longiClV->getEnergy()<<endl;
//cout<<"  Print 1DShower En in HalfClusterU: "<<endl;
//for(int i=0; i<m_longiClU->getCluster().size(); i++)
//  cout<<m_longiClU->getCluster()[i]->getDlayer()<<'\t'<<m_longiClU->getCluster()[i]->getEnergy()<<endl;
//cout<<"  Print 1DShower En in HalfClusterV: "<<endl;
//for(int i=0; i<m_longiClV->getCluster().size(); i++)
//  cout<<m_longiClV->getCluster()[i]->getDlayer()<<'\t'<<m_longiClV->getCluster()[i]->getEnergy()<<endl;


  std::vector<int> layerindex; layerindex.clear();
  std::map<int, std::vector<const Cyber::Calo1DCluster*> > map_showersUinlayer; map_showersUinlayer.clear();
  std::map<int, std::vector<const Cyber::Calo1DCluster*> > map_showersVinlayer; map_showersVinlayer.clear();

  for(int is=0; is<m_longiClU->getCluster().size(); is++){
    int m_layer = m_longiClU->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
    map_showersUinlayer[m_layer].push_back(m_longiClU->getCluster()[is]);
  }
  for(int is=0; is<m_longiClV->getCluster().size(); is++){
    int m_layer = m_longiClV->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
    map_showersVinlayer[m_layer].push_back(m_longiClV->getCluster()[is]);
  }

  for(int il=0; il<layerindex.size(); il++){
    std::vector<const Cyber::Calo1DCluster*> m_showerXcol = map_showersUinlayer[layerindex[il]];
    std::vector<const Cyber::Calo1DCluster*> m_showerYcol = map_showersVinlayer[layerindex[il]];

//cout<<"  Print 1D showers in layer "<<layerindex[il]<<endl;
//cout<<"  ShowerU size = "<<m_showerXcol.size()<<endl;
//for(int a=0; a<m_showerXcol.size(); a++) printf("    #%d shower: En %.3f, address %p \n", a, m_showerXcol[a]->getEnergy(), m_showerXcol[a]);
//cout<<"  ShowerV size = "<<m_showerYcol.size()<<endl;
//for(int a=0; a<m_showerYcol.size(); a++) printf("    #%d shower: En %.3f, address %p \n", a, m_showerYcol[a]->getEnergy(), m_showerYcol[a]);

    std::vector<Cyber::Calo2DCluster*> m_showerinlayer; m_showerinlayer.clear();

    if(m_showerXcol.size()==0 || m_showerYcol.size()==0) continue;
    //else if(m_showerXcol.size()==0){
    //  std::shared_ptr<Cyber::Calo2DCluster> tmp_shower = std::make_shared<Cyber::Calo2DCluster>();
    //  GetMatchedShowerFromEmpty(m_showerYcol[0], m_longiClU, tmp_shower.get());
    //  //m_showerinlayer.push_back(tmp_shower.get());
    //  //m_bkCol.map_2DCluster["bk2DCluster"].push_back(tmp_shower);
    //}
    //else if(m_showerYcol.size()==0){
    //  std::shared_ptr<Cyber::Calo2DCluster> tmp_shower = std::make_shared<Cyber::Calo2DCluster>();
    //  GetMatchedShowerFromEmpty(m_showerXcol[0], m_longiClU, tmp_shower.get());
    //}
    else if(m_showerXcol.size()==1 && m_showerYcol.size()==1){
      std::shared_ptr<Cyber::Calo2DCluster> tmp_shower = std::make_shared<Cyber::Calo2DCluster>();
      GetMatchedShowersL0(m_showerXcol[0], m_showerYcol[0], tmp_shower.get());
      m_showerinlayer.push_back(tmp_shower.get());
      m_bkCol.map_2DCluster["bk2DCluster"].push_back(tmp_shower);
    }
    else if(m_showerXcol.size()==1) GetMatchedShowersL1(m_showerXcol[0], m_showerYcol, m_showerinlayer );
    else if(m_showerYcol.size()==1) GetMatchedShowersL1(m_showerYcol[0], m_showerXcol, m_showerinlayer );
    else{ std::cout<<"CAUSION in XYClusterMatchingL0: HFCluster has ["<<m_showerXcol.size()<<", "<<m_showerYcol.size()<<"] showers in layer "<<layerindex[il]<<std::endl; }

//cout<<"    After matching: shower size = "<<m_showerinlayer.size()<<", Print showers: "<<endl;
//for(int is=0; is<m_showerinlayer.size(); is++) printf("  Pos+E (%.3f, %.3f, %.3f, %.3f) \t", m_showerinlayer[is]->getPos().x(), m_showerinlayer[is]->getPos().y(), m_showerinlayer[is]->getPos().z(), m_showerinlayer[is]->getEnergy() );
//cout<<endl;

    for(int is=0; is<m_showerinlayer.size(); is++) m_clus->addUnit(m_showerinlayer[is]);
  }

  m_clus->setCaloHitsFrom2DCluster();
  m_clus->addHalfClusterU( "LinkedLongiCluster", m_longiClU );
  m_clus->addHalfClusterV( "LinkedLongiCluster", m_longiClV );
  for(auto itrk : m_longiClU->getAssociatedTracks()){
    for(auto jtrk : m_longiClV->getAssociatedTracks()){
      if(itrk!=jtrk) continue;
      if( find(m_clus->getAssociatedTracks().begin(), m_clus->getAssociatedTracks().end(), itrk)==m_clus->getAssociatedTracks().end() )
        m_clus->addAssociatedTrack(itrk);
  }}
  m_clus->getLinkedMCPfromHFCluster("LinkedLongiCluster");

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL0( const Cyber::Calo1DCluster* barShowerU,
                                                       const Cyber::Calo1DCluster* barShowerV,
                                                       Cyber::Calo2DCluster* outsh )
{

  std::vector<const Cyber::CaloHit*> m_digiCol; m_digiCol.clear();
  int NbarsX = barShowerU->getBars().size();
  int NbarsY = barShowerV->getBars().size();
  if(NbarsX==0 || NbarsY==0){ std::cout<<"WARNING: empty DigiHitsCol returned!"<<std::endl; return StatusCode::SUCCESS; }
  if(barShowerU->getTowerID().size()==0) { std::cout<<"WARNING:GetMatchedShowersL0  No TowerID in 1DCluster!"<<std::endl; return StatusCode::SUCCESS; }
  //if(barShowerU->getTowerID().size()==0) { barShowerU->setIDInfo(); }

  int _layer = barShowerU->getDlayer();
  int _module = barShowerU->getTowerID()[0][1];
  int _system = barShowerU->getTowerID()[0][0];
  if(_system == Cyber::CaloUnit::System_Barrel){
    float rotAngle = -_module*TMath::TwoPi()/Cyber::CaloUnit::Nmodule;

    for(int ibar=0;ibar<NbarsX;ibar++){
      double En_x = barShowerU->getBars()[ibar]->getEnergy();
      TVector3 m_vecx = barShowerU->getBars()[ibar]->getPosition();
      m_vecx.RotateZ(rotAngle);

      for(int jbar=0;jbar<NbarsY;jbar++){
        double En_y = barShowerV->getBars()[jbar]->getEnergy();
        TVector3 m_vecy = barShowerV->getBars()[jbar]->getPosition();
        m_vecy.RotateZ(rotAngle);

        TVector3 p_hit(m_vecy.x(), (m_vecx.y()+m_vecy.y())/2., m_vecx.z() );
        p_hit.RotateZ(-rotAngle);
        double m_Ehit = En_x*En_y/barShowerV->getEnergy() + En_x*En_y/barShowerU->getEnergy();
        //Create new CaloHit
        std::shared_ptr<Cyber::CaloHit> hit = std::make_shared<Cyber::CaloHit>();
        hit->setcellID(_system, _module, _layer);
        hit->setPosition(p_hit);
        hit->setEnergy(m_Ehit);
        m_digiCol.push_back(hit.get());
        m_bkCol.map_CaloHit["bkHit"].push_back( hit );
      }
    }
  }
  else{

    for(int ibar=0;ibar<NbarsX;ibar++){
      double En_x = barShowerU->getBars()[ibar]->getEnergy();
      TVector3 m_vecx = barShowerU->getBars()[ibar]->getPosition();

      for(int jbar=0;jbar<NbarsY;jbar++){
        double En_y = barShowerV->getBars()[jbar]->getEnergy();
        TVector3 m_vecy = barShowerV->getBars()[jbar]->getPosition();

        TVector3 p_hit(m_vecx.x(), m_vecy.y(), (m_vecx.z()+m_vecy.z())/2. );
        double m_Ehit = En_x*En_y/barShowerV->getEnergy() + En_x*En_y/barShowerU->getEnergy();
        //Create new CaloHit
        std::shared_ptr<Cyber::CaloHit> hit = std::make_shared<Cyber::CaloHit>();
        hit->setcellID(_system, _module, _layer);
        hit->setPosition(p_hit);
        hit->setEnergy(m_Ehit);
        m_digiCol.push_back(hit.get());
        m_bkCol.map_CaloHit["bkHit"].push_back( hit );
      }
    }

  }

  outsh->addUnit( barShowerU );
  outsh->addUnit( barShowerV );
  outsh->setCaloHits( m_digiCol );
//cout<<"    End output shower"<<endl;

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL1( const Cyber::Calo1DCluster* shower1,
                                                       std::vector<const Cyber::Calo1DCluster*>& showerNCol,
                                                       std::vector<Cyber::Calo2DCluster*>& outshCol )
{
//cout<<"  GetMatchedShowersL1: input shower size: 1 * "<<showerNCol.size()<<endl;
  outshCol.clear();

  int _slayer = shower1->getBars()[0]->getSlayer();

  const int NshY = showerNCol.size();
  double totE_shY = 0;
  double EshY[NshY] = {0};
  for(int is=0;is<NshY;is++){ EshY[is] = showerNCol[is]->getEnergy(); totE_shY += EshY[is]; }
  for(int is=0;is<NshY;is++){
    double wi_E = EshY[is]/totE_shY;
    std::shared_ptr<Cyber::Calo1DCluster> m_splitshower1 = std::make_shared<Cyber::Calo1DCluster>();
    m_bkCol.map_1DCluster["bk1DCluster"].push_back( m_splitshower1 );

    std::shared_ptr<Cyber::CaloUnit> m_wiseed = nullptr;
    if(shower1->getSeeds().size()>0) m_wiseed = shower1->getSeeds()[0]->Clone();
    else{ cout<<"ERROR: Input shower has no seed! Check! Use the most energitic bar as seed. bar size: "<<shower1->getBars().size()<<endl;
      double m_maxE = -99;
      int index = -1;
      for(int ib=0; ib<shower1->getBars().size(); ib++){
        if(shower1->getBars()[ib]->getEnergy()>m_maxE) { m_maxE=shower1->getBars()[ib]->getEnergy(); index=ib; }
      }
      if(index>=0) m_wiseed = shower1->getBars()[index]->Clone();
    }
    m_wiseed->setQ( wi_E*m_wiseed->getQ1(), wi_E*m_wiseed->getQ2() );
    m_bkCol.map_BarCol["bkBar"].push_back( m_wiseed );

    std::vector<const Cyber::CaloUnit*> m_wibars; m_wibars.clear();
    for(int ib=0;ib<shower1->getBars().size();ib++){
      std::shared_ptr<Cyber::CaloUnit> m_wibar = shower1->getBars()[ib]->Clone();
      m_wibar->setQ(wi_E*m_wibar->getQ1(), wi_E*m_wibar->getQ2());
      m_wibars.push_back(m_wibar.get());
      m_bkCol.map_BarCol["bkBar"].push_back( m_wibar );
    }
    m_splitshower1->setBars( m_wibars );
    m_splitshower1->addSeed( m_wiseed.get() );
    m_splitshower1->setIDInfo();
    std::shared_ptr<Cyber::Calo2DCluster> m_shower = std::make_shared<Cyber::Calo2DCluster>();
    if(_slayer==0 ) GetMatchedShowersL0( m_splitshower1.get(), showerNCol[is], m_shower.get() );
    else            GetMatchedShowersL0( showerNCol[is], m_splitshower1.get(), m_shower.get() );

    outshCol.push_back( m_shower.get() );
    m_bkCol.map_2DCluster["bk2DCluster"].push_back( m_shower );
  }
  return StatusCode::SUCCESS;
}




StatusCode EnergyTimeMatchingAlg::ClusterReconnecting( std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_clusterCol ){
/*
cout<<"Print 3DClusters"<<endl;
for(int ic=0; ic<m_clusterCol.size(); ic++){
  printf("  Cluster #%d pos+En [%.3f, %.3f, %.3f, %.3f], hit size %d \n ", ic, m_clusterCol[ic]->getShowerCenter().x(), 
                                                          m_clusterCol[ic]->getShowerCenter().y(), 
                                                          m_clusterCol[ic]->getShowerCenter().z(), 
                                                          m_clusterCol[ic]->getEnergy(), 
                                                          m_clusterCol[ic]->getCaloHits().size() );
  std::vector<const CaloHalfCluster*> m_HFClusUInTower = m_clusterCol[ic]->getHalfClusterUCol("LinkedLongiCluster");
  std::vector<const CaloHalfCluster*> m_HFClusVInTower = m_clusterCol[ic]->getHalfClusterVCol("LinkedLongiCluster");  
  cout<<"  Linked half cluster size: "<<m_HFClusUInTower.size()<<", "<<m_HFClusVInTower.size()<<endl;
  cout<<"    Loop print HalfClusterU: "<<endl;
  for(int icl=0; icl<m_HFClusUInTower.size(); icl++){
    cout<<"      In HFClusU #"<<icl<<": shower size = "<<m_HFClusUInTower[icl]->getCluster().size()<<", En = "<<m_HFClusUInTower[icl]->getEnergy()<<", type "<<m_HFClusUInTower[icl]->getType();
    printf(", Position (%.3f, %.3f, %.3f), address %p ",m_HFClusUInTower[icl]->getPos().x(), m_HFClusUInTower[icl]->getPos().y(), m_HFClusUInTower[icl]->getPos().z(), m_HFClusUInTower[icl]);
    printf(", cousin size %d, address: ", m_HFClusUInTower[icl]->getHalfClusterCol("CousinCluster").size());
    for(int ics=0; ics<m_HFClusUInTower[icl]->getHalfClusterCol("CousinCluster").size(); ics++) printf("%p, ", m_HFClusUInTower[icl]->getHalfClusterCol("CousinCluster")[ics]);
    printf(", track size %d, address: ", m_HFClusUInTower[icl]->getAssociatedTracks().size());
    for(int itrk=0; itrk<m_HFClusUInTower[icl]->getAssociatedTracks().size(); itrk++) printf("%p, ", m_HFClusUInTower[icl]->getAssociatedTracks()[itrk]);
    cout<<endl;
    for(auto ish : m_HFClusUInTower[icl]->getCluster()){
      printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, Address %p \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(),  ish->getNseeds(), ish );
    }
  }
  cout<<endl;
  cout<<"    Loop print HalfClusterV: "<<endl;
  for(int icl=0; icl<m_HFClusVInTower.size(); icl++){
    cout<<"      In HFClusV #"<<icl<<": shower size = "<<m_HFClusVInTower[icl]->getCluster().size()<<", En = "<<m_HFClusVInTower[icl]->getEnergy()<<", type "<<m_HFClusVInTower[icl]->getType();
    printf(", Position (%.3f, %.3f, %.3f), address %p ",m_HFClusVInTower[icl]->getPos().x(), m_HFClusVInTower[icl]->getPos().y(), m_HFClusVInTower[icl]->getPos().z(), m_HFClusVInTower[icl]);
    printf(", cousin size %d, address: ", m_HFClusVInTower[icl]->getHalfClusterCol("CousinCluster").size());
    for(int ics=0; ics<m_HFClusVInTower[icl]->getHalfClusterCol("CousinCluster").size(); ics++) printf("%p, ", m_HFClusVInTower[icl]->getHalfClusterCol("CousinCluster")[ics]);
    printf(", track size %d, address: ", m_HFClusVInTower[icl]->getAssociatedTracks().size());
    for(int itrk=0; itrk<m_HFClusVInTower[icl]->getAssociatedTracks().size(); itrk++) printf("%p, ", m_HFClusVInTower[icl]->getAssociatedTracks()[itrk]);
    cout<<endl;
    for(auto ish : m_HFClusVInTower[icl]->getCluster()){
      printf("          Shower layer %d, Pos+E (%.3f, %.3f, %.3f, %.3f), Nbars %d, NSeed %d, Address %p \n", ish->getDlayer(), ish->getPos().x(), ish->getPos().y(), ish->getPos().z(), ish->getEnergy(), ish->getBars().size(), ish->getNseeds(), ish );
    }
  }
  cout<<endl;
}
*/
  //Remove the same clusters
  for(int ic=0; ic<m_clusterCol.size() && m_clusterCol.size()>1; ic++){
    for(int jc=ic+1; jc<m_clusterCol.size(); jc++){
      if(ic>m_clusterCol.size()) ic--;

      //From pos+E. Can confirm from ChildCluster info.
      float m_clusE_ic = m_clusterCol[ic].get()->getEnergy();
      TVector3 m_pos_ic = m_clusterCol[ic].get()->getShowerCenter();
      float m_clusE_jc = m_clusterCol[jc].get()->getEnergy();
      TVector3 m_pos_jc = m_clusterCol[jc].get()->getShowerCenter();
      if( fabs(m_clusE_ic-m_clusE_jc)<0.1 && (m_pos_ic-m_pos_jc).Mag()<10 ){
        m_clusterCol.erase(m_clusterCol.begin()+jc);
        jc--;
        if(ic>jc+1) ic--;
      }
    }
  }

//cout<<"  ClusterReconnecting: Cluster size after double-counting check: "<<m_clusterCol.size()<<endl;

  //  Type2: from longiCluster aspect
  for(int ic=0; ic<m_clusterCol.size() && m_clusterCol.size()>1; ic++){
    if( m_clusterCol[ic].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster").size()==0 ||
        m_clusterCol[ic].get()->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster").size()==0 ) continue;
//cout<<"  Cluster ic "<<ic<<" has parent"<<endl;
    for(int jc=ic+1; jc<m_clusterCol.size(); jc++){
      if( m_clusterCol[jc].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster").size()==0 ||
          m_clusterCol[jc].get()->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster").size()==0 ) continue;
//cout<<"  Cluster jc "<<jc<<" also has parent, start check"<<endl;
/*
      std::vector<const Cyber::CaloHalfCluster*> m_cousinU = m_clusterCol[ic].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster");
      std::vector<const Cyber::CaloHalfCluster*> m_cousinV = m_clusterCol[ic].get()->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster");
      if( find(m_cousinU.begin(), m_cousinU.end(), m_clusterCol[jc].get()->getHalfClusterUCol("LinkedLongiCluster")[0] )!=m_cousinU.end() &&
          find(m_cousinV.begin(), m_cousinV.end(), m_clusterCol[jc].get()->getHalfClusterVCol("LinkedLongiCluster")[0] )!=m_cousinV.end() &&
          m_clusterCol[ic].get()->getShowerCenter().Angle( m_clusterCol[jc].get()->getShowerCenter() )<0.05 ){ //3deg, ~10cm in ECAL.
        m_clusterCol[ic].get()->mergeCluster( m_clusterCol[jc].get() );
        m_clusterCol.erase(m_clusterCol.begin()+jc);
        jc--;
        if(jc<ic) jc=ic;
      }

      m_cousinU.clear();
      m_cousinV.clear();
      m_cousinU = m_clusterCol[jc].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster");
      m_cousinV = m_clusterCol[jc].get()->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster");
      if( find(m_cousinU.begin(), m_cousinU.end(), m_clusterCol[ic].get()->getHalfClusterUCol("LinkedLongiCluster")[0] )!=m_cousinU.end() &&
          find(m_cousinV.begin(), m_cousinV.end(), m_clusterCol[ic].get()->getHalfClusterVCol("LinkedLongiCluster")[0] )!=m_cousinV.end() &&
          m_clusterCol[ic].get()->getShowerCenter().Angle( m_clusterCol[jc].get()->getShowerCenter() )<0.05 ){ //3deg, ~10cm in ECAL.
        m_clusterCol[ic].get()->mergeCluster( m_clusterCol[jc].get() );
        m_clusterCol.erase(m_clusterCol.begin()+jc);
        jc--;
        if(jc<ic) jc=ic;
      }
*/
      if( m_clusterCol[ic].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster").size()>1 || 
          m_clusterCol[ic].get()->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster").size()>1 ||
          m_clusterCol[jc].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster").size()>1 ||
          m_clusterCol[jc].get()->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster").size()>1 ){

        cout<<"  ERROR: more than 1 Parent Cluster!"<<endl; 
        printf("    In #%d cluster: parent cluster size [%d, %d]. In  #%d cluster: parent cluster size [%d, %d] \n", ic, m_clusterCol[ic].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster").size(), m_clusterCol[ic].get()->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster").size(), jc, m_clusterCol[jc].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster").size(), m_clusterCol[jc].get()->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster").size() );
        continue;
      }      

      const Cyber::CaloHalfCluster* m_parentU1 = m_clusterCol[ic].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster")[0]; 
      const Cyber::CaloHalfCluster* m_parentV1 = m_clusterCol[ic].get()->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster")[0]; 
      const Cyber::CaloHalfCluster* m_parentU2 = m_clusterCol[jc].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster")[0];
      const Cyber::CaloHalfCluster* m_parentV2 = m_clusterCol[jc].get()->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("ParentCluster")[0];
      if(m_parentU1==m_parentU2 && m_parentV1==m_parentV2 && m_clusterCol[ic].get()->getShowerCenter().Angle( m_clusterCol[jc].get()->getShowerCenter() )<0.05){
//cout<<"  In pair ("<<ic<<", "<<jc<<"): merge "<<endl;
        m_clusterCol[ic].get()->mergeCluster( m_clusterCol[jc].get() );
        m_clusterCol.erase(m_clusterCol.begin()+jc);
        jc--;
        if(jc<ic) jc=ic;
      }
  }}

//cout<<"  ClusterReconnecting: Cluster size after neighbor check: "<<m_clusterCol.size()<<endl;

  //  Type3: merge clusters linked to the same track.
  for(int ic=0; ic<m_clusterCol.size() && m_clusterCol.size()>1; ic++){
    if(m_clusterCol[ic].get()->getAssociatedTracks().size()==0) continue;
    std::vector<const Cyber::Track*> m_trkCol = m_clusterCol[ic].get()->getAssociatedTracks();

    for(int jc=ic+1; jc<m_clusterCol.size(); jc++){
      if(m_clusterCol[jc].get()->getAssociatedTracks().size()==0) continue;
//cout<<"Check pair: ["<<ic<<", "<<jc<<"]. Both have tracks. "<<endl;

      for(int itrk=0; itrk<m_clusterCol[jc].get()->getAssociatedTracks().size(); itrk++){
        if( find(m_trkCol.begin(), m_trkCol.end(), m_clusterCol[jc].get()->getAssociatedTracks()[itrk])!= m_trkCol.end() ){
//cout<<"  Merge cluster pair: ["<<ic<<", "<<jc<<"]. "<<endl;
          m_clusterCol[ic].get()->mergeCluster( m_clusterCol[jc].get() );
          m_clusterCol.erase(m_clusterCol.begin()+jc);
          jc--;
          if(jc<ic) jc=ic;
        }
//cout<<"  After merge: ic="<<ic<<", jc="<<jc<<endl;
        break;
      }
    }
  }
  for(int ic=0; ic<m_clusterCol.size(); ic++) m_clusterCol[ic].get()->getLinkedMCPfromHFCluster("LinkedLongiCluster");

//cout<<"  ClusterReconnecting: Cluster size after track check: "<<m_clusterCol.size()<<endl;

/*  //Check the clusters with the same ChildHFCluster: re-matching. 
  //std::multimap<const CaloHalfCluster*, std::shared_ptr<Cyber::Calo3DCluster>> map_linkedClusters; map_linkedClusters.clear();
  std::vector<std::shared_ptr<Cyber::Calo3DCluster>> Cluster_needMatch; Cluster_needMatch.clear();
  for(int ic=0; ic<m_clusterCol.size(); ic++){
    std::vector<const CaloHalfCluster*> m_HFClusU = m_clusterCol[ic].get()->getHalfClusterUCol("LinkedLongiCluster");
    std::vector<const CaloHalfCluster*> m_HFClusV = m_clusterCol[ic].get()->getHalfClusterVCol("LinkedLongiCluster");

    if(m_HFClusU.size()!=1 || m_HFClusV.size()!=1){
      std::cout<<"  ERROR: 3Dcluster has more than 1 LinkedLongiCluster: ("<<m_HFClusU.size()<<", "<<m_HFClusV.size()<<")! Check! "<<std::endl;
      continue;
    }
    std::vector<const CaloHalfCluster*> m_ChildClU = m_HFClusU[0]->getHalfClusterCol("ChildCluster");
    std::vector<const CaloHalfCluster*> m_ChildClV = m_HFClusV[0]->getHalfClusterCol("ChildCluster");
    if(m_ChildClU.size()==0 && m_ChildClV.size()==0) continue; 

    //for(int icl=0; icl<m_ChildClU.size(); icl++) map_linkedClusters.insert( {m_HFClusU[icl], m_clusterCol[ic]} );
    //for(int icl=0; icl<m_ChildClV.size(); icl++) map_linkedClusters.insert( {m_HFClusV[icl], m_clusterCol[ic]} );
    Cluster_needMatch.push_back( m_clusterCol[ic] );
    m_clusterCol.erase( m_clusterCol.begin()+ic );
    ic--;
  }
cout<<"  ClusterReconnecting: Cluster size after child check: "<<m_clusterCol.size()<<endl;
cout<<"  Need-match cluster size: "<<Cluster_needMatch.size()<<endl;

  std::vector<const CaloHalfCluster*> m_ChildClU; 
  std::vector<const CaloHalfCluster*> m_ChildClV; 
  for(int ic=0; ic<Cluster_needMatch.size(); ic++){
    std::vector<const CaloHalfCluster*> tmp_ChildCl = Cluster_needMatch[ic].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("ChildCluster");
    m_ChildClU.insert( m_ChildClU.end(), tmp_ChildCl.begin(), tmp_ChildCl.end() );
    tmp_ChildCl.clear();
    tmp_ChildCl = Cluster_needMatch[ic].get()->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("ChildCluster");
    m_ChildClV.insert( m_ChildClV.end(), tmp_ChildCl.begin(), tmp_ChildCl.end() );
  }
cout<<"  Total child cluster size: ("<<m_ChildClU.size()<<", "<<m_ChildClV.size()<<endl;

  std::vector<const CaloHalfCluster*> m_mergedClusU; m_mergedClusU.clear();
  std::vector<const CaloHalfCluster*> m_eraseCl; m_eraseCl.clear(); 
  do{
    for(int icl=0; icl<m_ChildClU.size(); icl++){
      std::vector<const CaloHalfCluster*> tmp_cousin = m_ChildClU[icl]->getHalfClusterCol("CousinCluster");
      std::shared_ptr<Cyber::CaloHalfCluster> merged_HFClU = std::make_shared<Cyber::CaloHalfCluster>();

      merged_HFClU->mergeHalfCluster(m_ChildClU[icl]);
      merged_HFClU.get()->addHalfCluster("ChildCluster", m_ChildClU[icl]);
      for(auto cl_U : tmp_cousin){
        merged_HFClU.get()->mergeHalfCluster(cl_U);
        merged_HFClU.get()->addHalfCluster("ChildCluster", cl_U);
        m_eraseCl.push_back(cl_U);
      }
      m_eraseCl.push_back(m_ChildClU[icl]);

      m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(merged_HFClU);
      m_mergedClusU.push_back(merged_HFClU.get());
      break;
    }
cout<<"    waiting erasing cluster size: "<<m_eraseCl.size()<<endl;
    for(int icl=0; icl<m_eraseCl.size(); icl++){
      auto iter_find = find( m_ChildClU.begin(), m_ChildClU.end(), m_eraseCl[icl] );
      if(iter_find!=m_ChildClU.end()) m_ChildClU.erase( iter_find );
    }
    m_eraseCl.clear();
cout<<"    left child cluster size: "<<m_ChildClU.size()<<endl;
  }while(m_ChildClU.size()!=0);

  std::vector<const CaloHalfCluster*> m_mergedClusV; m_mergedClusV.clear();
  m_eraseCl.clear();
  do{
    for(int icl=0; icl<m_ChildClV.size(); icl++){
      std::vector<const CaloHalfCluster*> tmp_cousin = m_ChildClV[icl]->getHalfClusterCol("CousinCluster");
      std::shared_ptr<Cyber::CaloHalfCluster> merged_HFClV = std::make_shared<Cyber::CaloHalfCluster>();
      merged_HFClV->mergeHalfCluster(m_ChildClV[icl]);
      merged_HFClV.get()->addHalfCluster("ChildCluster", m_ChildClV[icl]);
      for(auto cl_V : tmp_cousin){
        merged_HFClV.get()->mergeHalfCluster(cl_V);
        merged_HFClV.get()->addHalfCluster("ChildCluster", cl_V);
        m_eraseCl.push_back(cl_V);
      }
      m_eraseCl.push_back(m_ChildClV[icl]);

      m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(merged_HFClV);
      m_mergedClusV.push_back(merged_HFClV.get());
      break;
    }
    for(int icl=0; icl<m_eraseCl.size(); icl++){
      auto iter_find = find( m_ChildClV.begin(), m_ChildClV.end(), m_eraseCl[icl] );
      if(iter_find!=m_ChildClV.end()) m_ChildClV.erase( iter_find );
    }
    m_eraseCl.clear();

  }while(m_ChildClV.size()!=0);
cout<<"  Merged HFCluster size: ("<<m_mergedClusU.size()<<", "<<m_mergedClusV.size()<<endl;

  std::vector<std::shared_ptr<Cyber::Calo3DCluster>> newClusters; newClusters.clear();
  if(m_mergedClusU.size()==0 || m_mergedClusV.size()==0) return StatusCode::SUCCESS;
  else if( m_mergedClusU.size()==1 && m_mergedClusV.size()==1 ){
    std::shared_ptr<Cyber::Calo3DCluster> tmp_clus = std::make_shared<Cyber::Calo3DCluster>();
    XYClusterMatchingL0(m_mergedClusU[0], m_mergedClusV[0], tmp_clus);
    newClusters.push_back(tmp_clus);
  }
  else if( m_mergedClusU.size()==1 ){
    std::vector<std::shared_ptr<Cyber::Calo3DCluster>> emptyCol;
    XYClusterMatchingL1(m_mergedClusU[0], m_mergedClusV, emptyCol);
    newClusters.insert( newClusters.end(), emptyCol.begin(), emptyCol.end() );
  }
  else if( m_mergedClusV.size()==1 ){
    std::vector<std::shared_ptr<Cyber::Calo3DCluster>> emptyCol;
    XYClusterMatchingL1(m_mergedClusV[0], m_mergedClusU, emptyCol);
    newClusters.insert( newClusters.end(), emptyCol.begin(), emptyCol.end() );
  }
  else if( m_mergedClusU.size()==m_mergedClusV.size() ){
    std::vector<std::shared_ptr<Cyber::Calo3DCluster>> emptyCol;
    XYClusterMatchingL2(m_mergedClusU, m_mergedClusV, emptyCol);
    newClusters.insert( newClusters.end(), emptyCol.begin(), emptyCol.end() );
  }
  else{
    std::vector<std::shared_ptr<Cyber::Calo3DCluster>> emptyCol;
    XYClusterMatchingL3(m_mergedClusU, m_mergedClusV, emptyCol);
    newClusters.insert( newClusters.end(), emptyCol.begin(), emptyCol.end() );
  }
cout<<"  Created new 3DCluster size: "<<newClusters.size()<<endl;

  for(int ic=0; ic<newClusters.size(); ic++){
    if( !newClusters[ic].get() ){
      newClusters.erase(newClusters.begin()+ic);
      ic--;
  }}
cout<<"  Created new 3DCluster size: "<<newClusters.size()<<endl;

  m_clusterCol.insert(m_clusterCol.end(), newClusters.begin(), newClusters.end() );
*/
  return StatusCode::SUCCESS;
}



#endif
