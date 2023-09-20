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
  if(settings.map_floatPars.find("th_UVdeltaE")==settings.map_floatPars.end())       settings.map_floatPars["th_UVdeltaE"] = 0.3;
  if(settings.map_floatPars.find("th_UVdeltaE1")==settings.map_floatPars.end())      settings.map_floatPars["th_UVdeltaE1"] = 0.6;
  if(settings.map_floatPars.find("th_ConeTheta")==settings.map_floatPars.end())      settings.map_floatPars["th_ConeTheta"] = TMath::Pi()/6.;
  if(settings.map_floatPars.find("th_ConeR")==settings.map_floatPars.end())          settings.map_floatPars["th_ConeR"] = 30.;
  if(settings.map_boolPars.find("TruthMatch")==settings.map_boolPars.end())          settings.map_boolPars["TruthMatch"] = false;
  if(settings.map_stringPars.find("ReadinHFClusterName")==settings.map_stringPars.end()) settings.map_stringPars["ReadinHFClusterName"] = "ESHalfCluster";
  if(settings.map_stringPars.find("ReadinTowerName")==settings.map_stringPars.end()) settings.map_stringPars["ReadinTowerName"] = "ESTower";
  

  return StatusCode::SUCCESS;
};


StatusCode EnergyTimeMatchingAlg::Initialize( PandoraPlusDataCol& m_datacol ){
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


StatusCode EnergyTimeMatchingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  if(settings.map_boolPars["TruthMatch"]){

    for(int it=0; it<m_towerCol.size(); it++){

      std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> tmp_clusters; tmp_clusters.clear();
      m_HFClusUCol = m_towerCol.at(it)->getHalfClusterUCol(settings.map_stringPars["ReadinHFClusterName"]+"U");
      m_HFClusVCol = m_towerCol.at(it)->getHalfClusterVCol(settings.map_stringPars["ReadinHFClusterName"]+"V");

printf("  In Tower [%d, %d, %d], ", m_towerCol[it]->getTowerID()[0][0], m_towerCol[it]->getTowerID()[0][1], m_towerCol[it]->getTowerID()[0][2] );
printf(" HalfCluster size (%d, %d), input total energy %.3f \n", m_HFClusUCol.size(), m_HFClusVCol.size(), m_towerCol[it]->getEnergy());

      TruthMatching(m_HFClusUCol, m_HFClusVCol, tmp_clusters);
      m_clusterCol.insert(m_clusterCol.end(), tmp_clusters.begin(), tmp_clusters.end());    

cout<<"  Reconstructed cluster energy: "<<endl;
for(int i=0; i<tmp_clusters.size(); i++) printf("  Pos+E: (%.3f, %.3f, %.3f). En %.3f \n", tmp_clusters[i]->getShowerCenter().x(), tmp_clusters[i]->getShowerCenter().y(), tmp_clusters[i]->getShowerCenter().z(), tmp_clusters[i]->getEnergy() );
    }


    //Merge clusters linked to the same MCP.
    for(int ic=0; ic<m_clusterCol.size() && m_clusterCol.size()>1; ic++){
      edm4hep::MCParticle mcp1 = m_clusterCol[ic].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getLinkedMCP()[0].first; 
      for(int jc=ic+1; jc<m_clusterCol.size(); jc++){
        if(ic>m_clusterCol.size()) ic--;
        edm4hep::MCParticle mcp2 = m_clusterCol[jc].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getLinkedMCP()[0].first;
        if(mcp1==mcp2){
          m_clusterCol[ic].get()->mergeCluster( m_clusterCol[jc].get() );
          m_clusterCol.erase(m_clusterCol.begin()+jc);
          jc--;
          if(jc<ic) jc=ic;
        }
      }
    }

cout<<"    After merge: cluster size "<<m_clusterCol.size()<<endl;
   for(int ic=0; ic<m_clusterCol.size(); ic++) m_clusterCol[ic].get()->getLinkedMCPfromHFCluster("LinkedLongiCluster");

cout<<"    After merge: cluster size "<<m_clusterCol.size()<<endl;
cout<<"    Print cluster MC truth info"<<endl;
for(int i=0; i<m_clusterCol.size(); i++){
  printf("  Cluster #%d: En = %.3f, Ntrk %d, Nmcp = %d \n ",i, m_clusterCol[i]->getLongiE(), m_clusterCol[i]->getAssociatedTracks().size(), m_clusterCol[i]->getLinkedMCP().size());
  for(auto imcp : m_clusterCol[i]->getLinkedMCP()){
    printf("    MCP mom (%.3f, %.3f, %.3f), PDG %d, weight %.3f \n", imcp.first.getMomentum().x, imcp.first.getMomentum().y, imcp.first.getMomentum().z, imcp.first.getPDG(), imcp.second);
  }

}

    m_datacol.map_CaloCluster["EcalCluster"] = m_clusterCol;

    //m_datacol.map_CaloHit["bkHit"].insert( m_datacol.map_CaloHit["bkHit"].end(), m_bkCol.map_CaloHit["bkHit"].begin(), m_bkCol.map_CaloHit["bkHit"].end() );
    m_datacol.map_BarCol["bkBar"].insert( m_datacol.map_BarCol["bkBar"].end(), m_bkCol.map_BarCol["bkBar"].begin(), m_bkCol.map_BarCol["bkBar"].end() );
    m_datacol.map_1DCluster["bk1DCluster"].insert( m_datacol.map_1DCluster["bk1DCluster"].end(), m_bkCol.map_1DCluster["bk1DCluster"].begin(), m_bkCol.map_1DCluster["bk1DCluster"].end() );
    m_datacol.map_2DCluster["bk2DCluster"].insert( m_datacol.map_2DCluster["bk2DCluster"].end(), m_bkCol.map_2DCluster["bk2DCluster"].begin(), m_bkCol.map_2DCluster["bk2DCluster"].end() );
    m_datacol.map_HalfCluster["bkHalfCluster"].insert( m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_bkCol.map_HalfCluster["bkHalfCluster"].begin(), m_bkCol.map_HalfCluster["bkHalfCluster"].end() );
    return StatusCode::SUCCESS;
  }

  m_clusterCol.clear();
  std::vector<const PandoraPlus::CaloHalfCluster*> m_leftHFClusterUCol; 
  std::vector<const PandoraPlus::CaloHalfCluster*> m_leftHFClusterVCol; 

  //Loop for towers:  
  for(int it=0; it<m_towerCol.size(); it++){
    m_HFClusUCol.clear(); m_HFClusVCol.clear();
cout<<"Check tower ID: ";
for(int i=0; i<m_towerCol[it]->getTowerID().size(); i++) printf("[%d, %d, %d], ", m_towerCol[it]->getTowerID()[i][0], m_towerCol[it]->getTowerID()[i][1], m_towerCol[it]->getTowerID()[i][2]);
cout<<endl;

    m_HFClusUCol = m_towerCol.at(it)->getHalfClusterUCol(settings.map_stringPars["ReadinHFClusterName"]+"U");
    m_HFClusVCol = m_towerCol.at(it)->getHalfClusterVCol(settings.map_stringPars["ReadinHFClusterName"]+"V");

    std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> tmp_clusters; tmp_clusters.clear();

printf("  In Tower [%d, %d, %d], ", m_towerCol[it]->getTowerID()[0][0], m_towerCol[it]->getTowerID()[0][1], m_towerCol[it]->getTowerID()[0][2] );
printf(" HalfCluster size (%d, %d), input total energy %.3f \n", m_HFClusUCol.size(), m_HFClusVCol.size(), m_towerCol[it]->getEnergy());

    //Check 1: pick out 1*1 from track info
    for(int icl=0; icl<m_HFClusUCol.size(); ++icl){
      std::vector<const PandoraPlus::Track*> m_linkedTrkU = m_HFClusUCol[icl]->getAssociatedTracks();
      if(m_linkedTrkU.size()==0) continue;

      std::vector<const PandoraPlus::CaloHalfCluster*> m_TrkMatchedClV; 
      for(int jcl=0; jcl<m_HFClusVCol.size(); ++jcl){
        if(m_HFClusVCol[jcl]->getAssociatedTracks().size()!=1) continue; //If HFClusterV linked to multiple tracks: leave to HFClusterV loop. 
        if( find(m_linkedTrkU.begin(), m_linkedTrkU.end(), m_HFClusVCol[jcl]->getAssociatedTracks()[0])!=m_linkedTrkU.end() ){
          m_TrkMatchedClV.push_back( m_HFClusVCol[jcl] );

          //Remove HFClusterV from collection
          m_HFClusVCol.erase(m_HFClusVCol.begin()+jcl);
          jcl--;
        }
      }

      if(m_TrkMatchedClV.size()==0) continue;
      else if(m_TrkMatchedClV.size()==1){ //1*1 matching
        std::shared_ptr<PandoraPlus::Calo3DCluster> tmp_clus = std::make_shared<PandoraPlus::Calo3DCluster>();
        XYClusterMatchingL0(m_HFClusUCol[icl], m_TrkMatchedClV[0], tmp_clus);
        tmp_clusters.push_back(tmp_clus);
      }
      else{ // 1*N matching
        std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> emptyCol;
        XYClusterMatchingL1(m_HFClusUCol[icl], m_TrkMatchedClV, emptyCol);
        tmp_clusters.insert(tmp_clusters.end(), emptyCol.begin(), emptyCol.end());
      }
      m_HFClusUCol.erase(m_HFClusUCol.begin()+icl );

    }

    for(int icl=0; icl<m_HFClusVCol.size(); ++icl){
      std::vector<const PandoraPlus::Track*> m_linkedTrkV = m_HFClusVCol[icl]->getAssociatedTracks();
      if(m_linkedTrkV.size()==0) continue;

      std::vector<const PandoraPlus::CaloHalfCluster*> m_TrkMatchedClU;
      for(int jcl=0; jcl<m_HFClusUCol.size(); ++jcl){
        if(m_HFClusUCol[jcl]->getAssociatedTracks().size()!=1) continue; //If HFClusterU linked to multiple tracks: leave to HFClusterV loop.
        if( find(m_linkedTrkV.begin(), m_linkedTrkV.end(), m_HFClusUCol[jcl]->getAssociatedTracks()[0])!=m_linkedTrkV.end() ){
          m_TrkMatchedClU.push_back( m_HFClusUCol[jcl] );

          //Remove HFClusterU from collection
          m_HFClusUCol.erase(m_HFClusUCol.begin()+jcl);
          jcl--;
        }
      }

      if(m_TrkMatchedClU.size()==0) continue;
      else if(m_TrkMatchedClU.size()==1){ //1*1 matching
        std::shared_ptr<PandoraPlus::Calo3DCluster> tmp_clus = std::make_shared<PandoraPlus::Calo3DCluster>();
        XYClusterMatchingL0(m_HFClusVCol[icl], m_TrkMatchedClU[0], tmp_clus);
        tmp_clusters.push_back(tmp_clus);
      }
      else{ // 1*N matching
        std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> emptyCol;
        XYClusterMatchingL1(m_HFClusVCol[icl], m_TrkMatchedClU, emptyCol);
        tmp_clusters.insert(tmp_clusters.end(), emptyCol.begin(), emptyCol.end());
      }
      m_HFClusVCol.erase(m_HFClusVCol.begin()+icl );

    }

    //In case: actually no solution. 
    //TODO: write a part for this?    

printf("  After Track check: HFCluster size (%d, %d) \n", m_HFClusUCol.size(), m_HFClusVCol.size() );

    //Check 2: pick out 1*1 from cousin info
    for(int icl=0; icl<m_HFClusUCol.size(); ++icl){
      if(m_HFClusUCol[icl]->getHalfClusterCol("CousinCluster").size()==0) continue;
      for(int jcl=0; jcl<m_HFClusVCol.size(); ++jcl){
        if(m_HFClusVCol[jcl]->getHalfClusterCol("CousinCluster").size()==0) continue;

        std::vector<const PandoraPlus::CaloHalfCluster*> m_cousinU = m_HFClusUCol[icl]->getHalfClusterCol("CousinCluster");
        std::vector<const PandoraPlus::CaloHalfCluster*> m_cousinV = m_HFClusVCol[jcl]->getHalfClusterCol("CousinCluster");

      //Option 1: Check deltaE in each tower. Link the HFClusters only when deltaE in each tower < threshold. 
        //Make map: <towerID, (E_clU - E_clV)/mean(E_clU, E_clV)  >
        std::map<std::vector<int>, pair<float, float>> map_pairE; map_pairE.clear();
        //map_pairE[m_towerCol[it]->getTowerID()[0]] = 2.*fabs(m_HFClusUCol[icl]->getEnergy()-m_HFClusVCol[jcl]->getEnergy())/(m_HFClusUCol[icl]->getEnergy()+m_HFClusVCol[jcl]->getEnergy());
        map_pairE[m_towerCol[it]->getTowerID()[0]] = make_pair(m_HFClusUCol[icl]->getEnergy(), m_HFClusVCol[jcl]->getEnergy());
        for(int ics=0; ics<m_cousinU.size(); ++ics){
        for(int jcs=0; jcs<m_cousinV.size(); ++jcs){
          if(m_cousinU[ics]->getTowerID()[0] == m_cousinV[jcs]->getTowerID()[0])
            //map_pairE[ m_cousinU[ics]->getTowerID()[0] ] = 2.*fabs(m_cousinU[ics]->getEnergy()-m_cousinV[jcs]->getEnergy())/(m_cousinU[ics]->getEnergy()+m_cousinV[jcs]->getEnergy());
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
          if(iter.second.first>1 && iter.second.second>1 && deltaE > settings.map_floatPars["th_UVdeltaE1"]){ isLink = false; break; }
          else if(iter.second.first/iter.second.second > 10 || iter.second.first/iter.second.second<0.1) { isLink = false; break; }
        }
        if( 2*fabs(totE_U-totE_V)/(totE_U+totE_V)<settings.map_floatPars["th_UVdeltaE"] ) isLink=true;

printf("  In pair (%d, %d): totE (%.3f, %.3f), deltaE = %.3f, isLink = %d \n", icl, jcl, totE_U, totE_V, 2*fabs(totE_U-totE_V)/(totE_U+totE_V), isLink);
        if(isLink){

          std::vector<const PandoraPlus::CaloHalfCluster*> m_parentHFClU = m_HFClusUCol[icl]->getHalfClusterCol("ParentCluster");
          std::vector<const PandoraPlus::CaloHalfCluster*> m_parentHFClV = m_HFClusVCol[jcl]->getHalfClusterCol("ParentCluster");

          if(m_parentHFClU.size()==0) m_parentHFClU.push_back(m_HFClusUCol[icl]);
          if(m_parentHFClV.size()==0) m_parentHFClV.push_back(m_HFClusVCol[jcl]);
          if(m_parentHFClU.size()!=1 || m_parentHFClV.size()!=1) cout<<"ERROR: more than 1 parent HFCluster! Check this. "<<endl;

          std::shared_ptr<PandoraPlus::Calo3DCluster> tmp_clus = std::make_shared<PandoraPlus::Calo3DCluster>();
          XYClusterMatchingL0(m_parentHFClU[0], m_parentHFClV[0], tmp_clus);
          tmp_clusters.push_back(tmp_clus);

          m_HFClusUCol.erase(m_HFClusUCol.begin()+icl );
          m_HFClusVCol.erase(m_HFClusVCol.begin()+jcl );
          icl--; jcl--;
          break;

        }
      }
    }
printf("  After cousin check: HFCluster size (%d, %d) \n", m_HFClusUCol.size(), m_HFClusVCol.size() );

    //Check 3: pick out 1*N from totEn and longitudinal profile info
    std::vector<const PandoraPlus::CaloHalfCluster*> m_HFnoCousinU; m_HFnoCousinU.clear();
    std::vector<const PandoraPlus::CaloHalfCluster*> m_HFnoCousinV; m_HFnoCousinV.clear();
    float totE_U_noCousin = 0.;
    float totE_V_noCousin = 0.;
    for(int icl=0; icl<m_HFClusUCol.size(); ++icl){
      if(m_HFClusUCol[icl]->getHalfClusterCol("CousinCluster").size()!=0) continue;
      m_HFnoCousinU.push_back(m_HFClusUCol[icl]);
      totE_U_noCousin += m_HFClusUCol[icl]->getEnergy();
    }
    for(int jcl=0; jcl<m_HFClusVCol.size(); ++jcl){
      if(m_HFClusVCol[jcl]->getHalfClusterCol("CousinCluster").size()!=0) continue;

      m_HFnoCousinV.push_back(m_HFClusVCol[jcl]);
      totE_V_noCousin += m_HFClusVCol[jcl]->getEnergy();
    }

    float deltaE_noCousin = fabs(2*(totE_U_noCousin-totE_V_noCousin)/(totE_U_noCousin+totE_V_noCousin));
printf("  No cousin cluster size: [%d, %d], totE: [%.3f, %.3f], deltaE = %.3f \n", m_HFnoCousinU.size(), m_HFnoCousinV.size(), totE_U_noCousin, totE_V_noCousin, deltaE_noCousin);
    if(deltaE_noCousin<settings.map_floatPars["th_UVdeltaE"]){
      if(m_HFnoCousinU.size()==1 && m_HFnoCousinV.size()==1 ){
        std::shared_ptr<PandoraPlus::Calo3DCluster> tmp_clus = std::make_shared<PandoraPlus::Calo3DCluster>();
        XYClusterMatchingL0(m_HFnoCousinU[0], m_HFnoCousinV[0], tmp_clus);
        tmp_clusters.push_back(tmp_clus);        

        auto iter = find(m_HFClusUCol.begin(), m_HFClusUCol.end(), m_HFnoCousinU[0] );
        if(iter!=m_HFClusUCol.end()) m_HFClusUCol.erase(iter);
        iter = find(m_HFClusVCol.begin(), m_HFClusVCol.end(), m_HFnoCousinV[0]);
        if(iter!=m_HFClusVCol.end()) m_HFClusVCol.erase(iter);
      }
      else if(m_HFnoCousinU.size()==1){
        std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> emptyCol;
        XYClusterMatchingL1(m_HFClusUCol[0], m_HFClusVCol, emptyCol);
        tmp_clusters.insert(tmp_clusters.end(), emptyCol.begin(), emptyCol.end());

        //Erase HFClusters from origin collection
        auto iter = find(m_HFClusUCol.begin(), m_HFClusUCol.end(), m_HFnoCousinU[0] );
        if(iter!=m_HFClusUCol.end()) m_HFClusUCol.erase(iter);
        for(int iclv=0; iclv<m_HFnoCousinV.size(); iclv++){
          iter = find(m_HFClusVCol.begin(), m_HFClusVCol.end(), m_HFnoCousinV[iclv]);
          if(iter!=m_HFClusVCol.end()) m_HFClusVCol.erase(iter);
        }
      }
      else if(m_HFnoCousinV.size()==1){
        std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> emptyCol;
        XYClusterMatchingL1(m_HFClusVCol[0], m_HFClusUCol, emptyCol);
        tmp_clusters.insert(tmp_clusters.end(), emptyCol.begin(), emptyCol.end());

        //Erase HFClusters from origin collection
        auto iter = find(m_HFClusVCol.begin(), m_HFClusVCol.end(), m_HFnoCousinV[0] );
        if(iter!=m_HFClusVCol.end()) m_HFClusVCol.erase(iter);
        for(int iclu=0; iclu<m_HFnoCousinU.size(); iclu++){
          iter = find(m_HFClusUCol.begin(), m_HFClusUCol.end(), m_HFnoCousinU[iclu]);
          if(iter!=m_HFClusUCol.end()) m_HFClusUCol.erase(iter);
        }
      }
      else if(m_HFnoCousinU.size()==m_HFnoCousinV.size()){ 
        std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> emptyCol;
        XYClusterMatchingL2(m_HFnoCousinU, m_HFnoCousinV, emptyCol);
        tmp_clusters.insert(tmp_clusters.end(), emptyCol.begin(), emptyCol.end());

        for(int iclu=0; iclu<m_HFnoCousinU.size(); iclu++){
          auto iter = find(m_HFClusUCol.begin(), m_HFClusUCol.end(), m_HFnoCousinU[iclu]);
          if(iter!=m_HFClusUCol.end()) m_HFClusUCol.erase(iter);
        }
        for(int iclv=0; iclv<m_HFnoCousinV.size(); iclv++){
          auto iter = find(m_HFClusVCol.begin(), m_HFClusVCol.end(), m_HFnoCousinV[iclv]);
          if(iter!=m_HFClusVCol.end()) m_HFClusVCol.erase(iter);
        }
      }
      else{
        std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> emptyCol;
        XYClusterMatchingL3(m_HFnoCousinU, m_HFnoCousinV, emptyCol);
        tmp_clusters.insert(tmp_clusters.end(), emptyCol.begin(), emptyCol.end());

        for(int iclu=0; iclu<m_HFnoCousinU.size(); iclu++){
          auto iter = find(m_HFClusUCol.begin(), m_HFClusUCol.end(), m_HFnoCousinU[iclu]);
          if(iter!=m_HFClusUCol.end()) m_HFClusUCol.erase(iter);
        }
        for(int iclv=0; iclv<m_HFnoCousinV.size(); iclv++){
          auto iter = find(m_HFClusVCol.begin(), m_HFClusVCol.end(), m_HFnoCousinV[iclv]);
          if(iter!=m_HFClusVCol.end()) m_HFClusVCol.erase(iter);
        }
      }

    }
printf("  After deltaE check: HFCluster size (%d, %d) \n", m_HFClusUCol.size(), m_HFClusVCol.size() );
 

    //Left HFClusters: 
    //Option 1: match in tower. 
    std::vector<const PandoraPlus::CaloHalfCluster*> m_HFClusUColParent; m_HFClusUColParent.clear();
    std::vector<const PandoraPlus::CaloHalfCluster*> m_HFClusVColParent; m_HFClusVColParent.clear();

    //Option 1.1: Use parent HFClusters
    //Option 1.2: Use left HFClusters.     
    for(int icl=0; icl<m_HFClusUCol.size(); icl++){
      std::vector<const PandoraPlus::CaloHalfCluster*> tmp_parent = m_HFClusUCol[icl]->getHalfClusterCol("ParentCluster");
      if(tmp_parent.size()==0)
        m_HFClusUColParent.push_back(m_HFClusUCol[icl]);
      else
        m_HFClusUColParent.insert(m_HFClusUColParent.end(), tmp_parent.begin(), tmp_parent.end());
    }
    for(int icl=0; icl<m_HFClusVCol.size(); icl++){
      std::vector<const PandoraPlus::CaloHalfCluster*> tmp_parent = m_HFClusVCol[icl]->getHalfClusterCol("ParentCluster");
      if(tmp_parent.size()==0)
        m_HFClusVColParent.push_back(m_HFClusVCol[icl]);
      else
        m_HFClusVColParent.insert(m_HFClusVColParent.end(), tmp_parent.begin(), tmp_parent.end());
    }

    //Use original HFCluster or parent HFClusters: depend on UV deltaE. 
    float totE_U = 0; 
    float totE_V = 0;
    for(auto icl: m_HFClusUColParent) totE_U += icl->getEnergy();
    for(auto icl: m_HFClusVColParent) totE_V += icl->getEnergy();
    bool isLink_parent = 2*fabs(totE_U-totE_V)/(totE_U+totE_V)<settings.map_floatPars["th_UVdeltaE"];

    totE_U=0; totE_V=0; 
    for(auto icl: m_HFClusUCol) totE_U += icl->getEnergy();
    for(auto icl: m_HFClusVCol) totE_V += icl->getEnergy();
    bool isLink_inTower = 2*fabs(totE_U-totE_V)/(totE_U+totE_V)<settings.map_floatPars["th_UVdeltaE"];

printf("  In Link: inTower = %d, parent = %d \n", isLink_inTower, isLink_parent);
    //Only when (isLink_inTower=False && isLink_parent=True) we use parent for matching. 
    if( isLink_inTower || !isLink_parent){    
      m_HFClusUColParent = m_HFClusUCol;
      m_HFClusVColParent = m_HFClusVCol;
    }

   
    if(m_HFClusUColParent.size()==1 && m_HFClusVColParent.size()==1){
      std::shared_ptr<PandoraPlus::Calo3DCluster> tmp_clus = std::make_shared<PandoraPlus::Calo3DCluster>();
      XYClusterMatchingL0(m_HFClusUColParent[0], m_HFClusVColParent[0], tmp_clus);
      tmp_clusters.push_back(tmp_clus);

      m_HFClusUColParent.erase(m_HFClusUColParent.begin() );
      m_HFClusVColParent.erase(m_HFClusVColParent.begin() );
    }
    else if( m_HFClusUColParent.size()==1 && m_HFClusVColParent.size()!=0 ){ 
      std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> emptyCol; 
      XYClusterMatchingL1(m_HFClusUColParent[0], m_HFClusVColParent, emptyCol);
      tmp_clusters.insert(tmp_clusters.end(), emptyCol.begin(), emptyCol.end()); 
    }
    else if( m_HFClusUColParent.size()!=0 && m_HFClusVColParent.size()==1 ){ 
      std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> emptyCol; 
      XYClusterMatchingL1(m_HFClusVColParent[0], m_HFClusUColParent, emptyCol); 
      tmp_clusters.insert(tmp_clusters.end(), emptyCol.begin(), emptyCol.end()); 
    }

    //Left cases: N*0, M*N.  
    else if(m_HFClusUColParent.size()==m_HFClusVColParent.size()){
      std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> emptyCol;
      XYClusterMatchingL2(m_HFClusUColParent, m_HFClusVColParent, emptyCol);
      tmp_clusters.insert(tmp_clusters.end(), emptyCol.begin(), emptyCol.end());
    }
    else{
      std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> emptyCol;
      XYClusterMatchingL3(m_HFClusUColParent, m_HFClusVColParent, emptyCol);
      tmp_clusters.insert(tmp_clusters.end(), emptyCol.begin(), emptyCol.end());
      //cout<<"  Left HFCluster pair: ("<<m_HFClusUColParent.size()<<", "<<m_HFClusVColParent.size()<<"), match to 1 3DCluster. "<<endl;
      //TODO: add an algorithm here. 
    }

/*
    //Option 2: save out the parents, doing match outside of tower.
    std::vector<const PandoraPlus::CaloHalfCluster*> m_HFClusUColParent; m_HFClusUColParent.clear();
    std::vector<const PandoraPlus::CaloHalfCluster*> m_HFClusVColParent; m_HFClusVColParent.clear();
    for(int icl=0; icl<m_HFClusUCol.size(); icl++){
      std::vector<const PandoraPlus::CaloHalfCluster*> tmp_parent = m_HFClusUCol[icl]->getHalfClusterCol("ParentCluster");
      if(tmp_parent.size()==0)
        m_HFClusUColParent.push_back(m_HFClusUCol[icl]);
      else
        m_HFClusUColParent.insert(m_HFClusUColParent.end(), tmp_parent.begin(), tmp_parent.end());
    }
    for(int icl=0; icl<m_HFClusVCol.size(); icl++){
      std::vector<const PandoraPlus::CaloHalfCluster*> tmp_parent = m_HFClusVCol[icl]->getHalfClusterCol("ParentCluster");
      if(tmp_parent.size()==0)
        m_HFClusVColParent.push_back(m_HFClusVCol[icl]);
      else
        m_HFClusVColParent.insert(m_HFClusVColParent.end(), tmp_parent.begin(), tmp_parent.end());
    }


    for(int icl=0; icl<m_HFClusUColParent.size(); icl++){
      if(find(m_leftHFClusterUCol.begin(), m_leftHFClusterUCol.end(), m_HFClusUColParent[icl])==m_leftHFClusterUCol.end())
        m_leftHFClusterUCol.push_back(m_HFClusUColParent[icl]);
    }
    for(int icl=0; icl<m_HFClusVColParent.size(); icl++){
      if(find(m_leftHFClusterVCol.begin(), m_leftHFClusterVCol.end(), m_HFClusVColParent[icl])==m_leftHFClusterVCol.end())
        m_leftHFClusterVCol.push_back(m_HFClusVColParent[icl]);
    }
*/

    //Clean empty Calo3DClusters
    for(int ic=0; ic<tmp_clusters.size(); ic++){
      if( !tmp_clusters[ic].get() ){
        tmp_clusters.erase(tmp_clusters.begin()+ic);
        ic--;
    }}

    //Save Calo3DClusters and Calo2DClusters into dataCol backupCol.
cout<<"  Reconstructed cluster energy: "<<endl;
for(int i=0; i<tmp_clusters.size(); i++) printf("  Pos+E: (%.3f, %.3f, %.3f). En %.3f \n", tmp_clusters[i]->getShowerCenter().x(), tmp_clusters[i]->getShowerCenter().y(), tmp_clusters[i]->getShowerCenter().z(), tmp_clusters[i]->getEnergy() );

    m_clusterCol.insert( m_clusterCol.end(), tmp_clusters.begin(), tmp_clusters.end() );
    tmp_clusters.clear();

    //for(int ic=0; ic<tmp_clusters.size(); ic++){
    //  m_clusterCol.push_back( tmp_clusters[ic] );
    //  std::vector<const PandoraPlus::Calo2DCluster*> m_showersinclus = tmp_clusters[ic]->getCluster();
    //  for(int is=0; is<m_showersinclus.size(); is++){
    //    m_transhowerCol.push_back( const_cast<PandoraPlus::Calo2DCluster *>(m_showersinclus[is]) );
    //  }
    //}   

  }//End loop towers

//Match for the left clusters:
  

cout<<"  Print reconstructed cluster energy and trk: "<<endl;
//for(int i=0; i<m_clusterCol.size(); i++) cout<<"En = "<<m_clusterCol[i].get()->getEnergy()<<", trk size "<<m_clusterCol[i].get()->getAssociatedTracks().size()<<endl;
//for(int ic=0; ic<m_clusterCol.size(); ic++){
//  cout<<"    Cluster #"<<ic<<": En = "<<m_clusterCol[ic]->getEnergy()<<", Nhit = "<<m_clusterCol[ic]->getCluster().size()<<", Ntrk = "<<m_clusterCol[ic].get()->getAssociatedTracks().size();
//  printf(", pos (%.2f, %.2f, %.2f) \n", m_clusterCol[ic]->getShowerCenter().x(), m_clusterCol[ic]->getShowerCenter().y(), m_clusterCol[ic]->getShowerCenter().z());
//
//  cout<<"    Print hits: "<<endl;
//  for(int ihit=0; ihit<m_clusterCol[ic]->getCluster().size(); ihit++){
//    printf("      Pos+E (%.3f, %.3f, %.3f, %.3f), address %p \n", m_clusterCol[ic]->getCluster()[ihit]->getPos().x(), m_clusterCol[ic]->getCluster()[ihit]->getPos().y(), m_clusterCol[ic]->getCluster()[ihit]->getPos().z(), m_clusterCol[ic]->getCluster()[ihit]->getEnergy(), m_clusterCol[ic]->getCluster()[ihit]);
//  }
//  cout<<endl;
//}



  //Clean the same clusters
  ClusterReconnecting( m_clusterCol );

/*
cout<<"Cluster Merge Type2: from HalfCluster aspect"<<endl;
  //  Type2: from longiCluster aspect
  for(int ic=0; ic<m_clusterCol.size() && m_clusterCol.size()>1; ic++){
    if( m_clusterCol[ic].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster").size()==0 && 
        m_clusterCol[ic].get()->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster").size()==0 ) continue;
//cout<<"  Cluster ic has cousin"<<endl;
    for(int jc=ic+1; jc<m_clusterCol.size(); jc++){
      if( m_clusterCol[jc].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster").size()==0 && 
          m_clusterCol[jc].get()->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster").size()==0 ) continue;
//cout<<"  Cluster jc also has cousin, start check"<<endl;

      std::vector<const PandoraPlus::CaloHalfCluster*> m_cousinU = m_clusterCol[ic].get()->getHalfClusterUCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster"); 
      std::vector<const PandoraPlus::CaloHalfCluster*> m_cousinV = m_clusterCol[ic].get()->getHalfClusterVCol("LinkedLongiCluster")[0]->getHalfClusterCol("CousinCluster"); 
      if( find(m_cousinU.begin(), m_cousinU.end(), m_clusterCol[jc].get()->getHalfClusterUCol("LinkedLongiCluster")[0] )!=m_cousinU.end() &&
          find(m_cousinV.begin(), m_cousinV.end(), m_clusterCol[jc].get()->getHalfClusterVCol("LinkedLongiCluster")[0] )!=m_cousinV.end() && 
          m_clusterCol[ic].get()->getShowerCenter().Angle( m_clusterCol[jc].get()->getShowerCenter() )<0.05 ){ //3deg, ~10cm in ECAL. 
        m_clusterCol[ic].get()->mergeCluster( m_clusterCol[jc].get() );
        m_clusterCol.erase(m_clusterCol.begin()+jc);
        jc--;
        if(jc<ic) jc=ic;
      }        
  }}
cout<<"  After Type2 merge: size = "<<m_clusterCol.size()<<endl;
cout<<"  Print reconstructed cluster energy and trk: "<<endl;
for(int i=0; i<m_clusterCol.size(); i++) cout<<"En = "<<m_clusterCol[i].get()->getEnergy()<<", trk size "<<m_clusterCol[i].get()->getAssociatedTracks().size()<<endl;
*/


//cout<<"Cluster Merge Tyep3: merge clusters linked to the same track. Cluster size: "<<m_clusterCol.size()<<endl;
  //  Type3: merge clusters linked to the same track. 
  for(int ic=0; ic<m_clusterCol.size() && m_clusterCol.size()>1; ic++){
    if(m_clusterCol[ic].get()->getAssociatedTracks().size()==0) continue;
    std::vector<const PandoraPlus::Track*> m_trkCol = m_clusterCol[ic].get()->getAssociatedTracks();

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

cout<<"    After merge: cluster size "<<m_clusterCol.size()<<endl;
cout<<"    Print cluster MC truth info"<<endl;
for(int i=0; i<m_clusterCol.size(); i++){
  printf("  Cluster #%d: En = %.3f, Ntrk %d, Nmcp = %d \n ",i, m_clusterCol[i]->getLongiE(), m_clusterCol[i]->getAssociatedTracks().size(), m_clusterCol[i]->getLinkedMCP().size());
  for(auto imcp : m_clusterCol[i]->getLinkedMCP()){
    printf("    MCP mom (%.3f, %.3f, %.3f), PDG %d, weight %.3f \n", imcp.first.getMomentum().x, imcp.first.getMomentum().y, imcp.first.getMomentum().z, imcp.first.getPDG(), imcp.second);
  }

}

  m_datacol.map_CaloCluster["EcalCluster"] = m_clusterCol;

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

StatusCode EnergyTimeMatchingAlg::ClusterReconnecting( std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>& m_clusterCol ){

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

cout<<"  ClusterReconnecting: Cluster size after double-counting check: "<<m_clusterCol.size()<<endl;

/*  //Check the clusters with the same ChildHFCluster: re-matching. 
  //std::multimap<const CaloHalfCluster*, std::shared_ptr<PandoraPlus::Calo3DCluster>> map_linkedClusters; map_linkedClusters.clear();
  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> Cluster_needMatch; Cluster_needMatch.clear();
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
      std::shared_ptr<PandoraPlus::CaloHalfCluster> merged_HFClU = std::make_shared<PandoraPlus::CaloHalfCluster>();

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
      std::shared_ptr<PandoraPlus::CaloHalfCluster> merged_HFClV = std::make_shared<PandoraPlus::CaloHalfCluster>();
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

  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> newClusters; newClusters.clear();
  if(m_mergedClusU.size()==0 || m_mergedClusV.size()==0) return StatusCode::SUCCESS;
  else if( m_mergedClusU.size()==1 && m_mergedClusV.size()==1 ){
    std::shared_ptr<PandoraPlus::Calo3DCluster> tmp_clus = std::make_shared<PandoraPlus::Calo3DCluster>();
    XYClusterMatchingL0(m_mergedClusU[0], m_mergedClusV[0], tmp_clus);
    newClusters.push_back(tmp_clus);
  }
  else if( m_mergedClusU.size()==1 ){
    std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> emptyCol;
    XYClusterMatchingL1(m_mergedClusU[0], m_mergedClusV, emptyCol);
    newClusters.insert( newClusters.end(), emptyCol.begin(), emptyCol.end() );
  }
  else if( m_mergedClusV.size()==1 ){
    std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> emptyCol;
    XYClusterMatchingL1(m_mergedClusV[0], m_mergedClusU, emptyCol);
    newClusters.insert( newClusters.end(), emptyCol.begin(), emptyCol.end() );
  }
  else if( m_mergedClusU.size()==m_mergedClusV.size() ){
    std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> emptyCol;
    XYClusterMatchingL2(m_mergedClusU, m_mergedClusV, emptyCol);
    newClusters.insert( newClusters.end(), emptyCol.begin(), emptyCol.end() );
  }
  else{
    std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> emptyCol;
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


StatusCode EnergyTimeMatchingAlg::TruthMatching( std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClUCol,
                                                 std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClVCol,
                                                 std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>& m_clusters )
{
  if(m_ClUCol.size()==0 || m_ClVCol.size()==0) return StatusCode::SUCCESS;
cout<<"TruthMatching: readin cluster size: "<<m_ClUCol.size()<<", "<<m_ClVCol.size()<<endl;
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_truthESClUCol; m_truthESClUCol.clear();
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_truthESClVCol; m_truthESClVCol.clear();


  //Truth split in HFClusterU
  for(int icl=0; icl<m_ClUCol.size(); icl++){
    auto truthMap = const_cast<CaloHalfCluster*>(m_ClUCol[icl])->getLinkedMCPfromUnit();
printf("  In HFClusterU #%d: En %.3f, truthMC map size %d \n", icl, m_ClUCol[icl]->getEnergy(), truthMap.size());
int counter = 0;
for(auto imcp : truthMap){
  printf("      MCParticle #%d: truth En %.3f, mom (%.3f, %.3f, %.3f), weight %.3f \n", counter, imcp.first.getEnergy(), imcp.first.getMomentum().x, imcp.first.getMomentum().y, imcp.first.getMomentum().z, imcp.second );
  counter++;
}
cout<<endl;

    for(auto iter: truthMap ){
      if(iter.second<0.05) continue;

      std::shared_ptr<PandoraPlus::CaloHalfCluster> newClus = std::make_shared<PandoraPlus::CaloHalfCluster>();
      //Copy and reweight 1D showers in HFCluster
      for(int ish=0; ish<m_ClUCol[icl]->getCluster().size(); ish++){
        const PandoraPlus::Calo1DCluster* p_shower = m_ClUCol[icl]->getCluster()[ish];

        std::vector<const PandoraPlus::CaloUnit*> Bars; Bars.clear();
        for(int ibar=0; ibar<p_shower->getCluster().size(); ibar++){
          auto bar = p_shower->getCluster()[ibar]->Clone();
          bar->setQ(bar->getQ1()*iter.second, bar->getQ2()*iter.second );

          std::vector<std::pair<edm4hep::MCParticle, float>> emptyMap; emptyMap.clear();
          emptyMap.push_back( make_pair(iter.first, 1.) );
          bar->setLinkedMCP(emptyMap);
          Bars.push_back(bar.get());
          m_bkCol.map_BarCol["bkBar"].push_back(bar);
        }

        std::shared_ptr<PandoraPlus::Calo1DCluster> shower = std::make_shared<PandoraPlus::Calo1DCluster>();
        shower->setBars(Bars);
        shower->setSeed();
        shower->setIDInfo();
        
        newClus->addUnit(shower.get());
        m_bkCol.map_1DCluster["bk1DCluster"].push_back( shower );
      }
      for(int itrk=0; itrk<m_ClUCol[icl]->getAssociatedTracks().size(); itrk++)  newClus->addAssociatedTrack( m_ClUCol[icl]->getAssociatedTracks()[itrk] );
      newClus->setHoughPars( m_ClUCol[icl]->getHoughAlpha(), m_ClUCol[icl]->getHoughRho() );
      newClus->setIntercept( m_ClUCol[icl]->getHoughIntercept() );
      newClus->fitAxis("");
      newClus->getLinkedMCPfromUnit();
      m_truthESClUCol.push_back(newClus);
      //m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(newClus);
    }
  }

  //Merge HFClusters linked to the same MCP
  for(int icl=0; icl<m_truthESClUCol.size() && m_truthESClUCol.size()>1; icl++){
    for(int jcl=icl+1; jcl<m_truthESClUCol.size(); jcl++){
      if(icl>m_truthESClUCol.size()) icl--;

      if(m_truthESClUCol[icl]->getLinkedMCP()[0].first==m_truthESClUCol[jcl]->getLinkedMCP()[0].first){
        m_truthESClUCol[icl].get()->mergeHalfCluster(m_truthESClUCol[jcl].get());
        m_truthESClUCol.erase(m_truthESClUCol.begin()+jcl);
        jcl--;
        if(jcl<icl) jcl=icl;
      }
    }
  }
  for(int icl=0; icl<m_truthESClUCol.size(); icl++) m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(m_truthESClUCol[icl]);

cout<<"  New ES clusterU size: "<<m_truthESClUCol.size()<<endl;
for(int icl=0; icl<m_truthESClUCol.size(); icl++){
  auto truthMap = m_truthESClUCol[icl]->getLinkedMCP();
  printf("  In HFClusterU #%d: En %.3f, truthMC map size %d \n", icl, m_truthESClUCol[icl]->getEnergy(), truthMap.size());
  int counter = 0;
  for(auto imcp : truthMap){
    printf("      MCParticle #%d: truth En %.3f, mom (%.3f, %.3f, %.3f), weight %.3f \n", counter, imcp.first.getEnergy(), imcp.first.getMomentum().x, imcp.first.getMomentum().y, imcp.first.getMomentum().z, imcp.second );
    counter++;
  }
}
cout<<endl;


  //Truth split in HFClusterV
  for(int icl=0; icl<m_ClVCol.size(); icl++){
    auto truthMap = const_cast<CaloHalfCluster*>(m_ClVCol[icl])->getLinkedMCPfromUnit();
printf("  In HFClusterV #%d: En %.3f, truthMC map size %d \n", icl, m_ClVCol[icl]->getEnergy(), truthMap.size());
int counter = 0;
for(auto imcp : truthMap){
  printf("      MCParticle #%d: truth En %.3f, mom (%.3f, %.3f, %.3f), weight %.3f \n", counter, imcp.first.getEnergy(), imcp.first.getMomentum().x, imcp.first.getMomentum().y, imcp.first.getMomentum().z, imcp.second );
  counter++;
}
cout<<endl;

    for(auto iter: truthMap ){
      if(iter.second<0.05) continue;

      std::shared_ptr<PandoraPlus::CaloHalfCluster> newClus = std::make_shared<PandoraPlus::CaloHalfCluster>();
      //Copy and reweight 1D showers in HFCluster
      for(int ish=0; ish<m_ClVCol[icl]->getCluster().size(); ish++){
        const PandoraPlus::Calo1DCluster* p_shower = m_ClVCol[icl]->getCluster()[ish];

        std::vector<const PandoraPlus::CaloUnit*> Bars; Bars.clear();
        for(int ibar=0; ibar<p_shower->getCluster().size(); ibar++){
          auto bar = p_shower->getCluster()[ibar]->Clone();
          bar->setQ(bar->getQ1()*iter.second, bar->getQ2()*iter.second );

          std::vector<std::pair<edm4hep::MCParticle, float>> emptyMap; emptyMap.clear();
          emptyMap.push_back( make_pair(iter.first, 1.) );
          bar->setLinkedMCP(emptyMap);
          Bars.push_back(bar.get());
          m_bkCol.map_BarCol["bkBar"].push_back(bar);
        }

        std::shared_ptr<PandoraPlus::Calo1DCluster> shower = std::make_shared<PandoraPlus::Calo1DCluster>();
        shower->setBars(Bars);
        shower->setSeed();
        shower->setIDInfo();

        newClus->addUnit(shower.get());
        m_bkCol.map_1DCluster["bk1DCluster"].push_back( shower );
      }
      for(int itrk=0; itrk<m_ClVCol[icl]->getAssociatedTracks().size(); itrk++)  newClus->addAssociatedTrack( m_ClVCol[icl]->getAssociatedTracks()[itrk] );
      newClus->setHoughPars( m_ClVCol[icl]->getHoughAlpha(), m_ClVCol[icl]->getHoughRho() );
      newClus->setIntercept( m_ClVCol[icl]->getHoughIntercept() );
      newClus->fitAxis("");
      newClus->getLinkedMCPfromUnit();
      m_truthESClVCol.push_back(newClus);
      //m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(newClus);
    }
  }

  //Merge HFClusters linked to the same MCP
  for(int icl=0; icl<m_truthESClVCol.size() && m_truthESClVCol.size()>1; icl++){
    for(int jcl=icl+1; jcl<m_truthESClVCol.size(); jcl++){
      if(icl>m_truthESClVCol.size()) icl--;

      if(m_truthESClVCol[icl]->getLinkedMCP()[0].first==m_truthESClVCol[jcl]->getLinkedMCP()[0].first){
        m_truthESClVCol[icl].get()->mergeHalfCluster(m_truthESClVCol[jcl].get());
        m_truthESClVCol.erase(m_truthESClVCol.begin()+jcl);
        jcl--;
        if(jcl<icl) jcl=icl;
      }
    }
  }
  for(int icl=0; icl<m_truthESClVCol.size(); icl++) m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(m_truthESClVCol[icl]);

cout<<"  New ES clusterV size: "<<m_truthESClVCol.size()<<endl;
for(int icl=0; icl<m_truthESClVCol.size(); icl++){
  auto truthMap = m_truthESClVCol[icl]->getLinkedMCP();
  printf("  In HFClusterV #%d: En %.3f, truthMC map size %d \n", icl, m_truthESClVCol[icl]->getEnergy(), truthMap.size());
  int counter = 0;
  for(auto imcp : truthMap){
    printf("      MCParticle #%d: truth En %.3f, mom (%.3f, %.3f, %.3f), weight %.3f \n", counter, imcp.first.getEnergy(), imcp.first.getMomentum().x, imcp.first.getMomentum().y, imcp.first.getMomentum().z, imcp.second );
    counter++;
  }
}

  //Doing matching. 
  for(int icl=0; icl<m_truthESClUCol.size(); icl++){
    for(int jcl=0; jcl<m_truthESClVCol.size(); jcl++){
      if(m_truthESClUCol[icl]->getLinkedMCP()[0].first == m_truthESClVCol[jcl]->getLinkedMCP()[0].first){
        std::shared_ptr<PandoraPlus::Calo3DCluster> tmp_clus = std::make_shared<PandoraPlus::Calo3DCluster>();
        XYClusterMatchingL0(m_truthESClUCol[icl].get(), m_truthESClVCol[jcl].get(), tmp_clus);
        m_clusters.push_back(tmp_clus);
      }
    }
  }

  return StatusCode::SUCCESS;
}


//Longitudinal cluster: 1*1
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL0( const PandoraPlus::CaloHalfCluster* m_longiClU, 
                                                       const PandoraPlus::CaloHalfCluster* m_longiClV, 
                                                       std::shared_ptr<PandoraPlus::Calo3DCluster>& m_clus )
{
//cout<<"  Cluster matching for case: 1 * 1. Input HalfCluster En: "<<m_longiClU->getEnergy()<<", "<<m_longiClV->getEnergy()<<endl;
//cout<<"  Print 1DShower En in HalfClusterU: "<<endl;
//for(int i=0; i<m_longiClU->getCluster().size(); i++)
//  cout<<m_longiClU->getCluster()[i]->getDlayer()<<'\t'<<m_longiClU->getCluster()[i]->getEnergy()<<endl;
//cout<<"  Print 1DShower En in HalfClusterV: "<<endl;
//for(int i=0; i<m_longiClV->getCluster().size(); i++)
//  cout<<m_longiClV->getCluster()[i]->getDlayer()<<'\t'<<m_longiClV->getCluster()[i]->getEnergy()<<endl;


  std::vector<int> layerindex; layerindex.clear(); 
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > map_showersUinlayer; map_showersUinlayer.clear();
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > map_showersVinlayer; map_showersVinlayer.clear();

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
    std::vector<const PandoraPlus::Calo1DCluster*> m_showerXcol = map_showersUinlayer[layerindex[il]];
    std::vector<const PandoraPlus::Calo1DCluster*> m_showerYcol = map_showersVinlayer[layerindex[il]];

//cout<<"  Print 1D showers in layer "<<layerindex[il]<<endl;
//cout<<"  ShowerU size = "<<m_showerXcol.size()<<endl;
//for(int a=0; a<m_showerXcol.size(); a++) printf("    #%d shower: En %.3f, address %p \n", a, m_showerXcol[a]->getEnergy(), m_showerXcol[a]);
//cout<<"  ShowerV size = "<<m_showerYcol.size()<<endl;
//for(int a=0; a<m_showerYcol.size(); a++) printf("    #%d shower: En %.3f, address %p \n", a, m_showerYcol[a]->getEnergy(), m_showerYcol[a]);

    std::vector<PandoraPlus::Calo2DCluster*> m_showerinlayer; m_showerinlayer.clear();

    if(m_showerXcol.size()==0 || m_showerYcol.size()==0) continue;
    //else if(m_showerXcol.size()==0){
    //  std::shared_ptr<PandoraPlus::Calo2DCluster> tmp_shower = std::make_shared<PandoraPlus::Calo2DCluster>();
    //  GetMatchedShowerFromEmpty(m_showerYcol[0], m_longiClU, tmp_shower.get());
    //  //m_showerinlayer.push_back(tmp_shower.get());
    //  //m_bkCol.map_2DCluster["bk2DCluster"].push_back(tmp_shower);
    //}
    //else if(m_showerYcol.size()==0){
    //  std::shared_ptr<PandoraPlus::Calo2DCluster> tmp_shower = std::make_shared<PandoraPlus::Calo2DCluster>();
    //  GetMatchedShowerFromEmpty(m_showerXcol[0], m_longiClU, tmp_shower.get());
    //}
    else if(m_showerXcol.size()==1 && m_showerYcol.size()==1){ 
      std::shared_ptr<PandoraPlus::Calo2DCluster> tmp_shower = std::make_shared<PandoraPlus::Calo2DCluster>();
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


//Longitudinal cluster: 1*N
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL1( const PandoraPlus::CaloHalfCluster* m_longiCl1, 
                                                       std::vector<const PandoraPlus::CaloHalfCluster*>& m_longiClN, 
                                                       std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>& m_clusters )
{
  m_clusters.clear(); 
  m_clusters.resize(m_longiClN.size());
  float totE_ClN = 0.;

  int slayer = m_longiCl1->getSlayer(); 
  std::vector<int> layerindex; layerindex.clear();
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > map_showersUinlayer; map_showersUinlayer.clear();
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > map_showersVinlayer; map_showersVinlayer.clear();

  for(int is=0; is<m_longiCl1->getCluster().size(); is++){
    int m_layer = m_longiCl1->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);

    if(slayer==0) map_showersUinlayer[m_layer].push_back(m_longiCl1->getCluster()[is]);
    else map_showersVinlayer[m_layer].push_back(m_longiCl1->getCluster()[is]);
  }
  for(int ic=0; ic<m_longiClN.size(); ic++){
  for(int is=0; is<m_longiClN[ic]->getCluster().size(); is++){
    int m_layer = m_longiClN[ic]->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);

    if(slayer==0) map_showersVinlayer[m_layer].push_back(m_longiClN[ic]->getCluster()[is]);
    else map_showersUinlayer[m_layer].push_back(m_longiClN[ic]->getCluster()[is]);
  }
  totE_ClN += m_longiClN[ic]->getEnergy();
  }
  
  std::map<int, std::vector<const PandoraPlus::Calo2DCluster*> > map_2Dshowersinlayer; map_2Dshowersinlayer.clear(); 

  sort(layerindex.begin(), layerindex.end());
  for(int il=0; il<layerindex.size(); il++){

    std::vector<const PandoraPlus::Calo1DCluster*> m_showerXcol = map_showersUinlayer[layerindex[il]];
    std::vector<const PandoraPlus::Calo1DCluster*> m_showerYcol = map_showersVinlayer[layerindex[il]];


    std::vector<PandoraPlus::Calo2DCluster*> m_showerinlayer; m_showerinlayer.clear();

    if(m_showerXcol.size()==0 || m_showerYcol.size()==0) continue;
    else if(m_showerXcol.size()==1 && m_showerYcol.size()==1){ 
      std::shared_ptr<PandoraPlus::Calo2DCluster> tmp_shower = std::make_shared<PandoraPlus::Calo2DCluster>();
      GetMatchedShowersL0(m_showerXcol[0], m_showerYcol[0], tmp_shower.get()); 
      m_showerinlayer.push_back(tmp_shower.get()); 
      m_bkCol.map_2DCluster["bk2DCluster"].push_back(tmp_shower);
    }
    else if(m_showerXcol.size()==1) GetMatchedShowersL1(m_showerXcol[0], m_showerYcol, m_showerinlayer );
    else if(m_showerYcol.size()==1) GetMatchedShowersL1(m_showerYcol[0], m_showerXcol, m_showerinlayer );
    else{ std::cout<<"CAUSION in XYClusterMatchingL1: HFCluster has ["<<m_showerXcol.size()<<", "<<m_showerYcol.size()<<"] showers in layer "<<layerindex[il]<<std::endl; }

    for(int is=0; is<m_showerinlayer.size(); is++) map_2Dshowersinlayer[layerindex[il]].push_back( m_showerinlayer[is] );
  }

  //Longitudinal linking for 2D clusters
  for(int il=0; il<layerindex.size(); il++){

  for(int is=0; is<map_2Dshowersinlayer[layerindex[il]].size(); is++){

    const PandoraPlus::Calo1DCluster* m_barshower; 
    if(slayer==1) m_barshower=map_2Dshowersinlayer[layerindex[il]][is]->getShowerUCol()[0];
    else          m_barshower=map_2Dshowersinlayer[layerindex[il]][is]->getShowerVCol()[0];

    int index_longiclus=-1; 
    bool fl_foundsh = false; 
    for(int ic=0; ic<m_longiClN.size() && !fl_foundsh; ic++){
    for(int jc=0; jc<m_longiClN[ic]->getCluster().size() && !fl_foundsh; jc++){
      if(m_barshower==m_longiClN[ic]->getCluster()[jc]) {index_longiclus=ic; fl_foundsh=true; break; }
    }}

    if(index_longiclus<0){ std::cout<<"WARNING: did not find properate longitudinal cluster! "<<std::endl; continue; }

    if( !m_clusters[index_longiclus].get() ){
      std::shared_ptr<PandoraPlus::Calo3DCluster> p_newclus = std::make_shared<PandoraPlus::Calo3DCluster>();
      p_newclus->addUnit(map_2Dshowersinlayer[layerindex[il]][is]);
      m_clusters[index_longiclus] = p_newclus;
    }
    else m_clusters[index_longiclus].get()->addUnit(map_2Dshowersinlayer[layerindex[il]][is]);
  }}

  //Split m_longiCl1 into N: 
  std::vector<const PandoraPlus::CaloHalfCluster*> m_splitHFClusN; m_splitHFClusN.clear();
  m_splitHFClusN.resize(m_longiClN.size());
  for(int icl=0; icl<m_longiClN.size(); icl++){
    float Eweight = m_longiClN[icl]->getEnergy()/totE_ClN;
    std::shared_ptr<PandoraPlus::CaloHalfCluster> newClus = std::make_shared<PandoraPlus::CaloHalfCluster>();
    for(int ish=0; ish<m_longiCl1->getCluster().size(); ish++){
      const PandoraPlus::Calo1DCluster* p_shower = m_longiCl1->getCluster()[ish];

      std::vector<const PandoraPlus::CaloUnit*> Bars; Bars.clear();
      for(int ibar=0; ibar<p_shower->getCluster().size(); ibar++){
        auto bar = p_shower->getCluster()[ibar]->Clone();
        bar->setQ(bar->getQ1()*Eweight, bar->getQ2()*Eweight );
        Bars.push_back(bar.get());
        m_bkCol.map_BarCol["bkBar"].push_back(bar);
      }

      std::shared_ptr<PandoraPlus::Calo1DCluster> shower = std::make_shared<PandoraPlus::Calo1DCluster>();
      shower->setBars(Bars);
      shower->setSeed();
      shower->setIDInfo();

      newClus->addUnit(shower.get());
      m_bkCol.map_1DCluster["bk1DCluster"].push_back( shower );
    }

    newClus->setHoughPars(m_longiCl1->getHoughAlpha(), m_longiCl1->getHoughRho());
    newClus->setIntercept(m_longiCl1->getHoughIntercept());
    for(auto iter: m_longiCl1->getHalfClusterMap()) newClus->setHalfClusters( iter.first, iter.second );    
    for(int itrk=0; itrk<m_longiCl1->getAssociatedTracks().size(); itrk++) newClus->addAssociatedTrack(m_longiCl1->getAssociatedTracks()[itrk]);
    newClus->setLinkedMCP(m_longiCl1->getLinkedMCP());
    newClus->setType(m_longiCl1->getType());

    m_splitHFClusN[icl]=newClus.get();
    m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(newClus);
  }

  //Match tracks 
  for(int ic=0; ic<m_clusters.size(); ic++){
    if(!m_clusters[ic].get() ) continue;
    if(slayer==0){ 
      m_clusters[ic].get()->addHalfClusterU( "LinkedLongiCluster", m_splitHFClusN[ic] );
      m_clusters[ic].get()->addHalfClusterV( "LinkedLongiCluster", m_longiClN[ic] );
      m_clusters[ic].get()->getLinkedMCPfromHFCluster("LinkedLongiCluster");

      for(auto itrk : m_longiCl1->getAssociatedTracks()){
        for(auto jtrk : m_longiClN[ic]->getAssociatedTracks()){
          if(itrk!=jtrk) continue; 

          if( find(m_clusters[ic].get()->getAssociatedTracks().begin(), m_clusters[ic].get()->getAssociatedTracks().end(), itrk)==m_clusters[ic].get()->getAssociatedTracks().end() )
            m_clusters[ic].get()->addAssociatedTrack(itrk);
      }}

    }
    else{          
      m_clusters[ic].get()->addHalfClusterU( "LinkedLongiCluster", m_longiClN[ic] );
      m_clusters[ic].get()->addHalfClusterV( "LinkedLongiCluster", m_splitHFClusN[ic] );
      m_clusters[ic].get()->getLinkedMCPfromHFCluster("LinkedLongiCluster");

      for(auto itrk : m_longiCl1->getAssociatedTracks()){
        for(auto jtrk : m_longiClN[ic]->getAssociatedTracks()){
          if(itrk!=jtrk) continue; 

          if( find(m_clusters[ic].get()->getAssociatedTracks().begin(), m_clusters[ic].get()->getAssociatedTracks().end(), itrk)==m_clusters[ic].get()->getAssociatedTracks().end() )
            m_clusters[ic].get()->addAssociatedTrack(itrk);
      }}
    }
  }

  return StatusCode::SUCCESS;
}

/*
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL2( std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClUCol,
                                                       std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClVCol,
                                                       std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>& m_clusters  )
{
  if( m_ClUCol.size()==0 || m_ClVCol.size()==0 ) return StatusCode::SUCCESS;


  std::vector<int> layerindex; layerindex.clear();
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > map_showersUinlayer; map_showersUinlayer.clear();
  std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> > map_showersVinlayer; map_showersVinlayer.clear();

  //Find layers need to match.
  for(int ic=0; ic<m_ClUCol.size(); ic++){
  for(int is=0; is<m_ClUCol[ic]->getCluster().size(); is++){
    int m_layer = m_ClUCol[ic]->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
    map_showersUinlayer[m_layer].push_back(m_ClUCol[ic]->getCluster()[is]);

  }}
  for(int ic=0; ic<m_ClVCol.size(); ic++){
  for(int is=0; is<m_ClVCol[ic]->getCluster().size(); is++){
    int m_layer = m_ClVCol[ic]->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
    map_showersVinlayer[m_layer].push_back(m_ClVCol[ic]->getCluster()[is]);
  }}
  sort(layerindex.begin(), layerindex.end());


  //Loop in layers to match. 
  std::map<int, std::vector<const PandoraPlus::Calo2DCluster*> > map_2Dshowersinlayer; map_2Dshowersinlayer.clear();
  for(int il=0; il<layerindex.size(); il++){

    std::vector<const PandoraPlus::Calo1DCluster*> m_showerXcol = map_showersUinlayer[layerindex[il]];
    std::vector<const PandoraPlus::Calo1DCluster*> m_showerYcol = map_showersVinlayer[layerindex[il]];
  
    std::vector<PandoraPlus::Calo2DCluster*> m_showerinlayer; m_showerinlayer.clear();

    if(m_showerXcol.size()==0 || m_showerYcol.size()==0) continue;
    else if(m_showerXcol.size()==1 && m_showerYcol.size()==1){
      std::shared_ptr<PandoraPlus::Calo2DCluster> tmp_shower = std::make_shared<PandoraPlus::Calo2DCluster>();
      GetMatchedShowersL0(m_showerXcol[0], m_showerYcol[0], tmp_shower.get());
      m_showerinlayer.push_back(tmp_shower.get());
      m_bkCol.map_2DCluster["bk2DCluster"].push_back(tmp_shower);
    }
    else if(m_showerXcol.size()==1) GetMatchedShowersL1(m_showerXcol[0], m_showerYcol, m_showerinlayer );
    else if(m_showerYcol.size()==1) GetMatchedShowersL1(m_showerYcol[0], m_showerXcol, m_showerinlayer );
    else if(m_showerXcol.size() == m_showerYcol.size()) GetMatchedShowersL2(m_showerXcol, m_showerYcol, m_showerinlayer );
    else GetMatchedShowersL3(m_showerXcol, m_showerYcol, m_showerinlayer);

    for(int is=0; is<m_showerinlayer.size(); is++) map_2Dshowersinlayer[layerindex[il]].push_back( m_showerinlayer[is] );
  }


  //Longitudinal linking for 2D clusters: 3D cone linking
  auto iter = map_2Dshowersinlayer.begin();
  std::vector<const PandoraPlus::Calo2DCluster*> HitsinFirstLayer = iter->second;
  for(int ihit=0; ihit<HitsinFirstLayer.size(); ihit++){
    std::shared_ptr<PandoraPlus::Calo3DCluster> p_newclus = std::make_shared<PandoraPlus::Calo3DCluster>();
    p_newclus->addUnit(HitsinFirstLayer[ihit]);
    m_clusters.push_back( p_newclus );
  }
  iter++;

  for(iter; iter!=map_2Dshowersinlayer.end(); iter++){
    std::vector<const PandoraPlus::Calo2DCluster*> HitsinLayer = iter->second;
    for(int is=0; is<HitsinLayer.size(); is++){
      const PandoraPlus::Calo2DCluster* m_hit = HitsinLayer[is];

      for(int ic=0; ic<m_clusters.size(); ic++ ){
        const PandoraPlus::Calo2DCluster* last_hit = m_clusters[ic]->getCluster().back();
        TVector3 relR_vec = m_hit->getPos() - last_hit->getPos();
        if( relR_vec.Angle(m_clusters[ic]->getShowerCenter())< settings.map_floatPars["th_ConeTheta"] && relR_vec.Mag()< settings.map_floatPars["th_ConeR"] ){
          m_clusters[ic].get()->addUnit(m_hit);
          HitsinLayer.erase(HitsinLayer.begin()+is);
          is--;
          break;
        }

      }
    }//end loop showers in layer.
    if(HitsinLayer.size()>0){
      for(int ihit=0;ihit<HitsinLayer.size(); ihit++){
        std::shared_ptr<PandoraPlus::Calo3DCluster> p_newclus = std::make_shared<PandoraPlus::Calo3DCluster>();
        p_newclus->addUnit(HitsinLayer[ihit]);
        m_clusters.push_back( p_newclus );
    }}//end new cluster
  }


  //Link HFClusters and tracks: 
  //  1. same ptr  2. deltaR < 1 cm. 
  for(int ic=0; ic<m_clusters.size(); ic++){
    const PandoraPlus::Calo1DCluster* p_showerU = m_clusters[ic].get()->getCluster()[0]->getShowerUCol()[0];
    bool fl_foundsh = false;
    int index_hfclusU = -1;
    for(int ihf=0; ihf<m_ClUCol.size() && !fl_foundsh; ihf++){
    for(int ish=0; ish<m_ClUCol[ihf]->getCluster().size() && !fl_foundsh; ish++){
      TVector3 delta_Pos = p_showerU->getPos() - m_ClUCol[ihf]->getCluster()[ish]->getPos(); 
      if(p_showerU == m_ClUCol[ihf]->getCluster()[ish] || delta_Pos.Mag()<1 ) { fl_foundsh=true; index_hfclusU = ihf; break; }
    }}

    if(index_hfclusU<0){ std::cout<<"ERROR: did not find properate HFClusterU! "<<std::endl; }
    else m_clusters[ic].get()->addHalfClusterU( "LinkedLongiCluster", m_ClUCol[index_hfclusU] );

    const PandoraPlus::Calo1DCluster* p_showerV = m_clusters[ic].get()->getCluster()[0]->getShowerVCol()[0];
    fl_foundsh = false;
    int index_hfclusV = -1;
    for(int ihf=0; ihf<m_ClVCol.size() && !fl_foundsh; ihf++){
    for(int ish=0; ish<m_ClVCol[ihf]->getCluster().size() && !fl_foundsh; ish++){
      TVector3 delta_Pos = p_showerV->getPos() - m_ClVCol[ihf]->getCluster()[ish]->getPos();
      if(p_showerV == m_ClVCol[ihf]->getCluster()[ish] || delta_Pos.Mag()<1 ) { fl_foundsh=true; index_hfclusV = ihf; break; }
    }}

    if(index_hfclusV<0){ std::cout<<"ERROR: did not find properate HFClusterV! "<<std::endl; }
    else m_clusters[ic].get()->addHalfClusterV( "LinkedLongiCluster", m_ClVCol[index_hfclusV] );

    for(auto itrk : m_ClUCol[index_hfclusU]->getAssociatedTracks()){
      for(auto jtrk : m_ClVCol[index_hfclusV]->getAssociatedTracks()){
        if(itrk!=jtrk) continue;

        if( find(m_clusters[ic].get()->getAssociatedTracks().begin(), m_clusters[ic].get()->getAssociatedTracks().end(), itrk)==m_clusters[ic].get()->getAssociatedTracks().end() )
          m_clusters[ic].get()->addAssociatedTrack(itrk);
    }}

  }

  return StatusCode::SUCCESS;
}
*/

//Longitudinal cluster: N*N
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL2( std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClUCol,
                                                       std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClVCol,
                                                       std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>& m_clusters  )
{
//cout<<"  XYClusterMatchingL2: Input cluster size: "<<m_ClUCol.size()<<'\t'<<m_ClVCol.size()<<endl;
  if(m_ClUCol.size()==0 || m_ClVCol.size()==0 || m_ClUCol.size()!=m_ClVCol.size()) return StatusCode::SUCCESS;

  std::vector<int> layerindex; layerindex.clear();
  std::map<int, std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> > map_showersXinlayer; map_showersXinlayer.clear();
  std::map<int, std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> > map_showersYinlayer; map_showersYinlayer.clear();


  //Find layers need to match.
  for(int ic=0; ic<m_ClUCol.size(); ic++){
  for(int is=0; is<m_ClUCol[ic]->getCluster().size(); is++){
    int m_layer = m_ClUCol[ic]->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
  }}
  for(int ic=0; ic<m_ClVCol.size(); ic++){
  for(int is=0; is<m_ClVCol[ic]->getCluster().size(); is++){
    int m_layer = m_ClVCol[ic]->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
  }}
  sort(layerindex.begin(), layerindex.end());

  //Fill shower maps
  for(int il=0; il<layerindex.size(); il++){
    for(int ic=0; ic<m_ClUCol.size(); ic++){
      map_showersXinlayer[layerindex[il]].push_back( m_ClUCol[ic]->getClusterInLayer(layerindex[il]) );
    }
    for(int ic=0; ic<m_ClVCol.size(); ic++){
      map_showersYinlayer[layerindex[il]].push_back( m_ClVCol[ic]->getClusterInLayer(layerindex[il]) );
    }
  }

  //Get the chi2 map for N*N
  const int Nclus = m_ClUCol.size();
  double  **map_chi2[PandoraPlus::CaloUnit::Nlayer] = {NULL};
  double sumchi2[Nclus][Nclus] = {0};

  for(int il=0; il<layerindex.size(); il++){
    std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> m_showerXcol = map_showersXinlayer[layerindex[il]];
    std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> m_showerYcol = map_showersYinlayer[layerindex[il]];

    //double **chi2inlayer = GetClusterChi2Map(m_showerXcol, m_showerYcol);
    //if(chi2inlayer!=nullptr) map_chi2[layerindex[il]-1] = chi2inlayer;
    map_chi2[layerindex[il]-1] = GetClusterChi2Map(m_showerXcol, m_showerYcol);
  }

//cout<<"  XYClusterMatchingL2: Print Sumchi2 matrix"<<endl;
  for(int ic=0; ic<Nclus; ic++){
  for(int jc=0; jc<Nclus; jc++){
    for(int il=0; il<PandoraPlus::CaloUnit::Nlayer; il++){
      if(map_chi2[il]==nullptr) continue;
      sumchi2[ic][jc] += map_chi2[il][ic][jc];
    }
//cout<<sumchi2[ic][jc]<<'\t';
  }
//cout<<endl;
  }

  //Get the chi2 of N! combinations
  int Ncomb=1;
  for(int i=Nclus; i>0; i--) Ncomb = Ncomb*i;

  map<double, vector<pair<int, int>> > matchingMap;
  int num[Nclus];
  int num_init[Nclus];
  for(int i=0;i<Nclus;i++){ num[i]=i; num_init[i]=i;}

  for(int icont=0;icont<Ncomb;icont++){
    vector<pair<int, int>> Index;
    for(int i=0;i<Nclus;i++){
       pair<int, int> p1(num_init[i], num[i]);
       Index.push_back(p1);
    }
    double chi2_tot=0;
    for(int i=0;i<Index.size();i++) chi2_tot += sumchi2[Index[i].first][Index[i].second];
    matchingMap[chi2_tot] = Index;

    Index.clear();
    if(!next_permutation(num, num+Nclus)) break;
  }

//cout<<"  XYClusterMatchingL2: Print chi2 map"<<endl;
//for(auto it : matchingMap){
//  vector<pair<int, int>> m_index = it.second;
//  cout<<"    In combination: ";
//  for(auto num : m_index) printf("[%d, %d] + ", num.first, num.second);
//  cout<<". Chi2 = "<<it.first<<endl;
//}

  //map is ordered with [double] value, first element has the smallest chi2 value.
  map<double, vector<pair<int, int>> >::iterator iter = matchingMap.begin();
  vector<pair<int, int>> Index = iter->second;

  for(int ii=0; ii<Index.size(); ii++){
    const PandoraPlus::CaloHalfCluster* m_clusX = m_ClUCol[Index[ii].first];
    const PandoraPlus::CaloHalfCluster* m_clusY = m_ClVCol[Index[ii].second];
    std::shared_ptr<PandoraPlus::Calo3DCluster> tmp_clus = std::make_shared<PandoraPlus::Calo3DCluster>();

    XYClusterMatchingL0(m_clusX, m_clusY, tmp_clus);
    m_clusters.push_back(tmp_clus);
  }

  for(int il=0; il<PandoraPlus::CaloUnit::Nlayer; il++) delete map_chi2[il];
  return StatusCode::SUCCESS;
}


//Longitudinal cluster: M*N
StatusCode EnergyTimeMatchingAlg::XYClusterMatchingL3( std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClUCol,
                                                       std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClVCol,
                                                       std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>& m_clusters )
{
cout<<"  In XYClusterMatchingL3: Input cluster size "<<m_ClUCol.size()<<'\t'<<m_ClVCol.size()<<endl;
  if( m_ClUCol.size()==0 || m_ClVCol.size()==0 ) return StatusCode::SUCCESS;

  std::vector<int> layerindex; layerindex.clear();
  std::map<int, std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> > map_showersXinlayer; map_showersXinlayer.clear();
  std::map<int, std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> > map_showersYinlayer; map_showersYinlayer.clear();

  //Find layers need to match.
  for(int ic=0; ic<m_ClUCol.size(); ic++){
  for(int is=0; is<m_ClUCol[ic]->getCluster().size(); is++){
    int m_layer = m_ClUCol[ic]->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
  }}
  for(int ic=0; ic<m_ClVCol.size(); ic++){
  for(int is=0; is<m_ClVCol[ic]->getCluster().size(); is++){
    int m_layer = m_ClVCol[ic]->getCluster()[is]->getDlayer();
    if( find( layerindex.begin(), layerindex.end(), m_layer )==layerindex.end() ) layerindex.push_back(m_layer);
  }}
  sort(layerindex.begin(), layerindex.end());

  //Fill shower maps
  for(int il=0; il<layerindex.size(); il++){
    for(int ic=0; ic<m_ClUCol.size(); ic++)
      map_showersXinlayer[layerindex[il]].push_back( m_ClUCol[ic]->getClusterInLayer(layerindex[il]) );
    for(int ic=0; ic<m_ClVCol.size(); ic++)
      map_showersYinlayer[layerindex[il]].push_back( m_ClVCol[ic]->getClusterInLayer(layerindex[il]) );
  }
cout<<"  XYClusterMatchingL3: map size X "<<map_showersXinlayer.size()<<", Y "<<map_showersYinlayer.size()<<endl;

  //Get the chi2 map
  const int NclusU = m_ClUCol.size();
  const int NclusV = m_ClVCol.size();
  double  **map_chi2[PandoraPlus::CaloUnit::Nlayer] = {NULL};
  double sumchi2[NclusU][NclusV] = {0};

  for(int il=0; il<layerindex.size(); il++){
    std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> m_showerXcol = map_showersXinlayer[layerindex[il]];
    std::vector<std::vector<const PandoraPlus::Calo1DCluster*>> m_showerYcol = map_showersYinlayer[layerindex[il]];
    map_chi2[layerindex[il]-1] = GetClusterChi2Map(m_showerXcol, m_showerYcol);
  }

cout<<"  XYClusterMatchingL3: Print Sumchi2 matrix"<<endl;
  map<double, pair<int, int> > m_chi2Map; m_chi2Map.clear();
  for(int ic=0; ic<NclusU; ic++){
  for(int jc=0; jc<NclusV; jc++){
    for(int il=0; il<PandoraPlus::CaloUnit::Nlayer; il++){
      if(map_chi2[il]==nullptr) continue;

      sumchi2[ic][jc] += map_chi2[il][ic][jc];
    }
    pair<int, int> p1(ic, jc);
    m_chi2Map[sumchi2[ic][jc]] = p1;
cout<<sumchi2[ic][jc]<<'\t';
  }
cout<<endl;
  }

  //pickout pairs <indexVec> with smallest chi2:
  pair<int, int> lastpair;
  vector<pair<int, int>> indexVec; indexVec.clear();
  auto iter = m_chi2Map.begin();
  for(iter; iter!=m_chi2Map.end(); iter++){
    pair<int, int> indexpair = iter->second;
    bool inLine = false;
    bool isLast = false;
    for(int i=0; i<indexVec.size(); i++){
      if( indexpair.first == indexVec[i].first || indexpair.second == indexVec[i].second ) { inLine=true; break; }
      if( indexVec.size() == min(NclusU, NclusV)-1 ) {lastpair = indexpair; isLast=true;  break; }
    }
    if(isLast) break;
    if(inLine) continue;
    indexVec.push_back(indexpair);
  }
  if( indexVec.size()!=min(NclusU, NclusV)-1 )
    std::cout<<"ERROR in XYClusterMatchingL3: found pair size "<<indexVec.size()<<" does not equal to min shower size -1 "<<min(NclusU, NclusV)-1<<endl;

  //For the left: 1*N
  std::vector<const PandoraPlus::CaloHalfCluster*> leftHFClusCol; leftHFClusCol.clear();
  for(int i=0; i<max(NclusU, NclusV); i++){
    bool fl_exist = false;
    for(int j=0; j<indexVec.size(); j++){
      int m_index = NclusU>NclusV ? indexVec[j].first : indexVec[j].second;
      if(i==m_index){ fl_exist = true; break; }
    }
    if(!fl_exist){
       const PandoraPlus::CaloHalfCluster* m_shower = NclusU>NclusV ? m_ClUCol[i] : m_ClVCol[i] ;
       leftHFClusCol.push_back(m_shower);
    }
  }
  if(leftHFClusCol.size() != fabs( NclusU-NclusV )+1 )
    std::cout<<"ERROR in XYClusterMatchingL3: Last pair number "<<leftHFClusCol.size()<<" does not equal to shower difference "<<fabs( NclusU-NclusV )+1<<std::endl;

  //Match the pairs in the indexVec:
cout<<"  XYClusterMatchingL3: Match the pairs in the indexVec: size = "<<indexVec.size()<<endl;
  for(int ip=0; ip<indexVec.size(); ip++){
    std::shared_ptr<PandoraPlus::Calo3DCluster> tmp_clus = std::make_shared<PandoraPlus::Calo3DCluster>();

printf("    Matching pair (%d, %d): HFClusterU Pos+E (%.3f, %.3f, %.3f, %.3f). HFClusterV Pos+E (%.3f, %.3f, %.3f, %.3f) \n", indexVec[ip].first, indexVec[ip].second, 
  m_ClUCol[indexVec[ip].first]->getPos().x(), m_ClUCol[indexVec[ip].first]->getPos().y(), m_ClUCol[indexVec[ip].first]->getPos().z(), m_ClUCol[indexVec[ip].first]->getEnergy(), 
  m_ClVCol[indexVec[ip].second]->getPos().x(), m_ClVCol[indexVec[ip].second]->getPos().y(), m_ClVCol[indexVec[ip].second]->getPos().z(), m_ClVCol[indexVec[ip].second]->getEnergy());

    XYClusterMatchingL0(m_ClUCol[indexVec[ip].first], m_ClVCol[indexVec[ip].second], tmp_clus);
printf("    Rec cluster Pos+E (%.3f, %.3f, %.3f, %.3f) \n", tmp_clus->getShowerCenter().x(), tmp_clus->getShowerCenter().y(), tmp_clus->getShowerCenter().z(), tmp_clus->getEnergy());
    m_clusters.push_back(tmp_clus);
  }

cout<<"  Last pair: "<<lastpair.first<<" "<<lastpair.second<<endl;
  //Match the left 1*N:
cout<<"  XYClusterMatchingL3: Match the left 1*"<<leftHFClusCol.size();
  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> tmp_3dclus; tmp_3dclus.clear();
  int ilast = NclusU<NclusV ? lastpair.first : lastpair.second;
cout<<" Index of 1: "<<ilast<<endl;
  const PandoraPlus::CaloHalfCluster* m_clus = NclusU<NclusV ? m_ClUCol[ilast] : m_ClVCol[ilast];
  XYClusterMatchingL1(m_clus, leftHFClusCol, tmp_3dclus );

  m_clusters.insert(m_clusters.end(), tmp_3dclus.begin(), tmp_3dclus.end());
  m_clus = nullptr;

cout<<"  Delete chi2 map ptr"<<endl;
  for(int il=0; il<PandoraPlus::CaloUnit::Nlayer; il++) delete map_chi2[il];
  return StatusCode::SUCCESS;
}

/*
StatusCode EnergyTimeMatchingAlg::GetMatchedShowerFromEmpty( const PandoraPlus::Calo1DCluster* barShower,
                                                             const PandoraPlus::CaloHalfCluster* m_longiCl,
                                                             PandoraPlus::Calo2DCluster* outsh )
{
  if(barShower->getBars().size()==0){
    std::cout<<"ERROR: empty 1DCluster in GetMatchedShowerFromEmpty. Empty DigiHitsCol returned! "<<std::endl;
    return StatusCode::SUCCESS;
  }

  int _module = barShower->getTowerID()[0][0];
  int _dlayer = barShowerU->getBars()[0]->getDlayer();
  int _slayer = barShower->getBars()[0]->getSlayer();

  float rotAngle = -_module*TMath::Pi()/4.;

  //Create a 1DCluster from 1DCluster in neighbor layers in HalfCluster: +-2 layers. 
  std::shared_ptr<CaloUnit> m_bar; //Clone from the seed of most energetic 1DCluster. 
  for(int il=-2; il<3; il++){
    if(m_longiCl->getClusterInLayer(_dlayer+il).size()!=0){
      float maxE = -99;
      int index = -1;
      for(int ish=0; ish<m_longiCl->getClusterInLayer(_dlayer+il).size(); ish++){
        if(m_longiCl->getClusterInLayer(_dlayer+il)[ish]->getEnergy()>maxE){
          maxE=m_longiCl->getClusterInLayer(_dlayer+il)[ish]->getEnergy();
          index = ish;
        }
      }
      if(index>=0){
        m_bar = m_longiCl->getClusterInLayer(_dlayer+il)[index]->getSeeds()[0]->Clone();
        m_bar.get()->SetQ(barShower->getEnergy()/2., barShower->getEnergy()/2.);
        m_bar.get()->setPosition();
        m_bar.get()->setcellID();
        break;
      }
    }
  }
  if(!m_bar.get()) {
cout<<"  No cluster in +-2 layers in origin HalfCluster. Skip! "<<endl;
    return StatusCode::SUCCESS;
  }



  return StatusCode::SUCCESS;
}
*/



StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL0( const PandoraPlus::Calo1DCluster* barShowerU,   
                                                       const PandoraPlus::Calo1DCluster* barShowerV,
                                                       PandoraPlus::Calo2DCluster* outsh )
{

  std::vector<const PandoraPlus::CaloHit*> m_digiCol; m_digiCol.clear();
  int NbarsX = barShowerU->getBars().size();
  int NbarsY = barShowerV->getBars().size();
  if(NbarsX==0 || NbarsY==0){ std::cout<<"WARNING: empty DigiHitsCol returned!"<<std::endl; return StatusCode::SUCCESS; }
  if(barShowerU->getTowerID().size()==0) { std::cout<<"WARNING:GetMatchedShowersL0  No TowerID in 1DCluster!"<<std::endl; return StatusCode::SUCCESS; }
  //if(barShowerU->getTowerID().size()==0) { barShowerU->setIDInfo(); }

  int _module = barShowerU->getTowerID()[0][0];
  float rotAngle = -_module*TMath::Pi()/4.;

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
      std::shared_ptr<PandoraPlus::CaloHit> hit = std::make_shared<PandoraPlus::CaloHit>();
      //hit.setCellID(0);
      hit->setPosition(p_hit);
      hit->setEnergy(m_Ehit);
      m_digiCol.push_back(hit.get());
      m_bkCol.map_CaloHit["bkHit"].push_back( hit );
    }
  }

  outsh->addUnit( barShowerU );
  outsh->addUnit( barShowerV );
  outsh->setCaloHits( m_digiCol );
//cout<<"    End output shower"<<endl;

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL1( const PandoraPlus::Calo1DCluster* shower1,
                                                       std::vector<const PandoraPlus::Calo1DCluster*>& showerNCol, 
                                                       std::vector<PandoraPlus::Calo2DCluster*>& outshCol )
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
    std::shared_ptr<PandoraPlus::Calo1DCluster> m_splitshower1 = std::make_shared<PandoraPlus::Calo1DCluster>();
    m_bkCol.map_1DCluster["bk1DCluster"].push_back( m_splitshower1 );

    std::shared_ptr<PandoraPlus::CaloUnit> m_wiseed = nullptr; 
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

    std::vector<const PandoraPlus::CaloUnit*> m_wibars; m_wibars.clear(); 
    for(int ib=0;ib<shower1->getBars().size();ib++){
      std::shared_ptr<PandoraPlus::CaloUnit> m_wibar = shower1->getBars()[ib]->Clone();
      m_wibar->setQ(wi_E*m_wibar->getQ1(), wi_E*m_wibar->getQ2());
      m_wibars.push_back(m_wibar.get());
      m_bkCol.map_BarCol["bkBar"].push_back( m_wibar );
    }
    m_splitshower1->setBars( m_wibars );
    m_splitshower1->addSeed( m_wiseed.get() );
    m_splitshower1->setIDInfo();
    std::shared_ptr<PandoraPlus::Calo2DCluster> m_shower = std::make_shared<PandoraPlus::Calo2DCluster>();
    if(_slayer==0 ) GetMatchedShowersL0( m_splitshower1.get(), showerNCol[is], m_shower.get() );
    else            GetMatchedShowersL0( showerNCol[is], m_splitshower1.get(), m_shower.get() );

    outshCol.push_back( m_shower.get() );
    m_bkCol.map_2DCluster["bk2DCluster"].push_back( m_shower );
  }
  return StatusCode::SUCCESS;
}

/*
StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL2( std::vector<const PandoraPlus::Calo1DCluster*>& barShowerUCol,
                                                       std::vector<const PandoraPlus::Calo1DCluster*>& barShowerVCol,
                                                       std::vector<PandoraPlus::Calo2DCluster*>& outshCol )
{
  outshCol.clear();
  if(barShowerUCol.size() != barShowerVCol.size() ) return StatusCode::FAILURE;

  const int Nshower = barShowerUCol.size();
  double chi2[Nshower][Nshower];
  double chi2_E[Nshower][Nshower];
  double chi2_tx[Nshower][Nshower];
  double chi2_ty[Nshower][Nshower];

  double wi_E = settings.map_floatPars["chi2Wi_E"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);
  double wi_T = settings.map_floatPars["chi2Wi_T"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);

  TVector3 m_vec(0,0,0);
  double rotAngle = -(barShowerUCol[0]->getBars())[0]->getModule()*TMath::Pi()/4.;
  TVector3 Cblock( (barShowerUCol[0]->getBars())[0]->getPosition().x(),
                   (barShowerUCol[0]->getBars())[0]->getPosition().y(),
                   (barShowerVCol[0]->getBars())[0]->getPosition().z() );
  Cblock.RotateZ(rotAngle);

  for(int ix=0;ix<Nshower;ix++){
  for(int iy=0;iy<Nshower;iy++){
    const PandoraPlus::Calo1DCluster* showerX = barShowerUCol[ix];
    const PandoraPlus::Calo1DCluster* showerY = barShowerVCol[iy];

    double Ex = showerX->getEnergy();
    double Ey = showerY->getEnergy();
    chi2_E[ix][iy] = pow(fabs(Ex-Ey)/settings.map_floatPars["sigmaE"], 2);
    double PosTx = C*(showerY->getT1()-showerY->getT2())/(2*settings.map_floatPars["nMat"]) + showerY->getPos().z();
    chi2_tx[ix][iy] = pow( fabs(PosTx-showerX->getPos().z())/settings.map_floatPars["sigmaPos"], 2 );

    double PosTy = C*(showerX->getT1()-showerX->getT2())/(2*settings.map_floatPars["nMat"]);
    m_vec.SetXYZ(showerY->getPos().x(), showerY->getPos().y(), showerY->getPos().z());
    m_vec.RotateZ(rotAngle);
    chi2_ty[ix][iy] = pow( fabs(PosTy - (m_vec-Cblock).x() )/settings.map_floatPars["sigmaPos"], 2);

    chi2[ix][iy] = chi2_E[ix][iy]*wi_E + (chi2_tx[ix][iy]+chi2_ty[ix][iy])*wi_T ;

  }}

  int Ncomb=1;
  for(int i=Nshower; i>0; i--) Ncomb = Ncomb*i;

  map<double, vector<pair<int, int>> > matchingMap;
  int num[Nshower];
  int num_init[Nshower];
  for(int i=0;i<Nshower;i++){ num[i]=i; num_init[i]=i;}

  for(int icont=0;icont<Ncomb;icont++){
    vector<pair<int, int>> Index;
    for(int i=0;i<Nshower;i++){
       pair<int, int> p1(num_init[i], num[i]);
       Index.push_back(p1);
    }
    double chi2_tot=0;
    for(int i=0;i<Index.size();i++) chi2_tot += chi2[Index[i].first][Index[i].second];
    matchingMap[chi2_tot] = Index;

    Index.clear();
    if(!next_permutation(num, num+Nshower)) break;
  }

  map<double, vector<pair<int, int>> >::iterator iter = matchingMap.begin();
  vector<pair<int, int>> Index = iter->second;

  for(int i=0;i<Index.size();i++){
    const PandoraPlus::Calo1DCluster* showerX = barShowerUCol[Index[i].first];
    const PandoraPlus::Calo1DCluster* showerY = barShowerVCol[Index[i].second];

    std::shared_ptr<PandoraPlus::Calo2DCluster> tmp_shower = std::make_shared<PandoraPlus::Calo2DCluster>();
    GetMatchedShowersL0(showerX, showerY, tmp_shower.get());
    outshCol.push_back(tmp_shower.get());
    m_bkCol.map_2DCluster["bk2DCluster"].push_back( tmp_shower );
  }

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL3(  std::vector<const PandoraPlus::Calo1DCluster*>& barShowerUCol,
                                                        std::vector<const PandoraPlus::Calo1DCluster*>& barShowerVCol,
                                                        std::vector<PandoraPlus::Calo2DCluster*>& outshCol )
{
  outshCol.clear();

  const int NshowerU = barShowerUCol.size();
  const int NshowerV = barShowerVCol.size();

  double chi2[NshowerU][NshowerV];
  double chi2_E[NshowerU][NshowerV];
  double chi2_tx[NshowerU][NshowerV];
  double chi2_ty[NshowerU][NshowerV];

  double wi_E = settings.map_floatPars["chi2Wi_E"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);
  double wi_T = settings.map_floatPars["chi2Wi_T"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);

  TVector3 m_vec(0,0,0);
  double rotAngle = -(barShowerUCol[0]->getBars())[0]->getModule()*TMath::Pi()/4.;
  TVector3 Cblock( (barShowerUCol[0]->getBars())[0]->getPosition().x(),
                   (barShowerUCol[0]->getBars())[0]->getPosition().y(),
                   (barShowerVCol[0]->getBars())[0]->getPosition().z() );
  Cblock.RotateZ(rotAngle);
  map<double, pair<int, int> > m_chi2Map; m_chi2Map.clear();

  for(int ix=0;ix<NshowerU;ix++){
  for(int iy=0;iy<NshowerV;iy++){
    const PandoraPlus::Calo1DCluster* showerX = barShowerUCol[ix];
    const PandoraPlus::Calo1DCluster* showerY = barShowerVCol[iy];

    double Ex = showerX->getEnergy();
    double Ey = showerY->getEnergy();
    chi2_E[ix][iy] = pow(fabs(Ex-Ey)/settings.map_floatPars["sigmaE"], 2);
    double PosTx = C*(showerY->getT1()-showerY->getT2())/(2*settings.map_floatPars["nMat"]) + showerY->getPos().z();
    chi2_tx[ix][iy] = pow( fabs(PosTx-showerX->getPos().z())/settings.map_floatPars["sigmaPos"], 2 );

    double PosTy = C*(showerX->getT1()-showerX->getT2())/(2*settings.map_floatPars["nMat"]);
    m_vec.SetXYZ(showerY->getPos().x(), showerY->getPos().y(), showerY->getPos().z());
    m_vec.RotateZ(rotAngle);
    chi2_ty[ix][iy] = pow( fabs(PosTy - (m_vec-Cblock).x() )/settings.map_floatPars["sigmaPos"], 2);

    chi2[ix][iy] = chi2_E[ix][iy]*wi_E + (chi2_tx[ix][iy]+chi2_ty[ix][iy])*wi_T ;

    pair<int, int> p1(ix, iy);
    m_chi2Map[chi2[ix][iy]] = p1;
  }}

  pair<int, int> lastpair;
  vector<pair<int, int>> indexVec; indexVec.clear();
  map<double, pair<int, int> >::iterator iter = m_chi2Map.begin();

  for(iter; iter!=m_chi2Map.end(); iter++){
    pair<int, int> indexpair = iter->second;
    bool inLine = false;
    bool isLast = false;
    for(int i=0; i<indexVec.size(); i++){
      if( indexVec.size() == min(NshowerU, NshowerV)-1 ) {lastpair = indexpair; isLast=true;  break; }
      if( indexpair.first == indexVec[i].first || indexpair.second == indexVec[i].second ) { inLine=true; break; }
    }
    if(isLast) break;
    if(inLine) continue;
    indexVec.push_back(indexpair);
  }
  if( indexVec.size()!=min(NshowerU, NshowerV)-1 )
    cout<<"ERROR in EnergyTimeMatchingAlg::GetMatchedShowersL3: found pair size "<<indexVec.size()<<" does not equal to min shower size -1 "<<min(NshowerU, NshowerV)-1<<endl;

  vector<const PandoraPlus::Calo1DCluster*> leftShowers; leftShowers.clear();
  for(int i=0; i<max(NshowerU, NshowerV); i++){
    bool fl_exist = false;
    for(int j=0; j<indexVec.size(); j++){
      int m_index = NshowerU>NshowerV ? indexVec[j].first : indexVec[j].second;
      if(i==m_index){ fl_exist = true; break; }
    }
    if(!fl_exist){
       const PandoraPlus::Calo1DCluster* m_shower = NshowerU>NshowerV ? barShowerUCol[i] : barShowerVCol[i] ;
       leftShowers.push_back(m_shower);
    }
  }
  if(leftShowers.size() != fabs( NshowerU-NshowerV )+1 )
    cout<<"ERROR in XYShowerChi2MatchingL1: Last pair number "<<leftShowers.size()<<" does not equal to shower difference "<<fabs( NshowerU-NshowerV )+1<<endl;
  for(int ip=0; ip<indexVec.size(); ip++){
    const PandoraPlus::Calo1DCluster* showerX = barShowerUCol[indexVec[ip].first];
    const PandoraPlus::Calo1DCluster* showerY = barShowerVCol[indexVec[ip].second];

    std::shared_ptr<PandoraPlus::Calo2DCluster> tmp_shower = std::make_shared<PandoraPlus::Calo2DCluster>();
    GetMatchedShowersL0(showerX, showerY, tmp_shower.get());
    outshCol.push_back(tmp_shower.get());
    m_bkCol.map_2DCluster["bk2DCluster"].push_back( tmp_shower );
  }


  std::vector<PandoraPlus::Calo2DCluster*> m_showerinlayer; m_showerinlayer.clear();
  int ilast = NshowerU<NshowerV ? lastpair.first : lastpair.second;
  const PandoraPlus::Calo1DCluster* m_shower = NshowerU<NshowerV ? barShowerUCol[ilast] : barShowerVCol[ilast] ;
  GetMatchedShowersL1( m_shower, leftShowers, m_showerinlayer);

  outshCol.insert(outshCol.end(), m_showerinlayer.begin(), m_showerinlayer.end());


  return StatusCode::SUCCESS;
}
*/

double** EnergyTimeMatchingAlg::GetClusterChi2Map( std::vector<std::vector<const PandoraPlus::Calo1DCluster*>>& barShowerUCol,
                                                   std::vector<std::vector<const PandoraPlus::Calo1DCluster*>>& barShowerVCol )
{

  const int NclusX = barShowerUCol.size();
  const int NclusY = barShowerVCol.size();

  if(NclusX==0 || NclusY==0) return nullptr;
  //If one longidutinal cluster is empty in this layer: skip this layer.
  //for(int icx=0; icx<NclusX; icx++) if(barShowerUCol[icx].size()==0) return nullptr;
  //for(int icy=0; icy<NclusY; icy++) if(barShowerVCol[icy].size()==0) return nullptr;

  double  **chi2map = new double*[NclusX];

  double chi2map_E[NclusX][NclusY]={0};
  double chi2map_tx[NclusX][NclusY]={0};
  double chi2map_ty[NclusX][NclusY]={0};

  double wi_E = settings.map_floatPars["chi2Wi_E"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);
  double wi_T = settings.map_floatPars["chi2Wi_T"]/(settings.map_floatPars["chi2Wi_E"] + settings.map_floatPars["chi2Wi_T"]);

  TVector3 m_vec(0,0,0);
  double rotAngle = -999;
  TVector3 Ctower(0,0,0);
  for(int ish=0; ish<barShowerUCol.size(); ish++){
    if(barShowerUCol[ish].size()==0) continue;
    rotAngle = -(barShowerUCol[ish][0]->getBars())[0]->getModule()*TMath::Pi()/4.;
    Ctower.SetX( (barShowerUCol[ish][0]->getBars())[0]->getPosition().x() );
    Ctower.SetY( (barShowerUCol[ish][0]->getBars())[0]->getPosition().y() );
  }
  for(int ish=0; ish<barShowerVCol.size(); ish++){
    if(barShowerVCol[ish].size()==0) continue;
    Ctower.SetZ( (barShowerVCol[ish][0]->getBars())[0]->getPosition().z() );
  }
  Ctower.RotateZ(rotAngle);

  for(int ix=0;ix<NclusX;ix++){
  chi2map[ix] = new double[NclusY];
  for(int iy=0;iy<NclusY;iy++){
    std::vector<const PandoraPlus::Calo1DCluster*> clusterU = barShowerUCol[ix];
    std::vector<const PandoraPlus::Calo1DCluster*> clusterV = barShowerVCol[iy];

    if(clusterU.size()==0){
      double totE_V = 0; 
      for(int icy=0; icy<clusterV.size(); icy++) totE_V += clusterV[icy]->getEnergy();
      chi2map[ix][iy] = wi_E*pow(totE_V/settings.map_floatPars["sigmaE"], 2); 
      continue;
    }
    if(clusterV.size()==0){
      double totE_U = 0;
      for(int icy=0; icy<clusterU.size(); icy++) totE_U += clusterU[icy]->getEnergy();
      chi2map[ix][iy] = wi_E*pow(totE_U/settings.map_floatPars["sigmaE"], 2);
      continue;
    }

    double min_chi2E = 999;
    double min_chi2tx = 999;
    double min_chi2ty = 999;

    for(int icx=0; icx<clusterU.size(); icx++){
    for(int icy=0; icy<clusterV.size(); icy++){
      const PandoraPlus::Calo1DCluster* showerX = clusterU[icx];
      const PandoraPlus::Calo1DCluster* showerY = clusterV[icy];

      double Ex = showerX->getEnergy();
      double Ey = showerY->getEnergy();
      double chi2_E = pow(fabs(Ex-Ey)/settings.map_floatPars["sigmaE"], 2);
      double PosTx = C*(showerY->getT1()-showerY->getT2())/(2*settings.map_floatPars["nMat"]) + showerY->getPos().z();
      double chi2_tx = pow( fabs(PosTx-showerX->getPos().z())/settings.map_floatPars["sigmaPos"], 2 );

      double PosTy = C*(showerX->getT1()-showerX->getT2())/(2*settings.map_floatPars["nMat"]);
      m_vec = showerY->getPos();
      m_vec.RotateZ(rotAngle);
      double chi2_ty = pow( fabs(PosTy - (m_vec-Ctower).x() )/settings.map_floatPars["sigmaPos"], 2);

      if(chi2_E<min_chi2E) min_chi2E=chi2_E;
      if(chi2_tx<min_chi2tx) min_chi2tx=chi2_tx;
      if(chi2_ty<min_chi2ty) min_chi2ty=chi2_ty;
    }}

    chi2map_E[ix][iy] = min_chi2E;
    chi2map_tx[ix][iy] = min_chi2tx;
    chi2map_ty[ix][iy] = min_chi2ty;
    chi2map[ix][iy] = chi2map_E[ix][iy]*wi_E + (chi2map_tx[ix][iy]+chi2map_ty[ix][iy])*wi_T ;

  }}

  return chi2map;

};

#endif
