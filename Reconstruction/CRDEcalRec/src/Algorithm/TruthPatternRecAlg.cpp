#ifndef _TRUTHPATTERNREC_ALG_C
#define _TRUTHPATTERNREC_ALG_C

#include "Algorithm/TruthPatternRecAlg.h"

StatusCode TruthPatternRecAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_stringPars.find("ReadinLocalMaxName")==settings.map_stringPars.end())
    settings.map_stringPars["ReadinLocalMaxName"] = "AllLocalMax";
  if(settings.map_stringPars.find("OutputLongiClusName")==settings.map_stringPars.end())
    settings.map_stringPars["OutputLongiClusName"] = "TruthAxis";    


  return StatusCode::SUCCESS;
};

StatusCode TruthPatternRecAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  p_HalfClusterU.clear();
  p_HalfClusterV.clear();

  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColU"].size(); ih++)
    p_HalfClusterU.push_back( m_datacol.map_HalfCluster["HalfClusterColU"][ih].get() );
  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColV"].size(); ih++)
    p_HalfClusterV.push_back( m_datacol.map_HalfCluster["HalfClusterColV"][ih].get() );

  return StatusCode::SUCCESS;
};


StatusCode TruthPatternRecAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

cout<<"  HFClusterU size: "<<p_HalfClusterU.size()<<endl;
  for(int ihc=0; ihc<p_HalfClusterU.size(); ihc++){
    std::map<edm4hep::MCParticle, std::vector<const Calo1DCluster*>> TruthAxesMap; 

    std::vector<const PandoraPlus::Calo1DCluster*> tmp_localMaxUCol = p_HalfClusterU[ihc]->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);
    for(int ilm=0; ilm<tmp_localMaxUCol.size(); ilm++){
      edm4hep::MCParticle mcp_lm = tmp_localMaxUCol[ilm]->getLeadingMCP();
      TruthAxesMap[mcp_lm].push_back(tmp_localMaxUCol[ilm]);
    }

cout<<"    In HFCluster #"<<ihc<<": print truth axes map"<<endl;
for(auto iter : TruthAxesMap)
  printf("      MCP id %d: localMax size %d \n", iter.first.getPDG(), iter.second.size());

    //Create axes from MCPMap. 
    std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> axisCol; axisCol.clear(); 
    for(auto& iter : TruthAxesMap){
      if(iter.second.size()==0) continue;
  
      std::shared_ptr<PandoraPlus::CaloHalfCluster> t_axis = std::make_shared<PandoraPlus::CaloHalfCluster>();
      for(int icl=0; icl<iter.second.size(); icl++) t_axis->addUnit(iter.second[icl]);
      t_axis->getLinkedMCPfromUnit();
      axisCol.push_back(t_axis);
      //m_datacol.map_HalfCluster["bkHalfCluster"].push_back(t_axis);
      //p_HalfClusterU[ihc]->addHalfCluster(settings.map_stringPars["OutputLongiClusName"], t_axis.get());    
    }
cout<<"    Axes size: "<<axisCol.size()<<endl;
    //Merge axes linked to the same MCP
    for(int iax=0; iax<axisCol.size(); iax++){
      for(int jax=iax+1; jax<axisCol.size(); jax++){
        if(axisCol[iax]->getLeadingMCP() == axisCol[jax]->getLeadingMCP()){
          axisCol[iax]->mergeHalfCluster( axisCol[jax].get() );
          axisCol.erase(axisCol.begin()+jax);
          jax--;
          if(iax>jax+1) iax--;
        }
      }
    }
cout<<"    Axes size after merge: "<<axisCol.size()<<endl;

    for(int iax=0; iax<axisCol.size(); iax++) p_HalfClusterU[ihc]->addHalfCluster(settings.map_stringPars["OutputLongiClusName"], axisCol[iax].get());
    m_datacol.map_HalfCluster["bkHalfCluster"].insert(m_datacol.map_HalfCluster["bkHalfCluster"].end(), axisCol.begin(), axisCol.end() );
  }


  //Get the linked axes in V
  for(int ihc=0; ihc<p_HalfClusterV.size(); ihc++){
    std::map<edm4hep::MCParticle, std::vector<const Calo1DCluster*>> TruthAxesMap; 

    std::vector<const PandoraPlus::Calo1DCluster*> tmp_localMaxVCol = p_HalfClusterV[ihc]->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);
    for(int ilm=0; ilm<tmp_localMaxVCol.size(); ilm++){
      edm4hep::MCParticle mcp_lm = tmp_localMaxVCol[ilm]->getLeadingMCP();
      TruthAxesMap[mcp_lm].push_back(tmp_localMaxVCol[ilm]);
    }

    //Create axes from MCPMap.
    std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> axisCol; axisCol.clear(); 
    for(auto& iter : TruthAxesMap){
      if(iter.second.size()==0) continue;
   
      std::shared_ptr<PandoraPlus::CaloHalfCluster> t_axis = std::make_shared<PandoraPlus::CaloHalfCluster>();
      for(int icl=0; icl<iter.second.size(); icl++) t_axis->addUnit(iter.second[icl]);
      t_axis->getLinkedMCPfromUnit();
      axisCol.push_back(t_axis);
      //m_datacol.map_HalfCluster["bkHalfCluster"].push_back(t_axis);
      //p_HalfClusterV[ihc]->addHalfCluster(settings.map_stringPars["OutputLongiClusName"], t_axis.get());
    }

    //Merge axes linked to the same MCP
    for(int iax=0; iax<axisCol.size(); iax++){
      for(int jax=iax+1; jax<axisCol.size(); jax++){
        if(axisCol[iax]->getLeadingMCP() == axisCol[jax]->getLeadingMCP()){
          axisCol[iax]->mergeHalfCluster( axisCol[jax].get() );
          axisCol.erase(axisCol.begin()+jax);
          jax--;
          if(iax>jax+1) iax--;
        }
      }
    }

    for(int iax=0; iax<axisCol.size(); iax++) p_HalfClusterV[ihc]->addHalfCluster(settings.map_stringPars["OutputLongiClusName"], axisCol[iax].get());
    m_datacol.map_HalfCluster["bkHalfCluster"].insert(m_datacol.map_HalfCluster["bkHalfCluster"].end(), axisCol.begin(), axisCol.end() );
  }

  return StatusCode::SUCCESS;
};


StatusCode TruthPatternRecAlg::ClearAlgorithm(){
  p_HalfClusterV.clear();
  p_HalfClusterU.clear();  

  return StatusCode::SUCCESS;
};


#endif
