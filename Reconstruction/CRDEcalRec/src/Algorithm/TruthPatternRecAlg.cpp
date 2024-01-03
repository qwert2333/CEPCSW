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

  if(settings.map_boolPars.find("DoAxisMerging")==settings.map_boolPars.end())
    settings.map_boolPars["DoAxisMerging"] = 1;
  if(settings.map_stringPars.find("ReadinAxisName")==settings.map_stringPars.end())
    settings.map_stringPars["ReadinAxisName"] = "TruthTrackAxis";
  if(settings.map_stringPars.find("OutputMergedAxisName")==settings.map_stringPars.end())
    settings.map_stringPars["OutputMergedAxisName"] = "TruthMergedAxis";
  if(settings.map_floatPars.find("th_overlap")==settings.map_floatPars.end()) 
    settings.map_floatPars["th_overlap"] = 0.8;

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
//cout<<"  TruthPatternRecAlg: Input HFCluster size ("<<p_HalfClusterU.size()<<", "<<p_HalfClusterV.size()<<") "<<endl;

  for(int ihc=0; ihc<p_HalfClusterU.size(); ihc++){
    std::map<edm4hep::MCParticle, std::vector<const Calo1DCluster*>> TruthAxesMap; 

    std::vector<const PandoraPlus::Calo1DCluster*> tmp_localMaxUCol = p_HalfClusterU[ihc]->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);
    for(int ilm=0; ilm<tmp_localMaxUCol.size(); ilm++){
      edm4hep::MCParticle mcp_lm = tmp_localMaxUCol[ilm]->getLeadingMCP();
      TruthAxesMap[mcp_lm].push_back(tmp_localMaxUCol[ilm]);
    }

//cout<<"    In HFClusterU #"<<ihc<<": localMax size "<<tmp_localMaxUCol.size()<<", print truth axes map"<<endl;
//for(auto iter : TruthAxesMap)
//  printf("      MCP id %d: localMax size %d \n", iter.first.getPDG(), iter.second.size());

    //Create axes from MCPMap. 
    std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> axisCol; axisCol.clear(); 
    for(auto& iter : TruthAxesMap){
      if(iter.second.size()==0) continue;
  
      std::shared_ptr<PandoraPlus::CaloHalfCluster> t_axis = std::make_shared<PandoraPlus::CaloHalfCluster>();
      for(int icl=0; icl<iter.second.size(); icl++) t_axis->addUnit(iter.second[icl]);
      t_axis->setType(100);
      t_axis->getLinkedMCPfromUnit();
      axisCol.push_back(t_axis);
    }
//cout<<"    Axes size: "<<axisCol.size()<<endl;

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

//cout<<"    In HFClusterV #"<<ihc<<": localMax size "<<tmp_localMaxVCol.size()<<", print truth axes map"<<endl;
//for(auto iter : TruthAxesMap)
//  printf("      MCP id %d: localMax size %d \n", iter.first.getPDG(), iter.second.size());

    //Create axes from MCPMap.
    std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> axisCol; axisCol.clear(); 
    for(auto& iter : TruthAxesMap){
      if(iter.second.size()==0) continue;
   
      std::shared_ptr<PandoraPlus::CaloHalfCluster> t_axis = std::make_shared<PandoraPlus::CaloHalfCluster>();
      for(int icl=0; icl<iter.second.size(); icl++) t_axis->addUnit(iter.second[icl]);
      t_axis->setType(100);
      t_axis->getLinkedMCPfromUnit();
      axisCol.push_back(t_axis);
    }

//cout<<"    Axes size: "<<axisCol.size()<<endl;

    for(int iax=0; iax<axisCol.size(); iax++) p_HalfClusterV[ihc]->addHalfCluster(settings.map_stringPars["OutputLongiClusName"], axisCol[iax].get());
    m_datacol.map_HalfCluster["bkHalfCluster"].insert(m_datacol.map_HalfCluster["bkHalfCluster"].end(), axisCol.begin(), axisCol.end() );
  }

  if(settings.map_boolPars["DoAxisMerging"]){

    for(int ihc=0; ihc<p_HalfClusterU.size(); ihc++){
      std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_newAxisCol; m_newAxisCol.clear();

      std::vector<const CaloHalfCluster*> m_recAxis = p_HalfClusterU[ihc]->getHalfClusterCol(settings.map_stringPars["OutputLongiClusName"]);
      std::vector<const CaloHalfCluster*> m_trkAxis = p_HalfClusterU[ihc]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"]);

      for(int iax=0; iax<m_recAxis.size(); iax++) m_newAxisCol.push_back(m_recAxis[iax]->Clone());
      for(int iax=0; iax<m_trkAxis.size(); iax++) m_newAxisCol.push_back(m_trkAxis[iax]->Clone());

//printf("  In HFClusterU #%d: rec axis %d, track axis %d. Print \n", ihc, m_recAxis.size(), m_trkAxis.size());
//for(int i=0; i<m_recAxis.size(); i++)
//  printf("    Rec axis %d: localMax size %d, track size %d, linked MC pid %d \n", i, m_recAxis[i]->getCluster().size(), m_recAxis[i]->getAssociatedTracks().size(), m_recAxis[i]->getLeadingMCP().getPDG());
//for(int i=0; i<m_trkAxis.size(); i++)
//  printf("    Trk axis %d: localMax size %d, track size %d, linked MC pid %d \n", i, m_trkAxis[i]->getCluster().size(), m_trkAxis[i]->getAssociatedTracks().size(), m_trkAxis[i]->getLeadingMCP().getPDG());

      if(m_newAxisCol.size()>1) OverlapMerging(m_newAxisCol);
//cout<<"  After merge: axis size "<<m_newAxisCol.size()<<endl;
//for(int i=0; i<m_newAxisCol.size(); i++)
//  printf("    Trk axis %d: localMax size %d, track size %d, linked MC pid %d \n", i, m_newAxisCol[i]->getCluster().size(), m_newAxisCol[i]->getAssociatedTracks().size(), m_newAxisCol[i]->getLeadingMCP().getPDG());

      for(int iax=0; iax<m_newAxisCol.size(); iax++) p_HalfClusterU[ihc]->addHalfCluster(settings.map_stringPars["OutputMergedAxisName"], m_newAxisCol[iax].get());
      m_datacol.map_HalfCluster["bkHalfCluster"].insert(m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_newAxisCol.begin(), m_newAxisCol.end());
    }

    for(int ihc=0; ihc<p_HalfClusterV.size(); ihc++){
      std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_newAxisCol; m_newAxisCol.clear();

      std::vector<const CaloHalfCluster*> m_recAxis = p_HalfClusterV[ihc]->getHalfClusterCol(settings.map_stringPars["OutputLongiClusName"]);
      std::vector<const CaloHalfCluster*> m_trkAxis = p_HalfClusterV[ihc]->getHalfClusterCol(settings.map_stringPars["ReadinAxisName"]);

      for(int iax=0; iax<m_recAxis.size(); iax++) m_newAxisCol.push_back(m_recAxis[iax]->Clone());
      for(int iax=0; iax<m_trkAxis.size(); iax++) m_newAxisCol.push_back(m_trkAxis[iax]->Clone());


//printf("  In HFClusterV #%d: rec axis %d, track axis %d. Print \n", ihc, m_recAxis.size(), m_trkAxis.size());
//for(int i=0; i<m_recAxis.size(); i++)
//  printf("    Rec axis %d: localMax size %d, track size %d, linked MC pid %d \n", i, m_recAxis[i]->getCluster().size(), m_recAxis[i]->getAssociatedTracks().size(), m_recAxis[i]->getLeadingMCP().getPDG());
//for(int i=0; i<m_trkAxis.size(); i++)
//  printf("    Trk axis %d: localMax size %d, track size %d, linked MC pid %d \n", i, m_trkAxis[i]->getCluster().size(), m_trkAxis[i]->getAssociatedTracks().size(), m_trkAxis[i]->getLeadingMCP().getPDG());

      if(m_newAxisCol.size()>1) OverlapMerging(m_newAxisCol);


//cout<<"  After merge: axis size "<<m_newAxisCol.size()<<endl;
//for(int i=0; i<m_newAxisCol.size(); i++)
//  printf("    Trk axis %d: localMax size %d, track size %d, linked MC pid %d \n", i, m_newAxisCol[i]->getCluster().size(), m_newAxisCol[i]->getAssociatedTracks().size(), m_newAxisCol[i]->getLeadingMCP().getPDG());
      for(int iax=0; iax<m_newAxisCol.size(); iax++) p_HalfClusterV[ihc]->addHalfCluster(settings.map_stringPars["OutputMergedAxisName"], m_newAxisCol[iax].get());
      m_datacol.map_HalfCluster["bkHalfCluster"].insert(m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_newAxisCol.begin(), m_newAxisCol.end());
    }

  }


  return StatusCode::SUCCESS;
};


StatusCode TruthPatternRecAlg::ClearAlgorithm(){
  p_HalfClusterV.clear();
  p_HalfClusterU.clear();  

  return StatusCode::SUCCESS;
};


StatusCode TruthPatternRecAlg::OverlapMerging( std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& m_axisCol ){
  if(m_axisCol.size()<2) return StatusCode::SUCCESS;

  for(int iax=0; iax<m_axisCol.size(); iax++){
    const PandoraPlus::CaloHalfCluster* m_axis = m_axisCol[iax].get();
    for(int jax=iax+1; jax<m_axisCol.size(); jax++){
      const PandoraPlus::CaloHalfCluster* p_axis = m_axisCol[jax].get();

      std::vector<const Calo1DCluster*> tmp_localMax = p_axis->getCluster();

      int nsharedHits = 0;
      for(int ihit=0; ihit<m_axis->getCluster().size(); ihit++)
        if( find(tmp_localMax.begin(), tmp_localMax.end(), m_axis->getCluster()[ihit])!=tmp_localMax.end() ) nsharedHits++;


      if( (m_axis->getCluster().size()<=p_axis->getCluster().size() && (float)nsharedHits/m_axis->getCluster().size()>settings.map_floatPars["th_overlap"] ) ||
          (p_axis->getCluster().size()<m_axis->getCluster().size() && (float)nsharedHits/p_axis->getCluster().size()>settings.map_floatPars["th_overlap"] ) ){

        m_axisCol[iax]->mergeHalfCluster( m_axisCol[jax].get() );

        // track axis + neutral axis = track axis
        int axis_type = m_axisCol[iax]->getType() + m_axisCol[jax]->getType();
        m_axisCol[iax]->setType(axis_type);

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
}

#endif
