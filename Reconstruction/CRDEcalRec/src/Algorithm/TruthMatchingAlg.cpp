#ifndef _TRUTHMATCHING_ALG_C
#define _TRUTHMATCHING_ALG_C

#include "Algorithm/TruthMatchingAlg.h"

StatusCode TruthMatchingAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_stringPars.find("ReadinHFClusterName")==settings.map_stringPars.end()) settings.map_stringPars["ReadinHFClusterName"] = "TruthESCluster";
  if(settings.map_stringPars.find("ReadinTowerName")==settings.map_stringPars.end()) settings.map_stringPars["ReadinTowerName"] = "TruthESTower";
  if(settings.map_stringPars.find("OutputClusterName")==settings.map_stringPars.end()) settings.map_stringPars["OutputClusterName"] = "TruthEcalCluster";
  
  return StatusCode::SUCCESS;
};

StatusCode TruthMatchingAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  std::cout<<"Initialize TruthMatchingAlg"<<std::endl;

  m_HFClusUCol.clear();
  m_HFClusVCol.clear();
  m_clusterCol.clear();
  m_towerCol.clear();
  m_bkCol.Clear();

  int ntower = m_datacol.map_CaloCluster[settings.map_stringPars["ReadinTowerName"]].size();
  for(int it=0; it<ntower; it++)
    m_towerCol.push_back( m_datacol.map_CaloCluster[settings.map_stringPars["ReadinTowerName"]][it].get() );

  return StatusCode::SUCCESS;
};

StatusCode TruthMatchingAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){
//cout<<"TruthMatchingAlg: readin tower size "<<m_towerCol.size()<<endl;

  for(int it=0; it<m_towerCol.size(); it++){

    std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> tmp_clusters; tmp_clusters.clear();
    m_HFClusUCol = m_towerCol.at(it)->getHalfClusterUCol(settings.map_stringPars["ReadinHFClusterName"]+"U");
    m_HFClusVCol = m_towerCol.at(it)->getHalfClusterVCol(settings.map_stringPars["ReadinHFClusterName"]+"V");
//cout<<"  In tower #"<<it<<": HFCluster size "<<m_HFClusUCol.size()<<", "<<m_HFClusVCol.size()<<endl;


    TruthMatching(m_HFClusUCol, m_HFClusVCol, tmp_clusters);
//cout<<"  After matching: 3DCluster size "<<tmp_clusters.size()<<endl;
    m_clusterCol.insert(m_clusterCol.end(), tmp_clusters.begin(), tmp_clusters.end());
  }

  //Merge clusters linked to the same MCP.
  for(int ic=0; ic<m_clusterCol.size() && m_clusterCol.size()>1; ic++){
    edm4hep::MCParticle mcp1 = m_clusterCol[ic].get()->getLeadingMCP();
    for(int jc=ic+1; jc<m_clusterCol.size(); jc++){
      if(ic>m_clusterCol.size()) ic--;
      edm4hep::MCParticle mcp2 = m_clusterCol[jc].get()->getLeadingMCP();
      if(mcp1==mcp2){
        m_clusterCol[ic].get()->mergeCluster( m_clusterCol[jc].get() );
        m_clusterCol.erase(m_clusterCol.begin()+jc);
        jc--;
        if(jc<ic) jc=ic;
      }
    }
  }
//cout<<"After cluster merging: size "<<m_clusterCol.size()<<endl;

  m_datacol.map_CaloCluster[settings.map_stringPars["OutputClusterName"]] = m_clusterCol;

  m_datacol.map_BarCol["bkBar"].insert( m_datacol.map_BarCol["bkBar"].end(), m_bkCol.map_BarCol["bkBar"].begin(), m_bkCol.map_BarCol["bkBar"].end() );
  m_datacol.map_1DCluster["bk1DCluster"].insert( m_datacol.map_1DCluster["bk1DCluster"].end(), m_bkCol.map_1DCluster["bk1DCluster"].begin(), m_bkCol.map_1DCluster["bk1DCluster"].end() );
  m_datacol.map_2DCluster["bk2DCluster"].insert( m_datacol.map_2DCluster["bk2DCluster"].end(), m_bkCol.map_2DCluster["bk2DCluster"].begin(), m_bkCol.map_2DCluster["bk2DCluster"].end() );
  m_datacol.map_HalfCluster["bkHalfCluster"].insert( m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_bkCol.map_HalfCluster["bkHalfCluster"].begin(), m_bkCol.map_HalfCluster["bkHalfCluster"].end() );



  return StatusCode::SUCCESS;
};

StatusCode TruthMatchingAlg::ClearAlgorithm(){
  std::cout<<"End run TruthMatchingAlg. Clean it."<<std::endl;

  m_HFClusUCol.clear();
  m_HFClusVCol.clear();
  m_clusterCol.clear();
  m_towerCol.clear();
  m_bkCol.Clear();

  return StatusCode::SUCCESS;
};


StatusCode TruthMatchingAlg::TruthMatching( std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClUCol,
                                            std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClVCol,
                                            std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>& m_clusters )
{
  if(m_ClUCol.size()==0 || m_ClVCol.size()==0) return StatusCode::SUCCESS;
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_truthESClUCol; m_truthESClUCol.clear();
  std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>> m_truthESClVCol; m_truthESClVCol.clear();

//cout<<"  Input HFCluster size "<<m_ClUCol.size()<<", "<<m_ClVCol.size()<<endl;

  //Truth split in HFClusterU
  for(int icl=0; icl<m_ClUCol.size(); icl++){
    auto truthMap = m_ClUCol[icl]->getLinkedMCP();
//cout<<"    In HFClusU #"<<icl<<": energy "<<m_ClUCol[icl]->getEnergy()<<", 1Dcluster size "<<m_ClUCol[icl]->getCluster().size()<<", truth link size "<<truthMap.size()<<endl;
//for(auto& iter: truthMap)
//  printf("    MC pid %d, weight %.3f \n", iter.first.getPDG(), iter.second);

    if(truthMap.size()==1){
      auto newClus = m_ClUCol[icl]->Clone();
      m_truthESClUCol.push_back(newClus);
      continue;
    }

    for(auto& iter: truthMap ){
      if(iter.second<0.05) continue;

      std::shared_ptr<PandoraPlus::CaloHalfCluster> newClus = std::make_shared<PandoraPlus::CaloHalfCluster>();
      //Copy and reweight 1D showers in HFCluster
      for(int ish=0; ish<m_ClUCol[icl]->getCluster().size(); ish++){
        const PandoraPlus::Calo1DCluster* p_shower = m_ClUCol[icl]->getCluster()[ish];

        std::vector<const PandoraPlus::CaloUnit*> Bars; Bars.clear();
        for(int ibar=0; ibar<p_shower->getCluster().size(); ibar++){
          auto bar = p_shower->getCluster()[ibar]->Clone();
          bar->setQ(bar->getQ1()*iter.second, bar->getQ2()*iter.second );

          bar->setLinkedMCP(std::vector<std::pair<edm4hep::MCParticle, float>>{std::make_pair(iter.first, 1.)});
          Bars.push_back(bar.get());
          m_bkCol.map_BarCol["bkBar"].push_back(bar);
        }

        std::shared_ptr<PandoraPlus::Calo1DCluster> shower = std::make_shared<PandoraPlus::Calo1DCluster>();
        shower->setBars(Bars);
        shower->setSeed();
        shower->setIDInfo();
        shower->getLinkedMCPfromUnit();
        newClus->addUnit(shower.get());
        m_bkCol.map_1DCluster["bk1DCluster"].push_back( shower );
      }
      for(int itrk=0; itrk<m_ClUCol[icl]->getAssociatedTracks().size(); itrk++)  newClus->addAssociatedTrack( m_ClUCol[icl]->getAssociatedTracks()[itrk] );
      newClus->setHoughPars( m_ClUCol[icl]->getHoughAlpha(), m_ClUCol[icl]->getHoughRho() );
      newClus->setIntercept( m_ClUCol[icl]->getHoughIntercept() );
      newClus->fitAxis("");
      newClus->getLinkedMCPfromUnit();
      newClus->mergeClusterInLayer();
      m_truthESClUCol.push_back(newClus);
    }
  }
//cout<<"    New HFClusterU size "<<m_truthESClUCol.size()<<endl;

  //Merge HFClusters linked to the same MCP
  for(int icl=0; icl<m_truthESClUCol.size() && m_truthESClUCol.size()>1; icl++){
    for(int jcl=icl+1; jcl<m_truthESClUCol.size(); jcl++){
      if(icl>m_truthESClUCol.size()) icl--;

      if(m_truthESClUCol[icl]->getLeadingMCP()==m_truthESClUCol[jcl]->getLeadingMCP()){
        m_truthESClUCol[icl].get()->mergeHalfCluster(m_truthESClUCol[jcl].get());
        m_truthESClUCol.erase(m_truthESClUCol.begin()+jcl);
        jcl--;
        if(jcl<icl) jcl=icl;
      }
    }
  }
//cout<<"    HFClusterU size after merge "<<m_truthESClUCol.size()<<endl;
  for(int icl=0; icl<m_truthESClUCol.size(); icl++) m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(m_truthESClUCol[icl]);


  //Truth split in HFClusterV
  for(int icl=0; icl<m_ClVCol.size(); icl++){
    auto truthMap = m_ClVCol[icl]->getLinkedMCP();
//cout<<"    In HFClusV #"<<icl<<": energy "<<m_ClVCol[icl]->getEnergy()<<", 1Dcluster size "<<m_ClVCol[icl]->getCluster().size()<<", truth link size "<<truthMap.size()<<endl;
//for(auto& iter: truthMap)
//  printf("    MC pid %d, weight %.3f \n", iter.first.getPDG(), iter.second);

    if(truthMap.size()==1){
      auto newClus = m_ClVCol[icl]->Clone();
      m_truthESClVCol.push_back(newClus);
      continue;
    }

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

          bar->setLinkedMCP(std::vector<std::pair<edm4hep::MCParticle, float>>{std::make_pair(iter.first, 1.)});
          Bars.push_back(bar.get());
          m_bkCol.map_BarCol["bkBar"].push_back(bar);
        }
        std::shared_ptr<PandoraPlus::Calo1DCluster> shower = std::make_shared<PandoraPlus::Calo1DCluster>();
        shower->setBars(Bars);
        shower->setSeed();
        shower->setIDInfo();
        shower->getLinkedMCPfromUnit();
        newClus->addUnit(shower.get());
        m_bkCol.map_1DCluster["bk1DCluster"].push_back( shower );
      }
      for(int itrk=0; itrk<m_ClVCol[icl]->getAssociatedTracks().size(); itrk++)  newClus->addAssociatedTrack( m_ClVCol[icl]->getAssociatedTracks()[itrk] );
      newClus->setHoughPars( m_ClVCol[icl]->getHoughAlpha(), m_ClVCol[icl]->getHoughRho() );
      newClus->setIntercept( m_ClVCol[icl]->getHoughIntercept() );
      newClus->fitAxis("");
      newClus->getLinkedMCPfromUnit();
      newClus->mergeClusterInLayer();
      m_truthESClVCol.push_back(newClus);
    }
  }
//cout<<"    New HFClusterV size "<<m_truthESClVCol.size()<<endl;

  //Merge HFClusters linked to the same MCP
  for(int icl=0; icl<m_truthESClVCol.size() && m_truthESClVCol.size()>1; icl++){
    for(int jcl=icl+1; jcl<m_truthESClVCol.size(); jcl++){
      if(icl>m_truthESClVCol.size()) icl--;

      if(m_truthESClVCol[icl]->getLeadingMCP()==m_truthESClVCol[jcl]->getLeadingMCP()){
        m_truthESClVCol[icl].get()->mergeHalfCluster(m_truthESClVCol[jcl].get());
        m_truthESClVCol.erase(m_truthESClVCol.begin()+jcl);
        jcl--;
        if(jcl<icl) jcl=icl;
      }
    }
  }
  for(int icl=0; icl<m_truthESClVCol.size(); icl++) m_bkCol.map_HalfCluster["bkHalfCluster"].push_back(m_truthESClVCol[icl]);
//cout<<"    HFClusterV size after merge "<<m_truthESClVCol.size()<<endl;


  //Doing matching.
  for(int icl=0; icl<m_truthESClUCol.size(); icl++){
    for(int jcl=0; jcl<m_truthESClVCol.size(); jcl++){
      if(m_truthESClUCol[icl]->getLeadingMCP() == m_truthESClVCol[jcl]->getLeadingMCP()){
        std::shared_ptr<PandoraPlus::Calo3DCluster> tmp_clus = std::make_shared<PandoraPlus::Calo3DCluster>();
        XYClusterMatchingL0(m_truthESClUCol[icl].get(), m_truthESClVCol[jcl].get(), tmp_clus);
        m_clusters.push_back(tmp_clus);
      }
    }
  }

  return StatusCode::SUCCESS;
}


//Longitudinal cluster: 1*1
StatusCode TruthMatchingAlg::XYClusterMatchingL0( const PandoraPlus::CaloHalfCluster* m_longiClU,
                                                  const PandoraPlus::CaloHalfCluster* m_longiClV,
                                                  std::shared_ptr<PandoraPlus::Calo3DCluster>& m_clus )
{
/*
cout<<"  XYClusterMatchingL0: print input HFClusterU. ";
printf("cluster size %d, linked MC pid %d, track size %d, Cover tower %d: ", m_longiClU->getCluster().size(), m_longiClU->getLeadingMCP().getPDG(), m_longiClU->getAssociatedTracks().size(), m_longiClU->getTowerID().size());
for(int it=0; it<m_longiClU->getTowerID().size(); it++) printf(" [%d, %d, %d] ", m_longiClU->getTowerID()[it][0], m_longiClU->getTowerID()[it][1], m_longiClU->getTowerID()[it][2]);
cout<<endl;
for(int ic=0; ic<m_longiClU->getCluster().size(); ic++){
  printf("    Cluster #%d: towerID [%d, %d, %d], layer %d, pos+E (%.3f, %.3f, %.3f, %.3f), bar size %d, seed size %d, MCP link size %d: ", 
    ic, m_longiClU->getCluster()[ic]->getTowerID()[0][0], m_longiClU->getCluster()[ic]->getTowerID()[0][1], m_longiClU->getCluster()[ic]->getTowerID()[0][2], m_longiClU->getCluster()[ic]->getDlayer(), 
    m_longiClU->getCluster()[ic]->getPos().x(), m_longiClU->getCluster()[ic]->getPos().y(), m_longiClU->getCluster()[ic]->getPos().z(),  m_longiClU->getCluster()[ic]->getEnergy(), 
    m_longiClU->getCluster()[ic]->getBars().size(), m_longiClU->getCluster()[ic]->getNseeds(), m_longiClU->getCluster()[ic]->getLinkedMCP().size() );
  for(int imc=0; imc<m_longiClU->getCluster()[ic]->getLinkedMCP().size(); imc++) cout<<m_longiClU->getCluster()[ic]->getLinkedMCP()[imc].first.getPDG()<<", ";
  cout<<endl;
}
cout<<"  XYClusterMatchingL0: print input HFClusterV. ";
printf("cluster size %d, linked MC pid %d, track size %d, Cover tower %d: ", m_longiClV->getCluster().size(), m_longiClV->getLeadingMCP().getPDG(),  m_longiClV->getAssociatedTracks().size(), m_longiClV->getTowerID().size());
for(int it=0; it<m_longiClV->getTowerID().size(); it++) printf(" [%d, %d, %d] ", m_longiClV->getTowerID()[it][0], m_longiClV->getTowerID()[it][1], m_longiClV->getTowerID()[it][2]);
cout<<endl;
for(int ic=0; ic<m_longiClV->getCluster().size(); ic++){
  printf("    Cluster #%d: towerID [%d, %d, %d], layer %d, pos+E (%.3f, %.3f, %.3f, %.3f), bar size %d, seed size %d, MCP link size %d: ", 
    ic, m_longiClV->getCluster()[ic]->getTowerID()[0][0], m_longiClV->getCluster()[ic]->getTowerID()[0][1], m_longiClV->getCluster()[ic]->getTowerID()[0][2], m_longiClV->getCluster()[ic]->getDlayer(), 
    m_longiClV->getCluster()[ic]->getPos().x(), m_longiClV->getCluster()[ic]->getPos().y(), m_longiClV->getCluster()[ic]->getPos().z(), m_longiClV->getCluster()[ic]->getEnergy(), 
    m_longiClV->getCluster()[ic]->getBars().size(), m_longiClV->getCluster()[ic]->getNseeds(), m_longiClV->getCluster()[ic]->getLinkedMCP().size() );
  for(int imc=0; imc<m_longiClV->getCluster()[ic]->getLinkedMCP().size(); imc++) cout<<m_longiClV->getCluster()[ic]->getLinkedMCP()[imc].first.getPDG()<<", ";
  cout<<endl;
}
cout<<endl;
*/

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

//cout<<"    in XYClusterMatchingL0: "<<layerindex.size()<<" layer need matching"<<endl;

  for(int il=0; il<layerindex.size(); il++){
    std::vector<const PandoraPlus::Calo1DCluster*> m_showerXcol = map_showersUinlayer[layerindex[il]];
    std::vector<const PandoraPlus::Calo1DCluster*> m_showerYcol = map_showersVinlayer[layerindex[il]];

    std::vector<PandoraPlus::Calo2DCluster*> m_showerinlayer; m_showerinlayer.clear();

//cout<<"    in layer "<<layerindex[il]<<": 1D shower size ("<<m_showerXcol.size()<<", "<<m_showerYcol.size()<<"). "<<endl;
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

    for(int is=0; is<m_showerinlayer.size(); is++) m_clus->addUnit(m_showerinlayer[is]);
//cout<<"    after matching: 2D shower size "<<m_showerinlayer.size()<<endl;
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


StatusCode TruthMatchingAlg::GetMatchedShowersL0( const PandoraPlus::Calo1DCluster* barShowerU,
                                                  const PandoraPlus::Calo1DCluster* barShowerV,
                                                  PandoraPlus::Calo2DCluster* outsh )
{
//cout<<"  In GetMatchedShowersL0"<<endl;

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


StatusCode TruthMatchingAlg::GetMatchedShowersL1( const PandoraPlus::Calo1DCluster* shower1,
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

#endif
