#ifndef PANDORAPLUS_ALG_C
#define PANDORAPLUS_ALG_C

#include "PandoraPlusPFAlg.h"


//#define C 299.79  // unit: mm/ns
//#define PI 3.141592653
using namespace std;
using namespace dd4hep;


DECLARE_COMPONENT( PandoraPlusPFAlg )

PandoraPlusPFAlg::PandoraPlusPFAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{
 
  // Output collections
  declareProperty("RecCaloHitCollection", w_RecCaloCol, "Handle of Reconstructed CaloHit collection");
  declareProperty("ClusterCollection", w_ClusterCollection, "Handle of Cluster collection");   
  declareProperty("RecoPFOCollection", w_ReconstructedParticleCollection, "Handle of Reconstructed PFO collection");   
  declareProperty("RecoVtxCollection", w_VertexCollection, "Handle of Reconstructed vertex collection");   
  declareProperty("MCRecoPFOAssociationCollection", w_MCRecoParticleAssociationCollection, "Handle of MC-RecoPFO association collection");   

}

StatusCode PandoraPlusPFAlg::initialize()
{
  //Initialize global settings
  m_GlobalSettings.map_floatPars["BField"] = m_BField;
  m_GlobalSettings.map_floatPars["Seed"] = m_seed;
  m_GlobalSettings.map_intPars["Debug"] = m_Debug;
  m_GlobalSettings.map_intPars["SkipEvt"] = m_Nskip;


  //Initialize Creator settings
  m_pMCParticleCreatorSettings.map_stringPars["MCParticleCollections"] = name_MCParticleCol.value();

  m_pTrackCreatorSettings.map_stringVecPars["trackCollections"] = name_TrackCol.value();
  m_pTrackCreatorSettings.map_floatPars["BField"] = m_BField; 

  std::vector<std::string> name_CaloHits = name_EcalHits; 
  std::vector<std::string> name_CaloReadout = name_EcalReadout;
  name_CaloHits.insert( name_CaloHits.end(), name_HcalHits.begin(), name_HcalHits.end() );
  name_CaloReadout.insert(name_CaloReadout.end(), name_HcalReadout.begin(), name_HcalReadout.end());

  m_CaloHitsCreatorSettings.map_stringVecPars["CaloHitCollections"] = name_CaloHits;
  m_CaloHitsCreatorSettings.map_stringPars["EcalType"] = m_EcalType.value();

  //Initialize Creators
  m_pMCParticleCreator = new MCParticleCreator( m_pMCParticleCreatorSettings );
  m_pTrackCreator      = new TrackCreator( m_pTrackCreatorSettings );
  m_pCaloHitsCreator   = new CaloHitsCreator( m_CaloHitsCreatorSettings );
  m_pOutputCreator     = new OutputCreator( m_OutputCreatorSettings );

  //Readin collections

  //---MC particle---
  if(!name_MCParticleCol.empty()) r_MCParticleCol = new DataHandle<edm4hep::MCParticleCollection> (name_MCParticleCol, Gaudi::DataHandle::Reader, this);

  //---Tracks---
  for(auto& _trk : name_TrackCol) if(!_trk.empty()) r_TrackCols.push_back( new TrackType(_trk, Gaudi::DataHandle::Reader, this) );

  //---Calo Hits---
  for(auto& _ecal : name_EcalHits){
    if(!_ecal.empty()){ 
      r_ECalHitCols.push_back( new CaloType(_ecal, Gaudi::DataHandle::Reader, this) );
      r_CaloHitCols.push_back( new CaloType(_ecal, Gaudi::DataHandle::Reader, this) );
  }}
  for(auto& _hcal : name_HcalHits){ 
    if(!_hcal.empty()){
      r_HCalHitCols.push_back( new CaloType(_hcal, Gaudi::DataHandle::Reader, this) );
      r_CaloHitCols.push_back( new CaloType(_hcal, Gaudi::DataHandle::Reader, this) );
  }}


  //Register Algorithms
  //--- Initialize algorithm maps ---
  m_algorithmManager.RegisterAlgorithmFactory("ExampleAlg",             new ExampleAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("GlobalClusteringAlg",    new GlobalClusteringAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("LocalMaxFindingAlg",     new LocalMaxFindingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("HoughClusteringAlg",     new HoughClusteringAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("EnergySplittingAlg",     new EnergySplittingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("EnergyTimeMatchingAlg",  new EnergyTimeMatchingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("ConeClusteringAlg",      new ConeClusteringAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("ConeClusteringAlgHCAL",  new ConeClusteringAlg::Factory);


  //--- Create algorithm from readin settings ---
  for(int ialg=0; ialg<name_Algs.value().size(); ialg++){
    Settings m_settings; 
    for(int ipar=0; ipar<name_AlgPars.value()[ialg].size(); ipar++){
      if(type_AlgPars.value()[ialg].at(ipar)=="double") m_settings.map_floatPars[name_AlgPars.value()[ialg].at(ipar)] = std::stod( (string)value_AlgPars.value()[ialg].at(ipar) );
      if(type_AlgPars.value()[ialg].at(ipar)=="string") m_settings.map_stringPars[name_AlgPars.value()[ialg].at(ipar)] = value_AlgPars.value()[ialg].at(ipar) ;
      if(type_AlgPars.value()[ialg].at(ipar)=="bool") m_settings.map_boolPars[name_AlgPars.value()[ialg].at(ipar)] = (bool)std::stoi( value_AlgPars.value()[ialg].at(ipar) );
    }

    m_algorithmManager.RegisterAlgorithm( name_Algs.value()[ialg], m_settings );
  }


  //Initialize services
  m_geosvc = service<IGeomSvc>("GeomSvc");
  if ( !m_geosvc )  throw "PandoraPlusPFAlg :Failed to find GeomSvc ...";

  for(unsigned int i=0; i<name_CaloReadout.size(); i++){
    if(name_CaloReadout[i].empty()) continue;
    dd4hep::DDSegmentation::BitFieldCoder* tmp_decoder = m_geosvc->getDecoder(name_CaloReadout[i]);
    if (!tmp_decoder) {
      error() << "Failed to get the decoder for: " << name_CaloReadout[i] << endmsg;
      return StatusCode::FAILURE;
    }
    map_readout_decoder[name_CaloHits[i]] = tmp_decoder;
  }

  rndm.SetSeed(m_seed);
  std::cout<<"PandoraPlusPFAlg::initialize"<<std::endl;


  //Output collections
  
  //Output ntuple for analysis.
  if(m_WriteAna){
    std::string s_outfile = m_filename;
    m_wfile = new TFile(s_outfile.c_str(), "recreate");
    t_SimBar = new TTree("SimBarHit", "SimBarHit");
    t_Layers = new TTree("RecLayers","RecLayers");
    t_Cluster = new TTree("RecClusters", "RecClusters");

    t_SimBar->Branch("simBar_x", &m_simBar_x);
    t_SimBar->Branch("simBar_y", &m_simBar_y);
    t_SimBar->Branch("simBar_z", &m_simBar_z);
    t_SimBar->Branch("simBar_T1", &m_simBar_T1);
    t_SimBar->Branch("simBar_T2", &m_simBar_T2);
    t_SimBar->Branch("simBar_Q1", &m_simBar_Q1);
    t_SimBar->Branch("simBar_Q2", &m_simBar_Q2);
    t_SimBar->Branch("simBar_module", &m_simBar_module);
    t_SimBar->Branch("simBar_dlayer", &m_simBar_dlayer);
    t_SimBar->Branch("simBar_part", &m_simBar_part);
    t_SimBar->Branch("simBar_stave", &m_simBar_stave);
    t_SimBar->Branch("simBar_slayer", &m_simBar_slayer);
   
    //Local maxima
    t_Layers->Branch("NshowerX", &m_NshowerX);
    t_Layers->Branch("NshowerY", &m_NshowerY);
    t_Layers->Branch("barShowerX_x", &m_barShowerX_x);
    t_Layers->Branch("barShowerX_y", &m_barShowerX_y);
    t_Layers->Branch("barShowerX_z", &m_barShowerX_z);
    t_Layers->Branch("barShowerX_E", &m_barShowerX_E);
    t_Layers->Branch("barShowerY_x", &m_barShowerY_x);
    t_Layers->Branch("barShowerY_y", &m_barShowerY_y);
    t_Layers->Branch("barShowerY_z", &m_barShowerY_z);
    t_Layers->Branch("barShowerY_E", &m_barShowerY_E);

    //Clusters
    t_Cluster->Branch("Nclus", &m_Nclus);
    t_Cluster->Branch("Clus_x", &m_Clus_x);
    t_Cluster->Branch("Clus_y", &m_Clus_y);
    t_Cluster->Branch("Clus_z", &m_Clus_z);
    t_Cluster->Branch("Clus_E", &m_Clus_E);
    t_Cluster->Branch("Nhit", &m_Nhit);

  }

  return GaudiAlgorithm::initialize();
}

StatusCode PandoraPlusPFAlg::execute()
{
  if(_nEvt==0) std::cout<<"PandoraPlusPFAlg::execute Start"<<std::endl;
  std::cout<<"Processing event: "<<_nEvt<<std::endl;

  if(_nEvt<m_Nskip){ _nEvt++;  return GaudiAlgorithm::initialize(); }


  //InitializeForNewEvent(); 
  m_DataCol.Clear();

  //Readin collections 
  cout<<"Readin MCParticle"<<endl;
  m_pMCParticleCreator->CreateMCParticle( m_DataCol, *r_MCParticleCol );
  cout<<"Readin Tracks"<<endl;
  m_pTrackCreator->CreateTracks( m_DataCol, r_TrackCols );
  cout<<"Readin CaloHits"<<endl;
  m_pCaloHitsCreator->CreateCaloHits( m_DataCol, r_CaloHitCols, map_readout_decoder );

  //Perform PFA algorithm
  cout<<"Run Algorithms"<<endl;
  m_algorithmManager.RunAlgorithm( m_DataCol );

/*
cout<<"Block Info: "<<endl;
for(int ib=0; ib<m_DataCol.BlockCol.size(); ib++){
  printf("  #%d Block: ID=(%d, %d, %d, %d). Bar shower size: X=%d, Y=%d \n", ib, m_DataCol.BlockCol[ib]->getModule(), m_DataCol.BlockCol[ib]->getStave(), m_DataCol.BlockCol[ib]->getPart(), m_DataCol.BlockCol[ib]->getDlayer(), m_DataCol.BlockCol[ib]->getShowerXCol().size(), m_DataCol.BlockCol[ib]->getShowerYCol().size() );

  cout<<"  Print BarShowersX "<<endl;
  for(int is=0 ;is<m_DataCol.BlockCol[ib]->getShowerXCol().size(); is++)
  printf("    #%d showerX: cellID=(%d, %d, %d, %d) \n", is, m_DataCol.BlockCol[ib]->getShowerXCol()[is]->getModule(), m_DataCol.BlockCol[ib]->getShowerXCol()[is]->getStave(), m_DataCol.BlockCol[ib]->getShowerXCol()[is]->getPart(), m_DataCol.BlockCol[ib]->getShowerXCol()[is]->getDlayer() );

  cout<<"  Print BarShowersY "<<endl;
  for(int is=0 ;is<m_DataCol.BlockCol[ib]->getShowerYCol().size(); is++)
  printf("    #%d showerX: cellID=(%d, %d, %d, %d) \n", is, m_DataCol.BlockCol[ib]->getShowerYCol()[is]->getModule(), m_DataCol.BlockCol[ib]->getShowerYCol()[is]->getStave(), m_DataCol.BlockCol[ib]->getShowerYCol()[is]->getPart(), m_DataCol.BlockCol[ib]->getShowerYCol()[is]->getDlayer() );
}

cout<<endl;
cout<<endl;

cout<<"Tower Info: size = "<<m_DataCol.TowerCol.size()<<endl;
for(int it=0; it<m_DataCol.TowerCol.size(); it++){
  printf("  #%d Tower: ID=(%d, %d, %d). Block size: %d \n", it, m_DataCol.TowerCol[it]->getModule(), m_DataCol.TowerCol[it]->getStave(), m_DataCol.TowerCol[it]->getPart(), m_DataCol.TowerCol[it]->getBlocks().size() );

for(int ib=0; ib<m_DataCol.TowerCol[it]->getBlocks().size(); ib++){
  printf("    #%d Block: ID=(%d, %d, %d, %d). Bar shower size: X=%d, Y=%d \n", ib, m_DataCol.TowerCol[it]->getBlocks()[ib]->getModule(), m_DataCol.TowerCol[it]->getBlocks()[ib]->getStave(), m_DataCol.TowerCol[it]->getBlocks()[ib]->getPart(), m_DataCol.TowerCol[it]->getBlocks()[ib]->getDlayer(), m_DataCol.TowerCol[it]->getBlocks()[ib]->getShowerXCol().size(), m_DataCol.TowerCol[it]->getBlocks()[ib]->getShowerYCol().size() );

  cout<<"    Print BarShowersX "<<endl;
  for(int is=0 ;is<m_DataCol.TowerCol[it]->getBlocks()[ib]->getShowerXCol().size(); is++)
  printf("      #%d showerX: cellID=(%d, %d, %d, %d) \n", is, m_DataCol.TowerCol[it]->getBlocks()[ib]->getShowerXCol()[is]->getModule(), m_DataCol.TowerCol[it]->getBlocks()[ib]->getShowerXCol()[is]->getStave(), m_DataCol.TowerCol[it]->getBlocks()[ib]->getShowerXCol()[is]->getPart(), m_DataCol.TowerCol[it]->getBlocks()[ib]->getShowerXCol()[is]->getDlayer() );

  cout<<"    Print BarShowersY "<<endl;
  for(int is=0 ;is<m_DataCol.TowerCol[it]->getBlocks()[ib]->getShowerYCol().size(); is++)
  printf("      #%d showerX: cellID=(%d, %d, %d, %d) \n", is, m_DataCol.TowerCol[it]->getBlocks()[ib]->getShowerYCol()[is]->getModule(), m_DataCol.TowerCol[it]->getBlocks()[ib]->getShowerYCol()[is]->getStave(), m_DataCol.TowerCol[it]->getBlocks()[ib]->getShowerYCol()[is]->getPart(), m_DataCol.TowerCol[it]->getBlocks()[ib]->getShowerYCol()[is]->getDlayer() );

}
}
*/

//for(int ic=0; ic<m_DataCol.ClusterCol.size(); ic++) m_DataCol.ClusterCol[ic]->Print();


  m_pOutputCreator->CreateRecCaloHits( m_DataCol, w_RecCaloCol );
  m_pOutputCreator->CreateCluster( m_DataCol, w_ClusterCollection );
  m_pOutputCreator->CreatePFO( m_DataCol, w_ReconstructedParticleCollection );

  //Write Ana tuples

  //Save Raw bars information
  ClearBar();
  std::vector<PandoraPlus::CaloBlock*> m_blockvec = m_DataCol.BlockCol;
  for(int ibl=0;ibl<m_blockvec.size();ibl++){
    PandoraPlus::CaloBlock tmp_barcol = *(m_blockvec[ibl]);
    for(int ibar=0;ibar<tmp_barcol.getBarXCol().size();ibar++){
      const PandoraPlus::CaloBar* p_hitbar = tmp_barcol.getBarXCol()[ibar];
      m_simBar_x.push_back(p_hitbar->getPosition().x());
      m_simBar_y.push_back(p_hitbar->getPosition().y());
      m_simBar_z.push_back(p_hitbar->getPosition().z());
      m_simBar_Q1.push_back(p_hitbar->getQ1());
      m_simBar_Q2.push_back(p_hitbar->getQ2());
      m_simBar_T1.push_back(p_hitbar->getT1());
      m_simBar_T2.push_back(p_hitbar->getT2());
      m_simBar_module.push_back(p_hitbar->getModule());
      m_simBar_dlayer.push_back(p_hitbar->getDlayer());
      m_simBar_part.push_back(p_hitbar->getPart());
      m_simBar_stave.push_back(p_hitbar->getStave());
      m_simBar_slayer.push_back(p_hitbar->getSlayer());
    }
    for(int ibar=0;ibar<tmp_barcol.getBarYCol().size();ibar++){
      const PandoraPlus::CaloBar* p_hitbar = tmp_barcol.getBarYCol()[ibar];
      m_simBar_x.push_back(p_hitbar->getPosition().x());
      m_simBar_y.push_back(p_hitbar->getPosition().y());
      m_simBar_z.push_back(p_hitbar->getPosition().z());
      m_simBar_Q1.push_back(p_hitbar->getQ1());
      m_simBar_Q2.push_back(p_hitbar->getQ2());
      m_simBar_T1.push_back(p_hitbar->getT1());
      m_simBar_T2.push_back(p_hitbar->getT2());
      m_simBar_module.push_back(p_hitbar->getModule());
      m_simBar_dlayer.push_back(p_hitbar->getDlayer());
      m_simBar_part.push_back(p_hitbar->getPart());
      m_simBar_stave.push_back(p_hitbar->getStave());
      m_simBar_slayer.push_back(p_hitbar->getSlayer());
    }
  }
  t_SimBar->Fill();

  //Save Layer info (local Max)
  ClearLayer();
  for(int ibl=0;ibl<m_blockvec.size();ibl++){
    PandoraPlus::CaloBlock tmp_barcol = *(m_blockvec[ibl]);
    m_NshowerX = tmp_barcol.getShowerXCol().size();
    m_NshowerY = tmp_barcol.getShowerYCol().size();
    for(int is=0; is<m_NshowerX; is++){
      m_barShowerX_x.push_back( tmp_barcol.getShowerXCol()[is]->getPos().x() );
      m_barShowerX_y.push_back( tmp_barcol.getShowerXCol()[is]->getPos().y() );
      m_barShowerX_z.push_back( tmp_barcol.getShowerXCol()[is]->getPos().z() );
      m_barShowerX_E.push_back( tmp_barcol.getShowerXCol()[is]->getE() );
    }
    for(int is=0; is<m_NshowerY; is++){
      m_barShowerY_x.push_back( tmp_barcol.getShowerYCol()[is]->getPos().x() );
      m_barShowerY_y.push_back( tmp_barcol.getShowerYCol()[is]->getPos().y() );
      m_barShowerY_z.push_back( tmp_barcol.getShowerYCol()[is]->getPos().z() );
      m_barShowerY_E.push_back( tmp_barcol.getShowerYCol()[is]->getE() );
    }
  }
  t_Layers->Fill();

  //Save Cluster info
  ClearCluster();
  std::vector<PandoraPlus::CaloCluster*> m_clusvec = m_DataCol.map_CaloCluster["EcalCluster"];
  m_Nclus = m_clusvec.size();
  for(int ic=0; ic<m_clusvec.size(); ic++){
    m_Clus_x.push_back( m_clusvec[ic]->getShowerCenter().x() );
    m_Clus_y.push_back( m_clusvec[ic]->getShowerCenter().y() );
    m_Clus_z.push_back( m_clusvec[ic]->getShowerCenter().z() );
    m_Clus_E.push_back( m_clusvec[ic]->getShowerE() );
    m_Nhit.push_back( m_clusvec[ic]->getShowers().size() );
  }
  t_Cluster->Fill();

  //Clean Events
  m_DataCol.Clean();

  std::cout<<"Event: "<<_nEvt<<" is done"<<std::endl;
  _nEvt ++ ;
  return StatusCode::SUCCESS;
}

StatusCode PandoraPlusPFAlg::finalize()
{

  m_wfile->cd();
  t_SimBar->Write();
  t_Layers->Write();
  t_Cluster->Write();
  m_wfile->Close();
  delete m_wfile, t_SimBar, t_Layers, t_Cluster;

  delete m_pMCParticleCreator;
  delete m_pTrackCreator; 
  delete m_pCaloHitsCreator;
  delete m_pOutputCreator;

  delete r_MCParticleCol;
  for(auto iter : r_TrackCols) delete iter; 
  for(auto iter : r_ECalHitCols) delete iter; 
  for(auto iter : r_HCalHitCols) delete iter; 
  for(auto iter : r_CaloHitCols) delete iter; 
  r_TrackCols.clear();
  r_ECalHitCols.clear();
  r_HCalHitCols.clear();
  r_CaloHitCols.clear();


  info() << "Processed " << _nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}

void PandoraPlusPFAlg::ClearBar(){
  m_simBar_x.clear();
  m_simBar_y.clear();
  m_simBar_z.clear();
  m_simBar_T1.clear();
  m_simBar_T2.clear();
  m_simBar_Q1.clear();
  m_simBar_Q2.clear();
  m_simBar_module.clear();
  m_simBar_dlayer.clear();
  m_simBar_part.clear();
  m_simBar_stave.clear();
  m_simBar_slayer.clear();
}

void PandoraPlusPFAlg::ClearLayer(){
  m_NshowerX=-99;
  m_NshowerY=-99;
  m_barShowerX_x.clear();
  m_barShowerX_y.clear();
  m_barShowerX_z.clear();
  m_barShowerX_E.clear();
  m_barShowerY_x.clear();
  m_barShowerY_y.clear();
  m_barShowerY_z.clear();
  m_barShowerY_E.clear();
}

void PandoraPlusPFAlg::ClearCluster(){
  m_Nclus = -99;
  m_Clus_x.clear();
  m_Clus_y.clear();
  m_Clus_z.clear();
  m_Clus_E.clear();
  m_Nhit.clear();
}
#endif
