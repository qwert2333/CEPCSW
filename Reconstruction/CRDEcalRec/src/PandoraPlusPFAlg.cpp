#ifndef PANDORAPLUS_ALG_C
#define PANDORAPLUS_ALG_C

#include "PandoraPlusPFAlg.h"


//#define C 299.79  // unit: mm/ns
//#define PI 3.141592653
using namespace std;
using namespace dd4hep;

int PandoraPlus::CaloUnit::Nmodule = 8;
int PandoraPlus::CaloUnit::Npart = 4;
int PandoraPlus::CaloUnit::Nstave = 11;
int PandoraPlus::CaloUnit::Nlayer = 14;
int PandoraPlus::CaloUnit::NbarPhi = 47;
int PandoraPlus::CaloUnit::NbarZ = 60;
int PandoraPlus::CaloUnit::over_module[28] = {13,15,16,18,19,21,22,24,25,26,28,29,30,32,33,35,36,38,39,41,42,43,45,46};
int PandoraPlus::CaloUnit::over_module_set = 2;

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
  m_algorithmManager.RegisterAlgorithmFactory("ConeClustering2DAlg",    new ConeClustering2DAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("AxisMergingAlg",         new AxisMergingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("EnergySplittingAlg",     new EnergySplittingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("EnergyTimeMatchingAlg",  new EnergyTimeMatchingAlg::Factory);
  //m_algorithmManager.RegisterAlgorithmFactory("ConeClusteringAlg",      new ConeClusteringAlg::Factory);
  //m_algorithmManager.RegisterAlgorithmFactory("ConeClusteringAlgHCAL",  new ConeClusteringAlg::Factory);


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
    t_Hough = new TTree("Hough","Hough");   // yyy: Hough
    t_LongiClus = new TTree("HalfClus", "HalfClus");
    t_Shower = new TTree("RecShowers", "RecShowers");
    t_Cluster = new TTree("RecClusters", "RecClusters");
    t_Clustering = new TTree("Clustering", "Clustering");
    t_Track = new TTree("RecTracks", "RecTracks");

    //Bar
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
    t_SimBar->Branch("simBar_bar", &m_simBar_bar);
   
    //BarShowers
    t_Layers->Branch("NshowerU", &m_NshowerU);
    t_Layers->Branch("NshowerV", &m_NshowerV);
    t_Layers->Branch("barShowerU_tag", &m_barShowerU_tag);
    t_Layers->Branch("barShowerU_x", &m_barShowerU_x);
    t_Layers->Branch("barShowerU_y", &m_barShowerU_y);
    t_Layers->Branch("barShowerU_z", &m_barShowerU_z);
    t_Layers->Branch("barShowerU_E", &m_barShowerU_E);
    //t_Layers->Branch("barShowerU_module", &m_barShowerU_module);
    //t_Layers->Branch("barShowerU_stave", &m_barShowerU_stave);
    //t_Layers->Branch("barShowerU_part", &m_barShowerU_part);
    //t_Layers->Branch("barShowerU_dlayer", &m_barShowerU_dlayer);
    //t_Layers->Branch("barShowerU_slayer", &m_barShowerU_slayer);
    //t_Layers->Branch("barShowerU_bar", &m_barShowerU_bar);
    t_Layers->Branch("barShowerV_tag", &m_barShowerV_tag);
    t_Layers->Branch("barShowerV_x", &m_barShowerV_x);
    t_Layers->Branch("barShowerV_y", &m_barShowerV_y);
    t_Layers->Branch("barShowerV_z", &m_barShowerV_z);
    t_Layers->Branch("barShowerV_E", &m_barShowerV_E);
    //t_Layers->Branch("barShowerV_module", &m_barShowerV_module);
    //t_Layers->Branch("barShowerV_stave", &m_barShowerV_stave);
    //t_Layers->Branch("barShowerV_part", &m_barShowerV_part);
    //t_Layers->Branch("barShowerV_dlayer", &m_barShowerV_dlayer);
    //t_Layers->Branch("barShowerV_slayer", &m_barShowerV_slayer);
    //t_Layers->Branch("barShowerV_bar", &m_barShowerV_bar);


    // yyy: Hough cluster
    t_Hough->Branch("NHoughAxisV", &m_NHoughAxisV);
    t_Hough->Branch("halfclusV_tag", &m_halfclusV_tag);
    t_Hough->Branch("houghV_x", &m_houghV_x);
    t_Hough->Branch("houghV_y", &m_houghV_y);
    t_Hough->Branch("houghV_z", &m_houghV_z);
    t_Hough->Branch("houghV_E", &m_houghV_E);
    //t_Hough->Branch("houghV_module", &m_houghV_module);
    //t_Hough->Branch("houghV_part", &m_houghV_part);
    //t_Hough->Branch("houghV_stave", &m_houghV_stave);
    //t_Hough->Branch("houghV_dlayer", &m_houghV_dlayer);
    //t_Hough->Branch("houghV_slayer", &m_houghV_slayer);
    t_Hough->Branch("NHoughAxisU", &m_NHoughAxisU);
    t_Hough->Branch("halfclusU_tag", &m_halfclusV_tag);
    t_Hough->Branch("houghU_x", &m_houghU_x);
    t_Hough->Branch("houghU_y", &m_houghU_y);
    t_Hough->Branch("houghU_z", &m_houghU_z);
    t_Hough->Branch("houghU_E", &m_houghU_E);
    //t_Hough->Branch("houghU_module", &m_houghU_module);
    //t_Hough->Branch("houghU_part", &m_houghU_part);
    //t_Hough->Branch("houghU_stave", &m_houghU_stave);
    //t_Hough->Branch("houghU_dlayer", &m_houghU_dlayer);
    //t_Hough->Branch("houghU_slayer", &m_houghU_slayer);
    t_Hough->Branch("NConeAxisV", &m_NConeAxisV);
    t_Hough->Branch("coneaxisV_tag", &m_coneaxisV_tag);
    t_Hough->Branch("coneaxisV_x", &m_coneaxisV_x);
    t_Hough->Branch("coneaxisV_y", &m_coneaxisV_y);
    t_Hough->Branch("coneaxisV_z", &m_coneaxisV_z);
    t_Hough->Branch("coneaxisV_E", &m_coneaxisV_E);
    t_Hough->Branch("NConeAxisU", &m_NConeAxisU);
    t_Hough->Branch("coneaxisU_tag", &m_coneaxisU_tag);
    t_Hough->Branch("coneaxisU_x", &m_coneaxisU_x);
    t_Hough->Branch("coneaxisU_y", &m_coneaxisU_y);
    t_Hough->Branch("coneaxisU_z", &m_coneaxisU_z);
    t_Hough->Branch("coneaxisU_E", &m_coneaxisU_E);
    t_Hough->Branch("NMergedAxisV", &m_NMergedAxisV);
    t_Hough->Branch("mergedaxisV_tag", &m_coneaxisV_tag);
    t_Hough->Branch("mergedaxisV_x",   &m_mergedaxisV_x);
    t_Hough->Branch("mergedaxisV_y",   &m_mergedaxisV_y);
    t_Hough->Branch("mergedaxisV_z",   &m_mergedaxisV_z);
    t_Hough->Branch("mergedaxisV_E",   &m_mergedaxisV_E);
    t_Hough->Branch("NMergedAxisU", &m_NMergedAxisU);
    t_Hough->Branch("mergedaxisU_tag", &m_mergedaxisU_tag);
    t_Hough->Branch("mergedaxisU_x",   &m_mergedaxisU_x);
    t_Hough->Branch("mergedaxisU_y",   &m_mergedaxisU_y);
    t_Hough->Branch("mergedaxisU_z",   &m_mergedaxisU_z);
    t_Hough->Branch("mergedaxisU_E",   &m_mergedaxisU_E);


    //CaloHalfClusters
    t_LongiClus->Branch("NHfClusU", &m_NHfClusU);
    t_LongiClus->Branch("NHfClusV", &m_NHfClusV);
    t_LongiClus->Branch("HfClusU_E", &m_HfClusU_E);
    t_LongiClus->Branch("HfClusU_Nhit", &m_HfClusU_Nhit);
    t_LongiClus->Branch("HfClusV_E", &m_HfClusV_E);
    t_LongiClus->Branch("HfClusV_Nhit", &m_HfClusV_Nhit);


    //2DShowers in each cluster
    t_Shower->Branch("Nshowers", &m_Nshowers);
    t_Shower->Branch("Eclus", &m_Eclus);
    t_Shower->Branch("shower2D_x", &m_shower2D_x);
    t_Shower->Branch("shower2D_y", &m_shower2D_y);
    t_Shower->Branch("shower2D_z", &m_shower2D_z);
    t_Shower->Branch("shower2D_E", &m_shower2D_E);
    t_Shower->Branch("shower2D_Module", &m_shower2D_Module);
    t_Shower->Branch("shower2D_Stave", &m_shower2D_Stave);
    t_Shower->Branch("shower2D_Part", &m_shower2D_Part);
    t_Shower->Branch("shower2D_Dlayer", &m_shower2D_Dlayer);


    //Clusters
    t_Cluster->Branch("Nclus", &m_Nclus);
    t_Cluster->Branch("Clus_x", &m_Clus_x);
    t_Cluster->Branch("Clus_y", &m_Clus_y);
    t_Cluster->Branch("Clus_z", &m_Clus_z);
    t_Cluster->Branch("Clus_E", &m_Clus_E);
    t_Cluster->Branch("Nhit", &m_Nhit);
    t_Cluster->Branch("Nmc", &m_Nmc);
    t_Cluster->Branch("mcPdgid",     &m_mcPdgid);
    t_Cluster->Branch("mcStatus",    &m_mcStatus);
    t_Cluster->Branch("mcPx", &m_mcPx);
    t_Cluster->Branch("mcPy", &m_mcPy);
    t_Cluster->Branch("mcPz", &m_mcPz);
    t_Cluster->Branch("mcEn", &m_mcEn);


    //neighbor clustering
    t_Clustering->Branch("m_3dcluster", &m_3dcluster);
    t_Clustering->Branch("m_2dcluster", &m_2dcluster);
    t_Clustering->Branch("m_1dcluster", &m_1dcluster);
    t_Clustering->Branch("m_barcluster", &m_barcluster);
    t_Clustering->Branch("m_bar", &m_bar);
    t_Clustering->Branch("m_units_2dcluster", &m_units_2dcluster);
    t_Clustering->Branch("m_slayer_2dcluster", &m_slayer_2dcluster);
    t_Clustering->Branch("m_E_3dcluster", &m_E_3dcluster);
    t_Clustering->Branch("m_E_2dcluster", &m_E_2dcluster);
    t_Clustering->Branch("m_E_1dcluster", &m_E_1dcluster);
    t_Clustering->Branch("m_E_barcluster", &m_E_barcluster);
    t_Clustering->Branch("m_E_bar", &m_E_bar);
    t_Clustering->Branch("m_bar_tag", &m_bar_tag);
    t_Clustering->Branch("m_bar_energy", &m_bar_energy);
    t_Clustering->Branch("m_bar_dlayer", &m_bar_dlayer);
    t_Clustering->Branch("m_bar_slayer", &m_bar_slayer);
    t_Clustering->Branch("m_bar_x", &m_bar_x);
    t_Clustering->Branch("m_bar_y", &m_bar_y);
    t_Clustering->Branch("m_bar_z", &m_bar_z);
    t_Clustering->Branch("m_bar_module", &m_bar_module);
    t_Clustering->Branch("m_bar_part", &m_bar_part);
    t_Clustering->Branch("m_bar_stave", &m_bar_stave);
    t_Clustering->Branch("m_bar_bar", &m_bar_bar);

    // Tracks
    t_Track->Branch("m_Ntrk", &m_Ntrk);
    t_Track->Branch("m_type", &m_type);
    t_Track->Branch("m_trkstate_d0", &m_trkstate_d0);
    t_Track->Branch("m_trkstate_z0", &m_trkstate_z0);
    t_Track->Branch("m_trkstate_phi", &m_trkstate_phi);
    t_Track->Branch("m_trkstate_tanL", &m_trkstate_tanL);
    t_Track->Branch("m_trkstate_kappa", &m_trkstate_kappa);
    t_Track->Branch("m_trkstate_omega", &m_trkstate_omega);
    t_Track->Branch("m_trkstate_refx", &m_trkstate_refx);
    t_Track->Branch("m_trkstate_refy", &m_trkstate_refy);
    t_Track->Branch("m_trkstate_refz", &m_trkstate_refz);
    t_Track->Branch("m_trkstate_location", &m_trkstate_location);
    t_Track->Branch("m_trkstate_tag", &m_trkstate_tag);

  }

  return GaudiAlgorithm::initialize();
}

StatusCode PandoraPlusPFAlg::execute()
{
  if(_nEvt==0) std::cout<<"PandoraPlusPFAlg::execute Start"<<std::endl;
  std::cout<<"Processing event: "<<_nEvt<<std::endl;

  if(_nEvt<m_Nskip){ _nEvt++;  return GaudiAlgorithm::initialize(); }

  //system("/cefs/higgs/songwz/winter22/CEPCSW/workarea/memory/memory_test.sh event_begin");
  //InitializeForNewEvent(); 
  m_DataCol.Clear();

  //Readin collections 
  std::cout<<"Readin MCParticle"<<std::endl;
  m_pMCParticleCreator->CreateMCParticle( m_DataCol, *r_MCParticleCol );
  cout<<"Readin Tracks"<<endl;
  m_pTrackCreator->CreateTracks( m_DataCol, r_TrackCols );
  cout<<"Readin CaloHits"<<endl;
  m_pCaloHitsCreator->CreateCaloHits( m_DataCol, r_CaloHitCols, map_readout_decoder );

  //system("/cefs/higgs/songwz/winter22/CEPCSW/workarea/memory/memory_test.sh read_data");
  //Perform PFA algorithm
  cout<<"Run Algorithms"<<endl;
  m_algorithmManager.RunAlgorithm( m_DataCol );
  //system("/cefs/higgs/songwz/winter22/CEPCSW/workarea/memory/memory_test.sh after_alg");

  m_pOutputCreator->CreateRecCaloHits( m_DataCol, w_RecCaloCol );
  m_pOutputCreator->CreateCluster( m_DataCol, w_ClusterCollection );
  m_pOutputCreator->CreatePFO( m_DataCol, w_ReconstructedParticleCollection );

  //Write Ana tuples
  cout<<"Writing tuples"<<endl;

cout<<"  Bar info"<<endl;
  //Save Raw bars information
  ClearBar();
  std::vector<PandoraPlus::CaloUnit*> m_barcol = m_DataCol.map_BarCol["BarCol"];
  for(int ibar=0;ibar<m_barcol.size();ibar++){
    const PandoraPlus::CaloUnit* p_hitbar = m_barcol[ibar];
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
    m_simBar_bar.push_back(p_hitbar->getBar());
  }
  t_SimBar->Fill();

cout<<"  LocalMax info"<<endl;
  //Save Layer info (local Max)
  ClearLayer();
  std::vector<PandoraPlus::CaloHalfCluster*> m_halfclusters = m_DataCol.map_HalfCluster["HalfClusterColU"];
cout<<"  HalfClusterU size: "<<m_halfclusters.size()<<endl;
  for(int ic=0;ic<m_halfclusters.size(); ic++){
    std::vector<const Calo1DCluster*> tmp_shower = m_halfclusters[ic]->getLocalMaxCol("AllLocalMax");
    m_NshowerU = tmp_shower.size(); 
    for(int is=0; is<tmp_shower.size(); is++)
    {
      m_barShowerU_tag.push_back( ic );
      m_barShowerU_x.push_back( tmp_shower[is]->getPos().x() );
      m_barShowerU_y.push_back( tmp_shower[is]->getPos().y() );
      m_barShowerU_z.push_back( tmp_shower[is]->getPos().z() );
      m_barShowerU_E.push_back( tmp_shower[is]->getEnergy() );
      //m_barShowerU_module.push_back( tmp_shower[is]->getSeeds().at(iseed)->getModule() );
      //m_barShowerU_stave.push_back( tmp_shower[is]->getSeeds().at(iseed)->getStave() );
      //m_barShowerU_part.push_back( tmp_shower[is]->getSeeds().at(iseed)->getPart() );
      //m_barShowerU_dlayer.push_back( tmp_shower[is]->getSeeds().at(iseed)->getDlayer() );
      //m_barShowerU_slayer.push_back( tmp_shower[is]->getSeeds().at(iseed)->getSlayer() );
      //m_barShowerU_bar.push_back( tmp_shower[is]->getSeeds().at(iseed)->getBar() );
      
    }
  }

  m_halfclusters.clear();
  m_halfclusters = m_DataCol.map_HalfCluster["HalfClusterColV"];
cout<<"  HalfClusterV size: "<<m_halfclusters.size();
  for(int ic=0;ic<m_halfclusters.size();ic++){
    std::vector<const Calo1DCluster*> tmp_shower = m_halfclusters[ic]->getLocalMaxCol("AllLocalMax");
    m_NshowerV = tmp_shower.size(); 
    for(int is=0; is<m_NshowerV; is++)
    {
        m_barShowerV_tag.push_back( ic );
        m_barShowerV_x.push_back( tmp_shower[is]->getPos().x() );
        m_barShowerV_y.push_back( tmp_shower[is]->getPos().y() );
        m_barShowerV_z.push_back( tmp_shower[is]->getPos().z() );
        m_barShowerV_E.push_back( tmp_shower[is]->getEnergy() );
        //m_barShowerV_module.push_back( tmp_shower[is]->getSeeds().at(iseed)->getModule() );
        //m_barShowerV_stave.push_back( tmp_shower[is]->getSeeds().at(iseed)->getStave() );
        //m_barShowerV_part.push_back( tmp_shower[is]->getSeeds().at(iseed)->getPart() );
        //m_barShowerV_dlayer.push_back( tmp_shower[is]->getSeeds().at(iseed)->getDlayer() );
        //m_barShowerV_slayer.push_back( tmp_shower[is]->getSeeds().at(iseed)->getSlayer() );
        //m_barShowerV_bar.push_back( tmp_shower[is]->getSeeds().at(iseed)->getBar() );
    }
  }
  t_Layers->Fill();

cout<<"  Hough axis info"<<endl;
  // yyy: check Hough cluster
  ClearHough();
  std::vector<PandoraPlus::CaloHalfCluster*> m_halfclusterV = m_DataCol.map_HalfCluster["HalfClusterColV"];
  std::vector<PandoraPlus::CaloHalfCluster*> m_halfclusterU = m_DataCol.map_HalfCluster["HalfClusterColU"];
  std::vector<const PandoraPlus::CaloHalfCluster*> m_houghaxisU; m_houghaxisU.clear(); 
  std::vector<const PandoraPlus::CaloHalfCluster*> m_houghaxisV; m_houghaxisV.clear(); 
  for(int i=0; i<m_halfclusterU.size(); i++){
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_clus = m_halfclusterU[i]->getHalfClusterCol("HoughAxis"); 
    m_houghaxisU.insert( m_houghaxisU.end(), tmp_clus.begin(), tmp_clus.end() );
  }
  for(int i=0; i<m_halfclusterV.size(); i++){
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_clus = m_halfclusterV[i]->getHalfClusterCol("HoughAxis");
    m_houghaxisV.insert( m_houghaxisV.end(), tmp_clus.begin(), tmp_clus.end() );
  }

  m_NHoughAxisU = m_houghaxisU.size();
  m_NHoughAxisV = m_houghaxisV.size();
/*  for(int ihc=0; ihc<m_houghaxisV.size(); ihc++){
    m_houghV_x.push_back( m_houghaxisV[ihc]->getPos().x() );
    m_houghV_y.push_back( m_houghaxisV[ihc]->getPos().y() );
    m_houghV_z.push_back( m_houghaxisV[ihc]->getPos().z() );
    m_houghV_E.push_back( m_houghaxisV[ihc]->getEnergy()  );
  }
  for(int ihc=0; ihc<m_houghaxisU.size(); ihc++){
    m_houghU_x.push_back( m_houghaxisU[ihc]->getPos().x() );
    m_houghU_y.push_back( m_houghaxisU[ihc]->getPos().y() );
    m_houghU_z.push_back( m_houghaxisU[ihc]->getPos().z() );
    m_houghU_E.push_back( m_houghaxisU[ihc]->getEnergy()  );
  }
*/
cout<<"  HoughAxisV size: "<<m_houghaxisV.size()<<endl;
  for(int ilc=0; ilc<m_houghaxisV.size(); ilc++){
    for(int ilm=0; ilm<m_houghaxisV[ilc]->getCluster().size(); ilm++){
      m_halfclusV_tag.push_back( m_houghaxisV[ilc]->getEnergy() );
      m_houghV_x.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getPos().x() );
      m_houghV_y.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getPos().y() );
      m_houghV_z.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getPos().z() );
      m_houghV_E.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getEnergy() );
      //m_houghV_module.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getTowerID()[0][0] );
      //m_houghV_part.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getTowerID()[0][1] );
      //m_houghV_stave.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getTowerID()[0][2] );
      //m_houghV_dlayer.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getDlayer() );
      //m_houghV_slayer.push_back( m_houghaxisV[ilc]->getCluster()[ilm]->getSlayer() );
    }
  }

cout<<"  HoughAxisU size: "<<m_houghaxisU.size()<<endl;
  for(int ilc=0; ilc<m_houghaxisU.size(); ilc++){
cout<<"    In Hough #"<<ilc<<": cluster size"<<m_houghaxisU[ilc]->getCluster().size()<<endl;
    for(int ilm=0; ilm<m_houghaxisU[ilc]->getCluster().size(); ilm++){
      m_halfclusU_tag.push_back( m_houghaxisU[ilc]->getEnergy() );
      m_houghU_x.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getPos().x() );
      m_houghU_y.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getPos().y() );
      m_houghU_z.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getPos().z() );
      m_houghU_E.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getEnergy() );
      //m_houghU_module.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getTowerID()[0][0] );
      //m_houghU_part.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getTowerID()[0][1] );
      //m_houghU_stave.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getTowerID()[0][2] );
      //m_houghU_dlayer.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getDlayer() );
      //m_houghU_slayer.push_back( m_houghaxisU[ilc]->getCluster()[ilm]->getSlayer() );
    }
  }

cout<<"  Cone axis info"<<endl;
  //Cone clustered axis
  std::vector<const PandoraPlus::CaloHalfCluster*> m_coneaxisU; m_coneaxisU.clear();
  std::vector<const PandoraPlus::CaloHalfCluster*> m_coneaxisV; m_coneaxisV.clear();
  for(int i=0; i<m_halfclusterU.size(); i++){
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_clus = m_halfclusterU[i]->getHalfClusterCol("ConeAxis");
    m_coneaxisU.insert( m_coneaxisU.end(), tmp_clus.begin(), tmp_clus.end() );
  }
  for(int i=0; i<m_halfclusterV.size(); i++){
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_clus = m_halfclusterV[i]->getHalfClusterCol("ConeAxis");
    m_coneaxisV.insert( m_coneaxisV.end(), tmp_clus.begin(), tmp_clus.end() );
  }  

  m_NConeAxisU = m_coneaxisU.size(); 
  m_NConeAxisV = m_coneaxisV.size(); 
  for(int ilc=0; ilc<m_coneaxisV.size(); ilc++){
    for(int ilm=0; ilm<m_coneaxisV[ilc]->getCluster().size(); ilm++){
      m_coneaxisV_tag.push_back( ilc );
      m_coneaxisV_x.push_back( m_coneaxisV[ilc]->getCluster()[ilm]->getPos().x() );
      m_coneaxisV_y.push_back( m_coneaxisV[ilc]->getCluster()[ilm]->getPos().y() );
      m_coneaxisV_z.push_back( m_coneaxisV[ilc]->getCluster()[ilm]->getPos().z() );
      m_coneaxisV_E.push_back( m_coneaxisV[ilc]->getCluster()[ilm]->getEnergy() );
    }
  }

  for(int ilc=0; ilc<m_coneaxisU.size(); ilc++){
    for(int ilm=0; ilm<m_coneaxisU[ilc]->getCluster().size(); ilm++){
      m_coneaxisU_tag.push_back( ilc );
      m_coneaxisU_x.push_back( m_coneaxisU[ilc]->getCluster()[ilm]->getPos().x() );
      m_coneaxisU_y.push_back( m_coneaxisU[ilc]->getCluster()[ilm]->getPos().y() );
      m_coneaxisU_z.push_back( m_coneaxisU[ilc]->getCluster()[ilm]->getPos().z() );
      m_coneaxisU_E.push_back( m_coneaxisU[ilc]->getCluster()[ilm]->getEnergy() );
    }
  }

  //Merged axes
  std::vector<const PandoraPlus::CaloHalfCluster*> m_mergedaxisU; m_mergedaxisU.clear();
  std::vector<const PandoraPlus::CaloHalfCluster*> m_mergedaxisV; m_mergedaxisV.clear();
  for(int i=0; i<m_halfclusterU.size(); i++){
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_clus = m_halfclusterU[i]->getHalfClusterCol("MergedAxis");
    m_mergedaxisU.insert( m_mergedaxisU.end(), tmp_clus.begin(), tmp_clus.end() );
  }
  for(int i=0; i<m_halfclusterV.size(); i++){
    std::vector<const PandoraPlus::CaloHalfCluster*> tmp_clus = m_halfclusterV[i]->getHalfClusterCol("MergedAxis");
    m_mergedaxisV.insert( m_mergedaxisV.end(), tmp_clus.begin(), tmp_clus.end() );
  }

  m_NMergedAxisU = m_mergedaxisU.size();
  m_NMergedAxisV = m_mergedaxisV.size();
  for(int ilc=0; ilc<m_mergedaxisV.size(); ilc++){
    for(int ilm=0; ilm<m_mergedaxisV[ilc]->getCluster().size(); ilm++){
      m_mergedaxisV_tag.push_back( ilc );
      m_mergedaxisV_x.push_back( m_mergedaxisV[ilc]->getCluster()[ilm]->getPos().x() );
      m_mergedaxisV_y.push_back( m_mergedaxisV[ilc]->getCluster()[ilm]->getPos().y() );
      m_mergedaxisV_z.push_back( m_mergedaxisV[ilc]->getCluster()[ilm]->getPos().z() );
      m_mergedaxisV_E.push_back( m_mergedaxisV[ilc]->getCluster()[ilm]->getEnergy() );
    }
  }

  for(int ilc=0; ilc<m_mergedaxisU.size(); ilc++){
    for(int ilm=0; ilm<m_mergedaxisU[ilc]->getCluster().size(); ilm++){
      m_mergedaxisU_tag.push_back( ilc );
      m_mergedaxisU_x.push_back( m_mergedaxisU[ilc]->getCluster()[ilm]->getPos().x() );
      m_mergedaxisU_y.push_back( m_mergedaxisU[ilc]->getCluster()[ilm]->getPos().y() );
      m_mergedaxisU_z.push_back( m_mergedaxisU[ilc]->getCluster()[ilm]->getPos().z() );
      m_mergedaxisU_E.push_back( m_mergedaxisU[ilc]->getCluster()[ilm]->getEnergy() );
    }
  }  
  
  t_Hough->Fill();



cout<<"  HalfCluster info"<<endl;
  //Half clusters after splitting.
  std::vector<PandoraPlus::CaloHalfCluster*> vec_LongiU = m_DataCol.map_HalfCluster["ESHalfClusterU"]; 
  std::vector<PandoraPlus::CaloHalfCluster*> vec_LongiV = m_DataCol.map_HalfCluster["ESHalfClusterV"]; 
  m_NHfClusU = vec_LongiU.size();
  m_NHfClusV = vec_LongiV.size();
  m_HfClusU_E.clear();
  m_HfClusU_Nhit.clear();
  m_HfClusV_E.clear();
  m_HfClusV_Nhit.clear();
  for(int ih=0; ih<m_NHfClusU; ih++){
    m_HfClusU_E.push_back(vec_LongiU[ih]->getEnergy());
    m_HfClusU_Nhit.push_back( vec_LongiU[ih]->getCluster().size() );
  }
  for(int ih=0; ih<m_NHfClusV; ih++){
    m_HfClusV_E.push_back(vec_LongiV[ih]->getEnergy());
    m_HfClusV_Nhit.push_back( vec_LongiV[ih]->getCluster().size() );
  }  
  t_LongiClus->Fill();


cout<<"  Cluster and MCtruth info"<<endl;
  //Save Cluster and MCP info
  ClearCluster();
  std::vector<PandoraPlus::Calo3DCluster*> m_clusvec = m_DataCol.map_CaloCluster["EcalCluster"];
  m_Nclus = m_clusvec.size();
  for(int ic=0; ic<m_clusvec.size(); ic++){
    m_Clus_x.push_back( m_clusvec[ic]->getShowerCenter().x() );
    m_Clus_y.push_back( m_clusvec[ic]->getShowerCenter().y() );
    m_Clus_z.push_back( m_clusvec[ic]->getShowerCenter().z() );
    m_Clus_E.push_back( m_clusvec[ic]->getEnergy() );
    m_Nhit.push_back( m_clusvec[ic]->getCluster().size() );
  }

  std::vector<edm4hep::MCParticle> m_MCPCol = m_DataCol.collectionMap_MC[name_MCParticleCol.value()];  
  m_Nmc = m_MCPCol.size(); 
  for(int imc=0; imc<m_MCPCol.size(); imc++){
    m_mcPdgid.push_back( m_MCPCol[imc].getPDG() );
    m_mcStatus.push_back( m_MCPCol[imc].getGeneratorStatus() );
    m_mcPx.push_back( m_MCPCol[imc].getMomentum()[0] );
    m_mcPy.push_back( m_MCPCol[imc].getMomentum()[1] );
    m_mcPz.push_back( m_MCPCol[imc].getMomentum()[2] );
    m_mcEn.push_back( m_MCPCol[imc].getEnergy() );
  }
  t_Cluster->Fill();

/*
  //Showers in cluster: 
  //ClearShower();
  for(int ic=0; ic<m_clusvec.size(); ic++){
    ClearShower();
    PandoraPlus::Calo3DCluster* p_cluster = m_clusvec[ic];

    m_Nshowers = p_cluster->getCluster().size();
    m_Eclus = p_cluster->getShowerE();
    for(int is=0; is<m_Nshowers; is++){
      m_shower2D_x.push_back( p_cluster->getCluster()[is]->getPos().x() );
      m_shower2D_y.push_back( p_cluster->getCluster()[is]->getPos().y() );
      m_shower2D_z.push_back( p_cluster->getCluster()[is]->getPos().z() );
      m_shower2D_E.push_back( p_cluster->getCluster()[is]->getShowerE() );
      m_shower2D_Module.push_back( p_cluster->getCluster()[is]->getModule() );
      m_shower2D_Part.push_back( p_cluster->getCluster()[is]->getPart() );
      m_shower2D_Stave.push_back( p_cluster->getCluster()[is]->getStave() );
      m_shower2D_Dlayer.push_back( p_cluster->getCluster()[is]->getDlayer() );
    }

    p_cluster = nullptr;
    t_Shower->Fill();
  }
  //t_Shower->Fill();

	//neighbor clustering
	ClearClustering();
	std::vector<PandoraPlus::Calo3DCluster*>  tmp_3dclusters = m_DataCol.Cluster3DCol;
	std::vector<PandoraPlus::CaloHalfCluster*>  tmp_halfclusters = m_DataCol.ClusterHalfCol;
	std::vector<PandoraPlus::Calo1DCluster*>  tmp_1dclusters = m_DataCol.Cluster1DCol;
	std::vector<PandoraPlus::CaloUnit*>  tmp_bars = m_DataCol.BarCol;
	
	m_3dcluster = tmp_3dclusters.size();
	m_2dcluster = tmp_halfclusters.size();
	m_1dcluster = tmp_1dclusters.size();
	m_bar = tmp_bars.size();

	for(int i=0; i<tmp_halfclusters.size(); i++)
	{
		PandoraPlus::CaloHalfCluster* tmp_halfcluster = tmp_halfclusters.at(i);
		m_E_2dcluster.push_back(tmp_halfcluster->getEnergy()); //
    m_units_2dcluster.push_back(tmp_halfcluster->getBars().size());
    m_slayer_2dcluster.push_back(tmp_halfcluster->getSlayer());
		std::vector<const PandoraPlus::CaloUnit*> allthebars = tmp_halfcluster->getBars(); //
		for(int n=0; n<allthebars.size(); n++)
		{
			const PandoraPlus::CaloUnit* the_bar = allthebars.at(n); //need to be fixed
			m_bar_tag.push_back(tmp_halfcluster->getEnergy());
			m_bar_energy.push_back(the_bar->getEnergy());
			m_bar_dlayer.push_back(the_bar->getDlayer());
			m_bar_slayer.push_back(the_bar->getSlayer());
			m_bar_x.push_back(the_bar->getPosition().x());
			m_bar_y.push_back(the_bar->getPosition().y());
			m_bar_z.push_back(the_bar->getPosition().z());
			m_bar_module.push_back(the_bar->getModule());
			m_bar_part.push_back(the_bar->getPart());
			m_bar_stave.push_back(the_bar->getStave());
      m_bar_bar.push_back(the_bar->getBar());
		}
	}
	for(int k=0; k<tmp_1dclusters.size(); k++)
	{
		PandoraPlus::Calo1DCluster* tmp_1dcluster = tmp_1dclusters.at(k);
		m_E_1dcluster.push_back(tmp_1dcluster->getEnergy());
	}
	for(int m=0; m<tmp_bars.size(); m++)
	{
		PandoraPlus::CaloUnit* tmp_bar = tmp_bars.at(m);
		m_E_bar.push_back(tmp_bar->getEnergy());
	}
	t_Clustering->Fill();
*/


  // Save Track info
  ClearTrack();
  std::vector<PandoraPlus::Track*> m_trkCol = m_DataCol.TrackCol;
  m_Ntrk = m_trkCol.size();
  for(int itrk=0; itrk<m_Ntrk; itrk++){
    m_type.push_back(m_trkCol[itrk]->getType());
    for(int istate=0; istate<m_trkCol[itrk]->trackStates_size(); istate++){
      m_trkstate_d0.push_back( m_trkCol[itrk]->getTrackStates(istate).D0 );
      m_trkstate_z0.push_back( m_trkCol[itrk]->getTrackStates(istate).Z0 );
      m_trkstate_phi.push_back( m_trkCol[itrk]->getTrackStates(istate).phi0 );
      m_trkstate_tanL.push_back( m_trkCol[itrk]->getTrackStates(istate).tanLambda );
      m_trkstate_kappa.push_back( m_trkCol[itrk]->getTrackStates(istate).Kappa);
      m_trkstate_omega.push_back( m_trkCol[itrk]->getTrackStates(istate).Omega );
      m_trkstate_refx.push_back( m_trkCol[itrk]->getTrackStates(istate).referencePoint.X() );
      m_trkstate_refy.push_back( m_trkCol[itrk]->getTrackStates(istate).referencePoint.Y() );
      m_trkstate_refz.push_back( m_trkCol[itrk]->getTrackStates(istate).referencePoint.Z() );
      m_trkstate_location.push_back( m_trkCol[itrk]->getTrackStates(istate).location );
      m_trkstate_tag.push_back(itrk);
  }}
  t_Track->Fill();


  //Clean Events
  //system("/cefs/higgs/songwz/winter22/CEPCSW/workarea/memory/memory_test.sh before_clean");
  m_DataCol.Clean();
  //system("/cefs/higgs/songwz/winter22/CEPCSW/workarea/memory/memory_test.sh event_end");

  std::cout<<"Event: "<<_nEvt<<" is done"<<std::endl;
  _nEvt ++ ;
  return StatusCode::SUCCESS;
}

StatusCode PandoraPlusPFAlg::finalize()
{
  //system("/cefs/higgs/songwz/winter22/CEPCSW/workarea/memory/memory_test.sh before_final");
  m_wfile->cd();
  t_SimBar->Write();
  t_Layers->Write();
  t_Hough->Write(); // yyy
  t_LongiClus->Write(); 
  t_Shower->Write();
  t_Cluster->Write();
  t_Clustering->Write();
  t_Track->Write();
  m_wfile->Close();
  delete m_wfile, t_SimBar, t_Layers, t_Hough, t_LongiClus, t_Cluster, t_Track, t_Clustering;


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
  //system("/cefs/higgs/songwz/winter22/CEPCSW/workarea/memory/memory_test.sh end_final");
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
  m_simBar_bar.clear();
}

void PandoraPlusPFAlg::ClearLayer(){
  m_NshowerU=-99;
  m_NshowerV=-99;
  m_barShowerU_tag.clear();
  m_barShowerU_x.clear();
  m_barShowerU_y.clear();
  m_barShowerU_z.clear();
  m_barShowerU_E.clear();
  m_barShowerU_module.clear();
  m_barShowerU_part.clear();
  m_barShowerU_stave.clear();
  m_barShowerU_dlayer.clear();
  m_barShowerU_slayer.clear();
  m_barShowerU_bar.clear();
  m_barShowerV_tag.clear();
  m_barShowerV_x.clear();
  m_barShowerV_y.clear();
  m_barShowerV_z.clear();
  m_barShowerV_E.clear();
  m_barShowerV_module.clear();
  m_barShowerV_part.clear();
  m_barShowerV_stave.clear();
  m_barShowerV_dlayer.clear();
  m_barShowerV_slayer.clear();
  m_barShowerV_bar.clear();
}

// yyy
void PandoraPlusPFAlg::ClearHough(){
  m_NHoughAxisU = -99;
  m_NHoughAxisV = -99;
  m_NConeAxisU = -99;
  m_NConeAxisV = -99;
  m_NMergedAxisU = -99;
  m_NMergedAxisV = -99;
  m_halfclusV_tag.clear();
  m_houghV_x.clear();
  m_houghV_y.clear();
  m_houghV_z.clear();
  m_houghV_E.clear();
  m_halfclusU_tag.clear();
  m_houghU_x.clear();
  m_houghU_y.clear();
  m_houghU_z.clear();
  m_houghU_E.clear();
  m_coneaxisV_tag.clear();
  m_coneaxisV_x.clear();
  m_coneaxisV_y.clear();
  m_coneaxisV_z.clear();
  m_coneaxisV_E.clear();
  m_coneaxisU_tag.clear();
  m_coneaxisU_x.clear();
  m_coneaxisU_y.clear();
  m_coneaxisU_z.clear();
  m_coneaxisU_E.clear();
  m_mergedaxisV_tag.clear();
  m_mergedaxisV_x.clear();
  m_mergedaxisV_y.clear();
  m_mergedaxisV_z.clear();
  m_mergedaxisV_E.clear();
  m_mergedaxisU_tag.clear();
  m_mergedaxisU_x.clear();
  m_mergedaxisU_y.clear();
  m_mergedaxisU_z.clear();
  m_mergedaxisU_E.clear();
}

void PandoraPlusPFAlg::ClearShower(){
  m_Nshowers = -99;
  m_Eclus = -99;
  m_shower2D_x.clear();
  m_shower2D_y.clear();
  m_shower2D_z.clear();
  m_shower2D_E.clear();
  m_shower2D_Module.clear();
  m_shower2D_Part.clear();
  m_shower2D_Stave.clear();
  m_shower2D_Dlayer.clear();
}

void PandoraPlusPFAlg::ClearCluster(){
  m_Nclus = -99;
  m_Nmc = -99;
  m_Clus_x.clear();
  m_Clus_y.clear();
  m_Clus_z.clear();
  m_Clus_E.clear();
  m_Nhit.clear();
  m_mcPdgid.clear();
  m_mcStatus.clear();
  m_mcPx.clear();
  m_mcPy.clear();
  m_mcPz.clear();
  m_mcEn.clear();

  //m_showerU_x.clear(); 
  //m_showerU_y.clear(); 
  //m_showerU_z.clear(); 
  //m_showerU_E.clear(); 
  //m_showerV_x.clear(); 
  //m_showerV_y.clear(); 
  //m_showerV_z.clear(); 
  //m_showerV_E.clear(); 
}


void PandoraPlusPFAlg::ClearClustering(){
	m_3dcluster=0;
	m_2dcluster=0;
	m_1dcluster=0;
	m_barcluster=0;
	m_bar=0;
  m_units_2dcluster.clear();
  m_slayer_2dcluster.clear();
	m_E_3dcluster.clear();
	m_E_2dcluster.clear();
	m_E_1dcluster.clear();
	m_E_barcluster.clear();
	m_E_bar.clear();
	m_bar_tag.clear();
	m_bar_energy.clear();
	m_bar_dlayer.clear();
	m_bar_slayer.clear(); 
	m_bar_x.clear();
	m_bar_y.clear();
	m_bar_z.clear();
  m_bar_module.clear();
  m_bar_part.clear();
  m_bar_stave.clear();
  m_bar_bar.clear();
}

void PandoraPlusPFAlg::ClearTrack(){
  m_Ntrk=-99;
  m_type.clear();
  m_trkstate_d0.clear(); 
  m_trkstate_z0.clear();
  m_trkstate_phi.clear();
  m_trkstate_tanL.clear();
  m_trkstate_kappa.clear();
  m_trkstate_omega.clear();
  m_trkstate_refx.clear();
  m_trkstate_refy.clear();
  m_trkstate_refz.clear();
  m_trkstate_location.clear();
  m_trkstate_tag.clear();
}
#endif
