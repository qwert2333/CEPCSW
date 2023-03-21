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
  m_algorithmManager.RegisterAlgorithmFactory("TrackMatchingAlg",       new TrackMatchingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("ConeClustering2DAlg",    new ConeClustering2DAlg::Factory);
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
    t_Match = new TTree("Match","Match");   // yyy: Match
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
    t_Layers->Branch("barShowerU_module", &m_barShowerU_module);
    t_Layers->Branch("barShowerU_stave", &m_barShowerU_stave);
    t_Layers->Branch("barShowerU_part", &m_barShowerU_part);
    t_Layers->Branch("barShowerU_dlayer", &m_barShowerU_dlayer);
    t_Layers->Branch("barShowerU_slayer", &m_barShowerU_slayer);
    t_Layers->Branch("barShowerU_bar", &m_barShowerU_bar);
    t_Layers->Branch("barShowerV_tag", &m_barShowerV_tag);
    t_Layers->Branch("barShowerV_x", &m_barShowerV_x);
    t_Layers->Branch("barShowerV_y", &m_barShowerV_y);
    t_Layers->Branch("barShowerV_z", &m_barShowerV_z);
    t_Layers->Branch("barShowerV_E", &m_barShowerV_E);
    t_Layers->Branch("barShowerV_module", &m_barShowerV_module);
    t_Layers->Branch("barShowerV_stave", &m_barShowerV_stave);
    t_Layers->Branch("barShowerV_part", &m_barShowerV_part);
    t_Layers->Branch("barShowerV_dlayer", &m_barShowerV_dlayer);
    t_Layers->Branch("barShowerV_slayer", &m_barShowerV_slayer);
    t_Layers->Branch("barShowerV_bar", &m_barShowerV_bar);


    // yyy: Hough cluster
    t_Hough->Branch("halfclusV_tag", &m_halfclusV_tag);
    t_Hough->Branch("longiclusV_tag", &m_longiclusV_tag);
    t_Hough->Branch("houghV_x", &m_houghV_x);
    t_Hough->Branch("houghV_y", &m_houghV_y);
    t_Hough->Branch("houghV_z", &m_houghV_z);
    t_Hough->Branch("houghV_E", &m_houghV_E);
    t_Hough->Branch("houghV_module", &m_houghV_module);
    t_Hough->Branch("houghV_part", &m_houghV_part);
    t_Hough->Branch("houghV_stave", &m_houghV_stave);
    t_Hough->Branch("houghV_dlayer", &m_houghV_dlayer);
    t_Hough->Branch("houghV_slayer", &m_houghV_slayer);
    t_Hough->Branch("halfclusU_tag", &m_halfclusU_tag);
    t_Hough->Branch("longiclusU_tag", &m_longiclusU_tag);
    t_Hough->Branch("houghU_x", &m_houghU_x);
    t_Hough->Branch("houghU_y", &m_houghU_y);
    t_Hough->Branch("houghU_z", &m_houghU_z);
    t_Hough->Branch("houghU_E", &m_houghU_E);
    t_Hough->Branch("houghU_module", &m_houghU_module);
    t_Hough->Branch("houghU_part", &m_houghU_part);
    t_Hough->Branch("houghU_stave", &m_houghU_stave);
    t_Hough->Branch("houghU_dlayer", &m_houghU_dlayer);
    t_Hough->Branch("houghU_slayer", &m_houghU_slayer);


    // yyy: Match
    t_Match->Branch("matchV_track_axis_tag", &matchV_track_axis_tag);
    t_Match->Branch("matchV_track_axis_x", &matchV_track_axis_x);
    t_Match->Branch("matchV_track_axis_y", &matchV_track_axis_y);
    t_Match->Branch("matchV_track_axis_z", &matchV_track_axis_z);
    t_Match->Branch("matchU_track_axis_tag", &matchU_track_axis_tag);
    t_Match->Branch("matchU_track_axis_x", &matchU_track_axis_x);
    t_Match->Branch("matchU_track_axis_y", &matchU_track_axis_y);
    t_Match->Branch("matchU_track_axis_z", &matchU_track_axis_z);
    


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

  //InitializeForNewEvent(); 
  m_DataCol.Clear();

  //Readin collections 
  std::cout<<"Readin MCParticle"<<std::endl;
  m_pMCParticleCreator->CreateMCParticle( m_DataCol, *r_MCParticleCol );
  cout<<"Readin Tracks"<<endl;
  m_pTrackCreator->CreateTracks( m_DataCol, r_TrackCols );
  cout<<"Readin CaloHits"<<endl;
  m_pCaloHitsCreator->CreateCaloHits( m_DataCol, r_CaloHitCols, map_readout_decoder );

  //Perform PFA algorithm
  cout<<"Run Algorithms"<<endl;
  m_algorithmManager.RunAlgorithm( m_DataCol );

  m_pOutputCreator->CreateRecCaloHits( m_DataCol, w_RecCaloCol );
  m_pOutputCreator->CreateCluster( m_DataCol, w_ClusterCollection );
  m_pOutputCreator->CreatePFO( m_DataCol, w_ReconstructedParticleCollection );

  //Write Ana tuples

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


  //Save Layer info (local Max)
  ClearLayer();
  std::vector<PandoraPlus::CaloHalfCluster*> m_halfclusters = m_DataCol.map_HalfCluster["HalfClusterColU"];
  for(int ic=0;ic<m_halfclusters.size();ic++){
    std::vector<const Calo1DCluster*> tmp_shower = m_halfclusters[ic]->getLocalMaxCol("AllLocalMax");
    m_NshowerU = tmp_shower.size(); 
    for(int is=0; is<m_NshowerU; is++)
    {
        m_barShowerU_tag.push_back( m_halfclusters[ic]->getEnergy() );
        m_barShowerU_x.push_back( tmp_shower[is]->getPos().x() );
        m_barShowerU_y.push_back( tmp_shower[is]->getPos().y() );
        m_barShowerU_z.push_back( tmp_shower[is]->getPos().z() );
        m_barShowerU_E.push_back( tmp_shower[is]->getEnergy() );
        m_barShowerU_module.push_back( tmp_shower[is]->getTowerID()[0][0] );
        m_barShowerU_stave.push_back( tmp_shower[is]->getTowerID()[0][2] );
        m_barShowerU_part.push_back( tmp_shower[is]->getTowerID()[0][1] );
        m_barShowerU_dlayer.push_back( tmp_shower[is]->getDlayer() );
        m_barShowerU_slayer.push_back( tmp_shower[is]->getSlayer() );
        // m_barShowerU_bar.push_back( tmp_shower[is]->getBar() ); 
    }
  }
  m_halfclusters.clear();
  m_halfclusters = m_DataCol.map_HalfCluster["HalfClusterColV"];
  for(int ic=0;ic<m_halfclusters.size();ic++){
    std::vector<const Calo1DCluster*> tmp_shower = m_halfclusters[ic]->getLocalMaxCol("AllLocalMax");
    m_NshowerV = tmp_shower.size(); 
    for(int is=0; is<m_NshowerV; is++)
    {
        m_barShowerV_tag.push_back( m_halfclusters[ic]->getEnergy() );
        m_barShowerV_x.push_back( tmp_shower[is]->getPos().x() );
        m_barShowerV_y.push_back( tmp_shower[is]->getPos().y() );
        m_barShowerV_z.push_back( tmp_shower[is]->getPos().z() );
        m_barShowerV_E.push_back( tmp_shower[is]->getEnergy() );
        m_barShowerV_module.push_back( tmp_shower[is]->getTowerID()[0][0] );
        m_barShowerV_stave.push_back( tmp_shower[is]->getTowerID()[0][2] );
        m_barShowerV_part.push_back( tmp_shower[is]->getTowerID()[0][1] );
        m_barShowerV_dlayer.push_back( tmp_shower[is]->getDlayer() );
        m_barShowerV_slayer.push_back( tmp_shower[is]->getSlayer() );
        // m_barShowerV_bar.push_back( tmp_shower[is]->getBar() );
    }
  }
  t_Layers->Fill();


  // yyy: check Hough cluster
  ClearHough();
  std::vector<PandoraPlus::CaloHalfCluster*> m_halfclusterV = m_DataCol.map_HalfCluster["HalfClusterColV"];
  std::vector<PandoraPlus::CaloHalfCluster*> m_halfclusterU = m_DataCol.map_HalfCluster["HalfClusterColU"];

  for(int i=0; i<m_halfclusterV.size(); i++){  // loop half cluster V
    std::vector<const PandoraPlus::CaloHalfCluster*> m_houghaxisV = m_halfclusterV[i]->getHalfClusterCol("HoughAxis"); 
    for(int iha=0; iha<m_houghaxisV.size(); iha++){ // loop Hough axis V
      for(int ilm=0; ilm<m_houghaxisV[iha]->getCluster().size(); ilm++){ // loop local max
        m_halfclusV_tag.push_back( m_halfclusterV[i]->getEnergy() );
        m_longiclusV_tag.push_back( m_houghaxisV[iha]->getEnergy() );
        m_houghV_x.push_back( m_houghaxisV[iha]->getCluster()[ilm]->getPos().x() );
        m_houghV_y.push_back( m_houghaxisV[iha]->getCluster()[ilm]->getPos().y() );
        m_houghV_z.push_back( m_houghaxisV[iha]->getCluster()[ilm]->getPos().z() );
        m_houghV_E.push_back( m_houghaxisV[iha]->getCluster()[ilm]->getEnergy() );
        m_houghV_module.push_back( m_houghaxisV[iha]->getCluster()[ilm]->getTowerID()[0][0] );
        m_houghV_part.push_back( m_houghaxisV[iha]->getCluster()[ilm]->getTowerID()[0][1] );
        m_houghV_stave.push_back( m_houghaxisV[iha]->getCluster()[ilm]->getTowerID()[0][2] );
        m_houghV_dlayer.push_back( m_houghaxisV[iha]->getCluster()[ilm]->getDlayer() );
        m_houghV_slayer.push_back( m_houghaxisV[iha]->getCluster()[ilm]->getSlayer() );
      }  // end loop local max
    }  // end loop Hough axis V
  }  // end loop half cluster V

  for(int i=0; i<m_halfclusterU.size(); i++){  // loop half cluster U
    std::vector<const PandoraPlus::CaloHalfCluster*> m_houghaxisU = m_halfclusterU[i]->getHalfClusterCol("HoughAxis"); 
    for(int iha=0; iha<m_houghaxisU.size(); iha++){ // loop Hough axis U
      for(int ilm=0; ilm<m_houghaxisU[iha]->getCluster().size(); ilm++){ // loop local max
        m_halfclusU_tag.push_back( m_halfclusterU[i]->getEnergy() );
        m_longiclusU_tag.push_back( m_houghaxisU[iha]->getEnergy() );
        m_houghU_x.push_back( m_houghaxisU[iha]->getCluster()[ilm]->getPos().x() );
        m_houghU_y.push_back( m_houghaxisU[iha]->getCluster()[ilm]->getPos().y() );
        m_houghU_z.push_back( m_houghaxisU[iha]->getCluster()[ilm]->getPos().z() );
        m_houghU_E.push_back( m_houghaxisU[iha]->getCluster()[ilm]->getEnergy() );
        m_houghU_module.push_back( m_houghaxisU[iha]->getCluster()[ilm]->getTowerID()[0][0] );
        m_houghU_part.push_back( m_houghaxisU[iha]->getCluster()[ilm]->getTowerID()[0][1] );
        m_houghU_stave.push_back( m_houghaxisU[iha]->getCluster()[ilm]->getTowerID()[0][2] );
        m_houghU_dlayer.push_back( m_houghaxisU[iha]->getCluster()[ilm]->getDlayer() );
        m_houghU_slayer.push_back( m_houghaxisU[iha]->getCluster()[ilm]->getSlayer() );
      }  // end loop local max
    }  // end loop Hough axis U
  }  // end loop half cluster U
  t_Hough->Fill();

  // yyy: TrackMatchingAlg
  ClearMatch();
  for(int i=0; i<m_halfclusterV.size(); i++){  // loop half cluster V
    std::vector<const PandoraPlus::CaloHalfCluster*> m_trackaxisV = m_halfclusterV[i]->getHalfClusterCol("TrackAxis");
    for(int ita=0; ita<m_trackaxisV.size(); ita++){ // loop track axis V
      for(int ilm=0; ilm<m_trackaxisV[ita]->getCluster().size(); ilm++){ // loop local max
        matchV_track_axis_tag.push_back( m_trackaxisV[ita]->getEnergy() );
        matchV_track_axis_x.push_back( m_trackaxisV[ita]->getCluster()[ilm]->getPos().x() );
        matchV_track_axis_y.push_back( m_trackaxisV[ita]->getCluster()[ilm]->getPos().y() );
        matchV_track_axis_z.push_back( m_trackaxisV[ita]->getCluster()[ilm]->getPos().z() );
      } 
    }
  }
  for(int i=0; i<m_halfclusterU.size(); i++){  // loop half cluster U
    std::vector<const PandoraPlus::CaloHalfCluster*> m_trackaxisU = m_halfclusterU[i]->getHalfClusterCol("TrackAxis");
    for(int ita=0; ita<m_trackaxisU.size(); ita++){ // loop track axis V
      for(int ilm=0; ilm<m_trackaxisU[ita]->getCluster().size(); ilm++){ // loop local max
        matchU_track_axis_tag.push_back( m_trackaxisU[ita]->getEnergy() );
        matchU_track_axis_x.push_back( m_trackaxisU[ita]->getCluster()[ilm]->getPos().x() );
        matchU_track_axis_y.push_back( m_trackaxisU[ita]->getCluster()[ilm]->getPos().y() );
        matchU_track_axis_z.push_back( m_trackaxisU[ita]->getCluster()[ilm]->getPos().z() );
      }
    }
  }
  t_Match->Fill();



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


/*
  //BarShower in each tower: 
  //ClearLayer();
  std::vector<PandoraPlus::Calo3DCluster*> m_3dclusters = m_DataCol.Cluster3DCol;
  for(int i=0; i<m_3dclusters.size(); i++){
  for(int itw=0; itw<m_3dclusters[i]->getTowers().size(); itw++){
  ClearLayer();
    const Calo3DCluster* p_tower = m_3dclusters[i]->getTowers()[itw];

    std::vector<const CaloHalfCluster*> m_longiClusU = p_tower->getCaloHalfClusterUCol("ESCaloHalfCluster");
    for(int il=0; il<m_longiClusU.size(); il++){
      for(int is=0; is<m_longiClusU[il]->getBarShowers().size(); is++){
        m_barShowerU_x.push_back( m_longiClusU[il]->getBarShowers()[is]->getPos().x() );
        m_barShowerU_y.push_back( m_longiClusU[il]->getBarShowers()[is]->getPos().y() );
        m_barShowerU_z.push_back( m_longiClusU[il]->getBarShowers()[is]->getPos().z() );
        m_barShowerU_E.push_back( m_longiClusU[il]->getBarShowers()[is]->getEnergy() );
      }
    }

    std::vector<const CaloHalfCluster*> m_longiClusV = p_tower->getCaloHalfClusterVCol("ESCaloHalfCluster");
    for(int il=0; il<m_longiClusV.size(); il++){
      for(int is=0; is<m_longiClusV[il]->getBarShowers().size(); is++){
        m_barShowerV_x.push_back( m_longiClusV[il]->getBarShowers()[is]->getPos().x() );
        m_barShowerV_y.push_back( m_longiClusV[il]->getBarShowers()[is]->getPos().y() );
        m_barShowerV_z.push_back( m_longiClusV[il]->getBarShowers()[is]->getPos().z() );
        m_barShowerV_E.push_back( m_longiClusV[il]->getBarShowers()[is]->getEnergy() );
      }
    }
  t_Layers->Fill();
  }}
  //t_Layers->Fill();
*/

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
*/


	//neighbor clustering
	ClearClustering();
	std::vector<PandoraPlus::CaloHalfCluster*>  tmp_halfclustersU = m_DataCol.map_HalfCluster["HalfClusterColU"];
  std::vector<PandoraPlus::CaloHalfCluster*>  tmp_halfclustersV = m_DataCol.map_HalfCluster["HalfClusterColV"];
  std::vector<PandoraPlus::CaloHalfCluster*>  tmp_halfclusters; tmp_halfclusters.clear();
  tmp_halfclusters.insert(tmp_halfclusters.end(), tmp_halfclustersU.begin(), tmp_halfclustersU.end());
  tmp_halfclusters.insert(tmp_halfclusters.end(), tmp_halfclustersV.begin(), tmp_halfclustersV.end());

	m_2dcluster = tmp_halfclusters.size();

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
	t_Clustering->Fill();



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
  m_wfile->cd();
  t_SimBar->Write();
  t_Layers->Write();
  t_Hough->Write(); // yyy
  t_Match->Write(); // yyy
  t_LongiClus->Write(); 
  t_Shower->Write();
  t_Cluster->Write();
  t_Clustering->Write();
  t_Track->Write();
  m_wfile->Close();
  delete m_wfile, t_SimBar, t_Layers, t_Hough, t_Match, t_LongiClus, t_Cluster, t_Track, t_Clustering;


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
  m_halfclusV_tag.clear();
  m_longiclusV_tag.clear();
  m_houghV_x.clear();
  m_houghV_y.clear();
  m_houghV_z.clear();
  m_houghV_E.clear();
  m_houghV_module.clear();
  m_houghV_part.clear();
  m_houghV_stave.clear();
  m_houghV_dlayer.clear();
  m_houghV_slayer.clear();
  m_halfclusU_tag.clear();
  m_longiclusU_tag.clear();
  m_houghU_x.clear();
  m_houghU_y.clear();
  m_houghU_z.clear();
  m_houghU_E.clear();
  m_houghU_module.clear();
  m_houghU_part.clear();
  m_houghU_stave.clear();
  m_houghU_dlayer.clear();
  m_houghU_slayer.clear();
}

// yyy
void PandoraPlusPFAlg::ClearMatch(){
  matchV_track_axis_tag.clear();
  matchV_track_axis_x.clear();
  matchV_track_axis_y.clear();
  matchV_track_axis_z.clear();
  matchU_track_axis_tag.clear();
  matchU_track_axis_x.clear();
  matchU_track_axis_y.clear();
  matchU_track_axis_z.clear();
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
