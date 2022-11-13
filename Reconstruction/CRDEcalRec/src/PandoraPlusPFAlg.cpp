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
  m_algorithmManager.RegisterAlgorithmFactory("EnergySplittingAlg",     new EnergySplittingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("EnergyTimeMatchingAlg",  new EnergyTimeMatchingAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("ConeClusteringAlg",      new ConeClusteringAlg::Factory);
  m_algorithmManager.RegisterAlgorithmFactory("ConeClusteringAlgHCAL",  new ConeClusteringAlg::Factory);


  //--- Create algorithm from readin settings ---
  for(int ialg=0; ialg<name_Algs.value().size(); ialg++){
    Settings m_settings; 
    for(int ipar=0; ipar<name_AlgPars.value()[ialg].size(); ipar++){
      if(type_AlgPars.value()[ialg].at(ipar)=="double")   m_settings.map_floatPars[  name_AlgPars.value()[ialg].at(ipar) ] = std::stod( (string)value_AlgPars.value()[ialg].at(ipar) );
      if(type_AlgPars.value()[ialg].at(ipar)=="string")   m_settings.map_stringPars[ name_AlgPars.value()[ialg].at(ipar) ] = value_AlgPars.value()[ialg].at(ipar) ;
      if(type_AlgPars.value()[ialg].at(ipar)=="bool")     m_settings.map_boolPars[   name_AlgPars.value()[ialg].at(ipar) ] = (bool)std::stoi( value_AlgPars.value()[ialg].at(ipar) );
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
    t_digiHit = new TTree("DigiHcalHit", "DigiHcalHit");
    t_Track = new TTree("RecTracks", "RecTracks");
    t_Cluster = new TTree("RecClusters", "RecClusters");

    //Hit
    t_digiHit->Branch("digiHit_x", &m_digiHit_x);
    t_digiHit->Branch("digiHit_y", &m_digiHit_y);
    t_digiHit->Branch("digiHit_z", &m_digiHit_z);
    t_digiHit->Branch("digiHit_E", &m_digiHit_E);
    t_digiHit->Branch("digiHit_layer", &m_digiHit_layer);

    //Track
    t_Track->Branch("m_Ntrk", &m_Ntrk);
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

    //Clusters
    t_Cluster->Branch("Nclus", &m_Nclus);
    t_Cluster->Branch("Clus_x", &m_Clus_x);
    t_Cluster->Branch("Clus_y", &m_Clus_y);
    t_Cluster->Branch("Clus_z", &m_Clus_z);
    t_Cluster->Branch("Clus_E", &m_Clus_E);
    t_Cluster->Branch("Nhit", &m_Clus_Nhit);

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


  m_pOutputCreator->CreateRecCaloHits( m_DataCol, w_RecCaloCol );
  m_pOutputCreator->CreateCluster( m_DataCol, w_ClusterCollection );
  m_pOutputCreator->CreatePFO( m_DataCol, w_ReconstructedParticleCollection );

  //Write Ana tuples

  //Save Raw bars information
  ClearHit();
  std::vector<PandoraPlus::CaloHit*> m_hitcol = m_DataCol.map_CaloHit["HCALBarrel"];
  std::vector<PandoraPlus::Calo2DCluster*> m_showerhit = m_DataCol.map_ShowerInLayer["EcalShowerInLayer"];
  for(int ihit=0; ihit<m_showerhit.size(); ihit++){
    const PandoraPlus::Calo2DCluster* p_hit = m_showerhit[ihit];
    m_digiHit_x.push_back(p_hit->getPos().x());
    m_digiHit_y.push_back(p_hit->getPos().y());
    m_digiHit_z.push_back(p_hit->getPos().z());
    m_digiHit_E.push_back(p_hit->getEnergy());
    m_digiHit_layer.push_back(p_hit->getDlayer());
  }
  for(int ihit=0; ihit<m_hitcol.size(); ihit++){
    const PandoraPlus::CaloHit* p_hit = m_hitcol[ihit];
    m_digiHit_x.push_back(p_hit->getPosition().x());
    m_digiHit_y.push_back(p_hit->getPosition().y());
    m_digiHit_z.push_back(p_hit->getPosition().z());
    m_digiHit_E.push_back(p_hit->getEnergy());
    m_digiHit_layer.push_back(p_hit->getLayer());
  }
  t_digiHit->Fill();


  // Save Track info
  ClearTrack();
  std::vector<PandoraPlus::Track*> m_trkCol = m_DataCol.TrackCol;
  m_Ntrk = m_trkCol.size();
  for(int itrk=0; itrk<m_Ntrk; itrk++){
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
  }}
  t_Track->Fill();


  //Save Cluster
  ClearCluster();
  std::vector<PandoraPlus::Calo3DCluster*> m_clusvec = m_DataCol.map_CaloCluster["HCALCluster"];
  m_Nclus = m_clusvec.size();
  for(int ic=0; ic<m_clusvec.size(); ic++){
    m_Clus_x.push_back( m_clusvec[ic]->getHitCenter().x() );
    m_Clus_y.push_back( m_clusvec[ic]->getHitCenter().y() );
    m_Clus_z.push_back( m_clusvec[ic]->getHitCenter().z() );
    m_Clus_E.push_back( m_clusvec[ic]->getHitsE() );
    m_Clus_Nhit.push_back( m_clusvec[ic]->getCaloHits().size() );
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
  t_digiHit->Write();
  t_Track->Write();
  t_Cluster->Write();
  m_wfile->Close();
  delete m_wfile, t_digiHit, t_Cluster, t_Track;


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

void PandoraPlusPFAlg::ClearHit(){
  m_digiHit_x.clear();
  m_digiHit_y.clear();
  m_digiHit_z.clear();
  m_digiHit_E.clear();
  m_digiHit_layer.clear();
}

void PandoraPlusPFAlg::ClearTrack(){
  m_Ntrk=-99;
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
}

void PandoraPlusPFAlg::ClearCluster(){
  m_Nclus = -99;
  m_Clus_x.clear();
  m_Clus_y.clear();
  m_Clus_z.clear();
  m_Clus_E.clear();
  m_Clus_Nhit.clear();
}

#endif
