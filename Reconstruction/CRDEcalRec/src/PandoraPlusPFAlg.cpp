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
    t_Clustering = new TTree("Clustering", "Clustering");
    t_Track = new TTree("RecTracks", "RecTracks");


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
    t_Layers->Branch("NshowerU", &m_NshowerU);
    t_Layers->Branch("NshowerV", &m_NshowerV);
    t_Layers->Branch("barShowerU_x", &m_barShowerU_x);
    t_Layers->Branch("barShowerU_y", &m_barShowerU_y);
    t_Layers->Branch("barShowerU_z", &m_barShowerU_z);
    t_Layers->Branch("barShowerU_E", &m_barShowerU_E);
    t_Layers->Branch("barShowerV_x", &m_barShowerV_x);
    t_Layers->Branch("barShowerV_y", &m_barShowerV_y);
    t_Layers->Branch("barShowerV_z", &m_barShowerV_z);
    t_Layers->Branch("barShowerV_E", &m_barShowerV_E);

    //Clusters
    t_Cluster->Branch("Nclus", &m_Nclus);
    t_Cluster->Branch("Clus_x", &m_Clus_x);
    t_Cluster->Branch("Clus_y", &m_Clus_y);
    t_Cluster->Branch("Clus_z", &m_Clus_z);
    t_Cluster->Branch("Clus_E", &m_Clus_E);
    t_Cluster->Branch("Nhit", &m_Nhit);


    //neighbor clustering
    t_Clustering->Branch("m_3dcluster", &m_3dcluster);
    t_Clustering->Branch("m_2dcluster", &m_2dcluster);
    t_Clustering->Branch("m_1dcluster", &m_1dcluster);
    t_Clustering->Branch("m_barcluster", &m_barcluster);
    t_Clustering->Branch("m_bar", &m_bar);
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
  ClearBar();
  std::vector<PandoraPlus::CaloUnit*> m_barcol = m_DataCol.BarCol;
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
  }
  t_SimBar->Fill();

  //Save Layer info (local Max)
  ClearLayer();
  std::vector<PandoraPlus::Calo2DCluster*> m_blockvec = m_DataCol.Cluster2DCol; 
  for(int ibl=0;ibl<m_blockvec.size();ibl++){
    PandoraPlus::Calo2DCluster tmp_barcol = *(m_blockvec[ibl]);
    m_NshowerU = tmp_barcol.getShowerUCol().size();
    m_NshowerV = tmp_barcol.getShowerVCol().size();
    for(int is=0; is<m_NshowerU; is++){
      m_barShowerU_x.push_back( tmp_barcol.getShowerUCol()[is]->getPos().x() );
      m_barShowerU_y.push_back( tmp_barcol.getShowerUCol()[is]->getPos().y() );
      m_barShowerU_z.push_back( tmp_barcol.getShowerUCol()[is]->getPos().z() );
      m_barShowerU_E.push_back( tmp_barcol.getShowerUCol()[is]->getE() );
    }
    for(int is=0; is<m_NshowerV; is++){
      m_barShowerV_x.push_back( tmp_barcol.getShowerVCol()[is]->getPos().x() );
      m_barShowerV_y.push_back( tmp_barcol.getShowerVCol()[is]->getPos().y() );
      m_barShowerV_z.push_back( tmp_barcol.getShowerVCol()[is]->getPos().z() );
      m_barShowerV_E.push_back( tmp_barcol.getShowerVCol()[is]->getE() );
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


	//neighbor clustering
	ClearClustering();
	std::vector<PandoraPlus::Calo3DCluster*>  tmp_3dclusters = m_DataCol.Cluster3DCol;
	std::vector<PandoraPlus::Calo2DCluster*>  tmp_2dclusters = m_DataCol.Cluster2DCol;
	std::vector<PandoraPlus::Calo1DCluster*>  tmp_1dclusters = m_DataCol.Cluster1DCol;
	std::vector<PandoraPlus::CaloUnit*>  tmp_bars = m_DataCol.BarCol;
	
	m_3dcluster = tmp_3dclusters.size();
	m_2dcluster = tmp_2dclusters.size();
	m_1dcluster = tmp_1dclusters.size();
	m_bar = tmp_bars.size();
	for(int i=0; i<tmp_3dclusters.size(); i++)
	{
		PandoraPlus::Calo3DCluster* tmp_3dcluster = tmp_3dclusters.at(i);
		m_E_3dcluster.push_back(tmp_3dcluster->getEnergy()); //
		std::vector<const PandoraPlus::CaloUnit*> allthebars = tmp_3dcluster->getBars(); //
		for(int n=0; n<allthebars.size(); n++)
		{
			const PandoraPlus::CaloUnit* the_bar = allthebars.at(n); //need to be fixed
			m_bar_tag.push_back(tmp_3dcluster->getEnergy());
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
	for(int j=0; j<tmp_2dclusters.size(); j++)
	{
		PandoraPlus::Calo2DCluster* tmp_2dcluster = tmp_2dclusters.at(j);
		m_E_2dcluster.push_back(tmp_2dcluster->getEnergy());
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
  t_Clustering->Write();
  t_Track->Write();
  m_wfile->Close();
  delete m_wfile, t_SimBar, t_Layers, t_Cluster, t_Track, t_Clustering;


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
  m_NshowerU=-99;
  m_NshowerV=-99;
  m_barShowerU_x.clear();
  m_barShowerU_y.clear();
  m_barShowerU_z.clear();
  m_barShowerU_E.clear();
  m_barShowerV_x.clear();
  m_barShowerV_y.clear();
  m_barShowerV_z.clear();
  m_barShowerV_E.clear();
}

void PandoraPlusPFAlg::ClearCluster(){
  m_Nclus = -99;
  m_Clus_x.clear();
  m_Clus_y.clear();
  m_Clus_z.clear();
  m_Clus_E.clear();
  m_Nhit.clear();
}


void PandoraPlusPFAlg::ClearClustering(){
	m_3dcluster=0;
	m_2dcluster=0;
	m_1dcluster=0;
	m_barcluster=0;
	m_bar=0;
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
#endif
