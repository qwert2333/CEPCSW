/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
// Unit in code: mm, ns. 
// NOTE: This digitialization highly matches detector geometry CRDEcalBarrel_v01.
// TODO: read geometry info automatically.  
#ifndef ARBORYOUTHALG_C
#define ARBORYOUTHALG_C
#include "ArborClusteringAlg.h"

#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Vector3f.h"
#include "edm4hep/Cluster.h"

#include "DD4hep/Detector.h"
#include <DD4hep/Objects.h>
#include <DDRec/CellIDPositionConverter.h>

#include "TVector3.h"
#include <math.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <map>

using namespace std;
using namespace dd4hep;

DECLARE_COMPONENT( ArborClusteringAlg )

ArborClusteringAlg::ArborClusteringAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{
  
	// Output collections
  declareProperty("ClusterCollection", w_ClusterCollection, "Handle of Cluster collection");
  declareProperty("RecoPFOCollection", w_ReconstructedParticleCollection, "Handle of Reconstructed PFO collection");
   
}

StatusCode ArborClusteringAlg::initialize()
{

  //Readin collections

  //---MC particle---
  if(!name_MCParticleCol.empty()) r_MCParticleCol = new DataHandle<edm4hep::MCParticleCollection> (name_MCParticleCol, Gaudi::DataHandle::Reader, this);

  //---Tracks---
  if(!name_TrackCol.empty()) r_TrackCols = new DataHandle<edm4hep::TrackCollection> (name_TrackCol, Gaudi::DataHandle::Reader, this);

  //---Calo Hits---
  if(!name_HcalHits.empty()) r_HCalHitCols = new DataHandle<edm4hep::CalorimeterHitCollection> (name_HcalHits, Gaudi::DataHandle::Reader, this);

  //Initialize services
  m_geosvc = service<IGeomSvc>("GeomSvc");
  if ( !m_geosvc )  throw "ArborClusteringAlg :Failed to find GeomSvc ...";

  m_decoder = m_geosvc->getDecoder(name_HcalReadout);
  if (m_decoder) {
    error() << "Failed to get the decoder for: " << name_HcalReadout << endmsg;
    return StatusCode::FAILURE;
  }

	rndm.SetSeed(_seed);

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

	std::cout<<"ArborClusteringAlg::initialize"<<std::endl;
	return GaudiAlgorithm::initialize();
}

StatusCode ArborClusteringAlg::execute()
{
	if(_nEvt==0) std::cout<<"ArborClusteringAlg::execute Start"<<std::endl;
	std::cout<<"Processing event: "<<_nEvt<<std::endl;
  if(_nEvt<_Nskip){ _nEvt++; return StatusCode::SUCCESS; }

  const edm4hep::TrackCollection* const_TrkCol = r_TrackCols->get();
  const edm4hep::CalorimeterHitCollection* const_HcalHitCol = r_HCalHitCols->get();

  //Convert CaloHit to ArborNodes
  std::vector<ArborNode> m_nodecol; m_nodecol.clear();
  for(int ihit=0; ihit<const_HcalHitCol->size(); ihit++){
    edm4hep::CalorimeterHit m_hit = const_HcalHitCol->at(ihit);
    ArborNode m_node(m_hit);
    m_nodecol.push_back(m_node);
  }

  //Build tree
  std::vector<ArborLink> m_fullLink; 
  BuildFullConnectionTree(m_nodecol, m_fullink);

  //Clean connection
  std::vector<ArborLink> m_cleanLink;
  LinkClean( m_nodecol, m_fullink, m_cleanLink ); 

  //Convert back to clusters

  //Match clusters and tracks
  MatchTracks()

  _nEvt ++ ;
  return StatusCode::SUCCESS;
}

StatusCode ArborClusteringAlg::finalize()
{
  m_wfile->cd();
  t_digiHit->Write();
  t_Track->Write();
  t_Cluster->Write();
  m_wfile->Close();
  delete m_wfile, t_digiHit, t_Cluster, t_Track;

  delete r_MCParticleCol;
  delete r_TrackCols;
  delete r_HCalHitCols;

  info() << "Processed " << _nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}

void ArborClusteringAlg::ClearHit(){
  m_digiHit_x.clear();
  m_digiHit_y.clear();
  m_digiHit_z.clear();
  m_digiHit_E.clear();
  m_digiHit_layer.clear();
}

void ArborClusteringAlg::ClearTrack(){
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

void ArborClusteringAlg::ClearCluster(){
  m_Nclus = -99;
  m_Clus_x.clear();
  m_Clus_y.clear();
  m_Clus_z.clear();
  m_Clus_E.clear();
  m_Clus_Nhit.clear();
}

#endif
