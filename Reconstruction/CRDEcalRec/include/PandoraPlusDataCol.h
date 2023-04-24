#ifndef _PANDORAPLUS_DATA_H
#define _PANDORAPLUS_DATA_H
#include <iostream>
#include <algorithm>
#include <map>

#include "Objects/CaloHit.h"
#include "Objects/CaloUnit.h"
#include "Objects/Calo1DCluster.h"
#include "Objects/Calo2DCluster.h"
#include "Objects/CaloHalfCluster.h"
#include "Objects/Calo3DCluster.h"
#include "Objects/HoughObject.h"
#include "Objects/HoughSpace.h"
//#include "Objects/PFObject.h"
#include "Objects/Track.h"

#include "k4FWCore/DataHandle.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/Vertex.h"
#include "edm4hep/VertexCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/MCRecoCaloAssociation.h"
#include "edm4hep/MCRecoTrackerAssociation.h"
#include "edm4hep/MCRecoParticleAssociationCollection.h"

#define PI 3.141592653
//#define C 299.79  // unit: mm/ns
using namespace std;
const double C = 299.79;
class PandoraPlusDataCol{
public:

  PandoraPlusDataCol() {}; 
  ~PandoraPlusDataCol() { Clear(); }
  StatusCode Clear(); 

  //Readin CollectionMap
  std::map<std::string, std::vector<edm4hep::MCParticle> >       collectionMap_MC;
  std::map<std::string, std::vector<edm4hep::CalorimeterHit> >   collectionMap_CaloHit;
  std::map<std::string, std::vector<edm4hep::Vertex> >           collectionMap_Vertex;
  std::map<std::string, std::vector<edm4hep::Track> >            collectionMap_Track;
  std::map<std::string, std::vector<edm4hep::MCRecoCaloAssociation> > collectionMap_CaloRel;
  std::map<std::string, std::vector<edm4hep::MCRecoTrackerAssociation> > collectionMap_TrkRel;

  //Self used objects
  //General objects for all PFA
  std::vector<std::shared_ptr<PandoraPlus::Track>>       TrackCol;
  std::map<std::string, std::vector<std::shared_ptr<PandoraPlus::CaloHit>>> map_CaloHit; //Hit
  std::map<std::string, std::vector<std::shared_ptr<PandoraPlus::CaloUnit>>> map_BarCol; 
  std::map<std::string, std::vector<std::shared_ptr<PandoraPlus::Calo1DCluster>>> map_1DCluster; 
  std::map<std::string, std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>> map_HalfCluster;
  std::map<std::string, std::vector<std::shared_ptr<PandoraPlus::Calo2DCluster>>> map_2DCluster; 
  std::map<std::string, std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>> map_CaloCluster; //Cluster


};
#endif
