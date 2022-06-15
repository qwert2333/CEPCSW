#ifndef _PANDORAPLUS_DATA_H
#define _PANDORAPLUS_DATA_H
#include <iostream>
#include <algorithm>
#include <map>

#include "Objects/CaloBar.h"
#include "Objects/CaloBlock.h"
#include "Objects/CaloTower.h"
#include "Objects/CaloBarShower.h"
#include "Objects/CaloBarCluster.h"
#include "Objects/HoughObject.h"
#include "Objects/HoughSpace.h"
#include "Objects/LongiCluster.h"
#include "Objects/TransShower.h"
#include "Objects/CaloCluster.h"
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
  std::vector<PandoraPlus::Track*>       TrackCol;

  std::vector<PandoraPlus::CaloBar*>     BarCol;  
  std::vector<PandoraPlus::CaloBlock*>   BlockCol; 
  std::vector<PandoraPlus::CaloTower*>   TowerCol;
  std::vector<PandoraPlus::TransShower*> TransShowerCol;
  std::vector<PandoraPlus::CaloCluster*> ClusterCol;

};
#endif
