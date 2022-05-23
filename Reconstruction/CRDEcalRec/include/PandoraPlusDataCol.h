#ifndef _PANDORAPLUS_DATA_SVC_
#define _PANDORAPLUS_DATA_SVC_
#include <iostream>
#include <algorithm>
#include <map>

#include "Objects/CRDCaloBar.h"
#include "Objects/CRDCaloBlock.h"
#include "Objects/CRDCaloTower.h"
#include "Objects/CRDCaloBarShower.h"
#include "Objects/CRDCaloBarCluster.h"
#include "Objects/CRDCaloLayer.h"
#include "Objects/CRDHoughObject.h"
#include "Objects/CRDHoughSpace.h"
#include "Objects/CRDShadowCluster.h"
#include "Objects/CRDCaloHitTransShower.h"
#include "Objects/CRDCaloHitLongiCluster.h"
#include "Objects/CRDCaloHit3DCluster.h"
#include "Objects/PFObject.h"
#include "Objects/Track.h"

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "Gaudi/Property.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/VertexCollection.h"
#include "edm4hep/MCRecoCaloAssociation.h"
#include "edm4hep/MCRecoTrackerAssociation.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/MCRecoParticleAssociationCollection.h"

#define PI 3.141592653
#define C 299.79  // unit: mm/ns

class PandoraPlusDataCol{
public:

  PandoraPlusDataCol() {}; 
  void Clear() {}; 

  //Readin CollectionMap
  std::map<std::string, std::vector<edm4hep::MCParticle> >                collectionMap_MC;
  std::map<std::string, std::vector<edm4hep::CalorimeterHit> >            collectionMap_CaloHit;
  std::map<std::string, std::vector<edm4hep::Vertex> >                    collectionMap_Vertex;
  std::map<std::string, std::vector<edm4hep::Track> >                     collectionMap_Track;
  std::map<std::string, std::vector<edm4hep::MCRecoCaloAssociation> >     collectionMap_CaloRel;
  std::map<std::string, std::vector<edm4hep::MCRecoTrackerAssociation> >  collectionMap_TrkRel;

  //Self used objects
  std::vector<CRDEcalEDM::CRDCaloBlock>   BlockVec; 
  std::vector<CRDEcalEDM::CRDCaloTower>   TowerCol;


};
#endif
