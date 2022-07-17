#ifndef _PANDORAPLUS_DATA_H
#define _PANDORAPLUS_DATA_H
#include <iostream>
#include <algorithm>
#include <map>

#include "Objects/CaloUnit.h"
#include "Objects/Calo1DCluster.h"
#include "Objects/Calo2DCluster.h"
#include "Objects/Calo3DCluster.h"
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
  ~PandoraPlusDataCol() { Clean(); }
  StatusCode Clear(); 
  StatusCode Clean();

  //Readin CollectionMap
  std::map<std::string, std::vector<edm4hep::MCParticle> >       collectionMap_MC;
  std::map<std::string, std::vector<edm4hep::CalorimeterHit> >   collectionMap_CaloHit;
  std::map<std::string, std::vector<edm4hep::Vertex> >           collectionMap_Vertex;
  std::map<std::string, std::vector<edm4hep::Track> >            collectionMap_Track;
  std::map<std::string, std::vector<edm4hep::MCRecoCaloAssociation> > collectionMap_CaloRel;
  std::map<std::string, std::vector<edm4hep::MCRecoTrackerAssociation> > collectionMap_TrkRel;

  //Self used objects
  //General objects for all PFA
  std::vector<PandoraPlus::Track*>       TrackCol;
  std::map<std::string, std::vector<PandoraPlus::CaloHit*>> map_CaloHit;
  std::map<std::string, std::vector<PandoraPlus::CaloCluster*>> map_CaloCluster;


  std::vector<PandoraPlus::CaloUnit*>       BarCol; 
  std::vector<PandoraPlus::Calo1DCluster*>  Cluster1DCol; 
  std::vector<PandoraPlus::Calo2DCluster*>  Cluster2DCol;  
  std::vector<PandoraPlus::Calo3DCluster*>  Cluster3DCol;

  std::vector<PandoraPlus::CaloBlock*>   BlockCol; 
  std::vector<PandoraPlus::CaloTower*>   TowerCol;
  std::vector<PandoraPlus::TransShower*> TransShowerCol;
  //std::vector<PandoraPlus::CaloCluster*> ClusterCol;


  //Backup collections, for memory clean. TODO: replace with object managers. 
  std::vector<PandoraPlus::Track*>          bk_TrackCol;

  std::vector<PandoraPlus::CaloHit*>        bk_HitCol;
  std::vector<PandoraPlus::CaloUnit*>        bk_BarCol;
  std::vector<PandoraPlus::Calo1DCluster*>     bk_Cluster1DCol; 
  std::vector<PandoraPlus::Calo2DCluster*>     bk_Cluster2DCol;  
  std::vector<PandoraPlus::Calo3DCluster*>     bk_Cluster3DCol;
  std::vector<PandoraPlus::CaloBlock*>      bk_BlockCol;
  std::vector<PandoraPlus::CaloTower*>      bk_TowerCol;
  std::vector<PandoraPlus::CaloBarCluster*> bk_BarClusCol;
  std::vector<PandoraPlus::CaloBarShower*>  bk_BarShowerCol;
  std::vector<PandoraPlus::TransShower*>    bk_TransShowerCol;
  std::vector<PandoraPlus::LongiCluster*>   bk_LongiClusCol;
  std::vector<PandoraPlus::CaloCluster*>    bk_ClusterCol;


};
#endif
