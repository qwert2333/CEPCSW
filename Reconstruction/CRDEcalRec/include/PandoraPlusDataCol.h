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
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackCollection.h"

#define PI 3.141592653
#define C 299.79  // unit: mm/ns

class PandoraPlusDataCol{
public:
 
  //MCParticle Collections
  std::vector<edm4hep::MCParticle>            MCParticleCol;

  //Track Collections
  //std::vector<edm4hep::Track>                 TrackCol; 
  std::vector<CRDEcalEDM::Track>              TrackCol; 

  //Vertex Collections

  //ECal Collections
  std::vector<CRDEcalEDM::CRDCaloBlock>       BlockVec;  //All fired crystal bars, grouped as blocks 
  std::vector<CRDEcalEDM::CRDCaloTower>       TowerCol;  
  std::vector<CRDEcalEDM::CRDCaloLayer>       LayerCol;  //Results of EnergySplittingAlg
  std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> LongiClusXCol; 
  std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> LongiClusYCol; 

    std::vector<CRDEcalEDM::CRDCaloLayer>       LayerCol_tmp;

  std::vector<CRDEcalEDM::CRDCaloHitTransShower> Shower2DCol;
  std::vector<CRDEcalEDM::CRDCaloHitTransShower> MIPShower2DCol; 
  std::vector<CRDEcalEDM::CRDCaloHitTransShower> EMShower2DCol; 
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> GoodClus3DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> BadClus3DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> Clus3DCol;


  //PFO Collections
  std::vector< CRDEcalEDM::PFObject >     PFOCol; 

  //Hcal Collections


  //Function
  void PrintTower(); 
  void PrintLayer();
  void PrintShower();
  void Print3DClus();
  void Clear();
  void ClearBlock()  { BlockVec.clear(); }
  void ClearTower()  { TowerCol.clear(); }
  void ClearLayer()  { LayerCol.clear(); }
  void ClearShower() { Shower2DCol.clear(); MIPShower2DCol.clear(); EMShower2DCol.clear(); }
  void ClearCluster(){ GoodClus3DCol.clear(); BadClus3DCol.clear(); Clus3DCol.clear(); }
  void ClearPFO()    { PFOCol.clear(); }
  void ClearTrack()  { TrackCol.clear(); }

};
#endif
