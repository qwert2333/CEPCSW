#ifndef _PANDORAPLUS_DATA_SVC_
#define _PANDORAPLUS_DATA_SVC_
#include <iostream>
#include <map>

#include "Objects/CRDCaloBar.h"
#include "Objects/CRDCaloBlock.h"
#include "Objects/CRDCaloBarShower.h"
#include "Objects/CRDCaloBarCluster.h"
#include "Objects/CRDCaloLayer.h"
#include "Objects/CRDShadowCluster.h"
#include "Objects/CRDCaloHit2DShower.h"
#include "Objects/CRDCaloHit3DCluster.h"
#include "Objects/CRDArborNode.h"
#include "Objects/CRDArborTree.h"
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
  std::vector<CRDEcalEDM::CRDCaloLayer>       LayerCol;  //Results of EnergySplittingAlg
  std::vector<CRDEcalEDM::CRDCaloHit2DShower> Shower2DCol;
  std::vector<CRDEcalEDM::CRDCaloHit2DShower> MIPShower2DCol; 
  std::vector<CRDEcalEDM::CRDCaloHit2DShower> EMShower2DCol; 
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> GoodClus3DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> BadClus3DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> Clus3DCol;

  //Another reconstruction way: ArborTrees
  std::vector<CRDEcalEDM::CRDArborTree>   ArborTreeCol;
  std::vector<CRDEcalEDM::CRDArborNode*>  IsoNodes;


    //Temporary collections in iteration
    std::vector<CRDEcalEDM::CRDCaloBlock>       BlockVec_raw;

    std::vector<CRDEcalEDM::CRDCaloBlock>       BlockVec_iter0;
    std::vector<CRDEcalEDM::CRDCaloLayer>       LayerCol_iter0; 
    std::vector<CRDEcalEDM::CRDCaloHit2DShower> Shower2DCol_iter0;
    std::vector<CRDEcalEDM::CRDCaloHit3DCluster> GoodClus3DCol_iter0;
    std::vector<CRDEcalEDM::CRDCaloHit3DCluster> BadClus3DCol_iter0;
    std::vector<CRDEcalEDM::CRDCaloHit3DCluster> Clus3DCol_iter0;

    std::vector<CRDEcalEDM::CRDCaloBlock>       BlockVec_iter1;
    std::vector<CRDEcalEDM::CRDCaloLayer>       LayerCol_iter1;
    std::vector<CRDEcalEDM::CRDCaloHit2DShower> Shower2DCol_iter1;
    std::vector<CRDEcalEDM::CRDCaloHit3DCluster> GoodClus3DCol_iter1;
    std::vector<CRDEcalEDM::CRDCaloHit3DCluster> BadClus3DCol_iter1;
    std::vector<CRDEcalEDM::CRDCaloHit3DCluster> Clus3DCol_iter1;

    std::vector<CRDEcalEDM::CRDCaloBlock>       BlockVec_iter2;
    std::vector<CRDEcalEDM::CRDCaloLayer>       LayerCol_iter2;
    std::vector<CRDEcalEDM::CRDCaloHit2DShower> Shower2DCol_iter2;
    std::vector<CRDEcalEDM::CRDCaloHit3DCluster> GoodClus3DCol_iter2;
    std::vector<CRDEcalEDM::CRDCaloHit3DCluster> BadClus3DCol_iter2;
    std::vector<CRDEcalEDM::CRDCaloHit3DCluster> Clus3DCol_iter2;


  //PFO Collections
  std::vector< CRDEcalEDM::PFObject >     PFOCol; 

  //Hcal Collections


  //Function
  void PrintLayer();
  void PrintShower();
  void Print3DClus();
  void PrintArborTree(); 
  void Clear();
  void ClearBlock()  { BlockVec.clear(); }
  void ClearLayer()  { LayerCol.clear(); }
  void ClearShower() { Shower2DCol.clear(); MIPShower2DCol.clear(); EMShower2DCol.clear(); }
  void ClearCluster(){ GoodClus3DCol.clear(); BadClus3DCol.clear(); Clus3DCol.clear(); }
  void ClearPFO()    { PFOCol.clear(); }
  void ClearTrack()  { TrackCol.clear(); }
  void ClearTempCol(); 
  void ClearArbor() { ArborTreeCol.clear(); IsoNodes.clear(); }

};
#endif
