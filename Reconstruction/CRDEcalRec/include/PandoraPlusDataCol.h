#ifndef _PANDORAPLUS_DATA_SVC_
#define _PANDORAPLUS_DATA_SVC_
#include <iostream>
#include <map>

#include "CRDEcalEDMSvc/CRDCaloBar.h"
#include "CRDEcalEDMSvc/CRDCaloBlock.h"
#include "CRDEcalEDMSvc/CRDCaloBarShower.h"
#include "CRDEcalEDMSvc/CRDCaloBarCluster.h"
#include "CRDEcalEDMSvc/CRDCaloLayer.h"
#include "CRDEcalEDMSvc/CRDShowerCandidate.h"
#include "CRDEcalEDMSvc/CRDCaloHit2DShower.h"
#include "CRDEcalEDMSvc/CRDCaloHit3DShower.h"
#include "CRDEcalEDMSvc/PFObject.h"
#include "CRDEcalEDMSvc/Track.h"

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
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> GoodClus3DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> BadClus3DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> Clus3DCol;


    //Temporary collections in iteration
    std::vector<CRDEcalEDM::CRDCaloBlock>       BlockVec_raw;

    std::vector<CRDEcalEDM::CRDCaloBlock>       BlockVec_iter0;
    std::vector<CRDEcalEDM::CRDCaloLayer>       LayerCol_iter0; 
    std::vector<CRDEcalEDM::CRDCaloHit2DShower> Shower2DCol_iter0;
    std::vector<CRDEcalEDM::CRDCaloHit3DShower> GoodClus3DCol_iter0;
    std::vector<CRDEcalEDM::CRDCaloHit3DShower> BadClus3DCol_iter0;
    std::vector<CRDEcalEDM::CRDCaloHit3DShower> Clus3DCol_iter0;

    std::vector<CRDEcalEDM::CRDCaloBlock>       BlockVec_iter1;
    std::vector<CRDEcalEDM::CRDCaloLayer>       LayerCol_iter1;
    std::vector<CRDEcalEDM::CRDCaloHit2DShower> Shower2DCol_iter1;
    std::vector<CRDEcalEDM::CRDCaloHit3DShower> GoodClus3DCol_iter1;
    std::vector<CRDEcalEDM::CRDCaloHit3DShower> BadClus3DCol_iter1;
    std::vector<CRDEcalEDM::CRDCaloHit3DShower> Clus3DCol_iter1;

    std::vector<CRDEcalEDM::CRDCaloBlock>       BlockVec_iter2;
    std::vector<CRDEcalEDM::CRDCaloLayer>       LayerCol_iter2;
    std::vector<CRDEcalEDM::CRDCaloHit2DShower> Shower2DCol_iter2;
    std::vector<CRDEcalEDM::CRDCaloHit3DShower> GoodClus3DCol_iter2;
    std::vector<CRDEcalEDM::CRDCaloHit3DShower> BadClus3DCol_iter2;
    std::vector<CRDEcalEDM::CRDCaloHit3DShower> Clus3DCol_iter2;


  //PFO Collections
  std::vector< CRDEcalEDM::PFObject >     PFOCol; 

  //Hcal Collections


  //Function
  void PrintLayer();
  void PrintShower();
  void Print3DClus();
  void Clear();
  void ClearBlock()  { BlockVec.clear(); }
  void ClearLayer()  { LayerCol.clear(); }
  void ClearShower() { Shower2DCol.clear(); }
  void ClearCluster(){ GoodClus3DCol.clear(); BadClus3DCol.clear(); Clus3DCol.clear(); }
  void ClearPFO()    { PFOCol.clear(); }
  void ClearTrack()  { TrackCol.clear(); }
  void ClearTempCol(); 


};
#endif
