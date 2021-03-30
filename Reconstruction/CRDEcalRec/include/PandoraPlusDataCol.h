#ifndef _PANDORAPLUS_DATA_SVC_
#define _PANDORAPLUS_DATA_SVC_
#include <iostream>
#include <map>

#include "CRDEcalEDMSvc/CRDCaloBar.h"
#include "CRDEcalEDMSvc/CRDCaloBlock.h"
#include "CRDEcalEDMSvc/CRDCaloBarShower.h"
#include "CRDEcalEDMSvc/CRDCaloBarCluster.h"
#include "CRDEcalEDMSvc/CRDCaloLayer.h"
#include "CRDEcalEDMSvc/CRDExpEMShower.h"
#include "CRDEcalEDMSvc/CRDCaloHit2DShower.h"
#include "CRDEcalEDMSvc/CRDCaloHit3DShower.h"
#include "CRDEcalEDMSvc/PFObject.h"

#include "k4FWCore/DataHandle.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCParticle.h"

#define PI 3.141592653
#define C 299.79  // unit: mm/ns

class PandoraPlusDataCol{
public:
 
  //MCParticle Collections
  std::vector<edm4hep::MCParticle>            MCParticleCol;

  //Track Collections
  
  //Vertex Collections

  //ECal Collections
  std::vector<CRDEcalEDM::CRDCaloBlock>       blockVec;  //All fired crystal bars, grouped as blocks 
  std::vector<CRDEcalEDM::CRDCaloLayer>       LayerCol;  //Results of EnergySplittingAlg
  std::vector<CRDEcalEDM::CRDCaloHit2DShower> shower2DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> GoodClus3DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> BadClus3DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> Clus3DCol;

  bool Flag_Iter;

  //PFO Collections
  std::vector< CRDEcalEDM::PFObject >     PFOCol; 

  //Hcal Collections


  //Function
  void PrintLayer();
  void PrintShower();
  void Print3DClus();
  void Clear();
  void ClearBlock()  { blockVec.clear(); }
  void ClearLayer()  { LayerCol.clear(); }
  void ClearShower() { shower2DCol.clear(); }
  void ClearCluster(){ GoodClus3DCol.clear(); BadClus3DCol.clear(); Clus3DCol.clear(); }
  void ClearPFO()    { PFOCol.clear(); }

  


};
#endif
