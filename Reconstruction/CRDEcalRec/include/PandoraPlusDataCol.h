#ifndef _PANDORAPLUS_DATA_SVC_
#define _PANDORAPLUS_DATA_SVC_

#include "CRDEcalEDMSvc/CRDCaloBar.h"
#include "CRDEcalEDMSvc/CRDCaloBarShower.h"
#include "CRDEcalEDMSvc/CRDCaloBarCluster.h"
#include "CRDEcalEDMSvc/CRDCaloLayer.h"
#include "CRDEcalEDMSvc/CRDCaloHit2DShower.h"
#include "CRDEcalEDMSvc/CRDCaloHit3DShower.h"

#define PI 3.141592653
#define C 299.79  // unit: mm/ns

class PandoraPlusDataCol{
public:
 
  //MCParticle Collections
 
  //Track Collections
  
  //Vertex Collections


  //ECal Collections
  std::vector<CRDEcalEDM::DigiBlock>          blockVec;  //All fired crystal bars, grouped as blocks 
  std::vector<CRDEcalEDM::CRDCaloLayer>       LayerCol;  //Results of EnergySplittingAlg
  std::vector<CRDEcalEDM::CRDCaloHit2DShower> shower2DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> GoodClus3DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> BadClus3DCol;
  std::vector<CRDEcalEDM::CRDCaloHit3DShower> Clus3DCol;


  void PrintLayer();
  void PrintShower();
  void Print3DClus();
  void Clear();
 
  //Hcal Collections


  //PFO Collections


};
#endif
