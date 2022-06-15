#ifndef _PANDORAPLUS_DATA_C
#define _PANDORAPLUS_DATA_C

#include "PandoraPlusDataCol.h"

StatusCode PandoraPlusDataCol::Clear(){
  collectionMap_MC.clear();
  collectionMap_CaloHit.clear();
  collectionMap_Vertex.clear();
  collectionMap_Track.clear();
  collectionMap_CaloRel.clear();
  collectionMap_TrkRel.clear();

  TrackCol.clear();
  BarCol.clear();
  BlockCol.clear();
  TowerCol.clear(); 

  return StatusCode::SUCCESS; 
};

#endif
