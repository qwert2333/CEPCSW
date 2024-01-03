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
  map_CaloHit.clear();

  map_BarCol.clear();
  map_1DCluster.clear();
  map_HalfCluster.clear();
  map_2DCluster.clear();
  map_CaloCluster.clear();
  map_PFObjects.clear();
  
  return StatusCode::SUCCESS; 
};

#endif
