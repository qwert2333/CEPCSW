//=============================================================
// CyberPFA: a PFA developed for CEPC referenece detector
// Ver. CyberPFA-5.0.1(2025.01.09)
//-------------------------------------------------------------
// Data Collection with CyberPFA EDM
//-------------------------------------------------------------
//  Author: Fangyi Guo, Yang Zhang, Weizheng Song, Shengsen Sun
//          (IHEP, CAS)
//  Contact: guofangyi@ihep.ac.cn,
//           sunss@ihep.ac.cn
//=============================================================
#ifndef _PANDORAPLUS_DATA_C
#define _PANDORAPLUS_DATA_C

#include "CyberDataCol.h"

StatusCode CyberDataCol::Clear(){
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
 
  //if(tofCol) delete tofCol;
  //if(dNdxCol) delete dNdxCol;
 
  return StatusCode::SUCCESS; 
};

#endif
