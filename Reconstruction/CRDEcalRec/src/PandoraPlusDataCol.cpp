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
  
  return StatusCode::SUCCESS; 
};

StatusCode PandoraPlusDataCol::Clean(){
  // for(vector<PandoraPlus::Track*>::iterator pObj = bk_TrackCol.begin(); pObj != bk_TrackCol.end(); ++pObj)
  // {
  //   delete *pObj; 
  //   *pObj =NULL;
  // }
  // for(vector<PandoraPlus::CaloUnit*>::iterator pObj = bk_BarCol.begin(); pObj != bk_BarCol.end(); ++pObj)
  // {
  //   delete *pObj; 
  //   *pObj =NULL;
  // }

  for(int i=0; i<bk_TrackCol.size(); i++) 
    if(!bk_TrackCol[i]) { delete bk_TrackCol[i]; bk_TrackCol[i]=NULL; }

  for(int i=0; i<bk_BarCol.size(); i++)
    if(!bk_BarCol[i]) { delete bk_BarCol[i]; bk_BarCol[i]=NULL; }

  //for(int i=0; i<bk_LongiClusCol.size(); i++)
  //  if(!bk_LongiClusCol[i]) { bk_LongiClusCol[i]->Clean(); delete bk_LongiClusCol[i]; bk_LongiClusCol[i]=NULL; }

  for(int i=0; i<bk_Cluster1DCol.size(); i++)
    if(!bk_Cluster1DCol[i]) { bk_Cluster1DCol[i]->Clean(); delete bk_Cluster1DCol[i]; bk_Cluster1DCol[i]=NULL; }

  for(int i=0; i<bk_Cluster2DCol.size(); i++)
    if(!bk_Cluster2DCol[i]) { bk_Cluster2DCol[i]->Clean(); delete bk_Cluster2DCol[i]; bk_Cluster2DCol[i]=NULL; }

  for(int i=0; i<bk_ClusterHalfCol.size(); i++)
    if(!bk_ClusterHalfCol[i]) { bk_ClusterHalfCol[i]->Clean(); delete bk_ClusterHalfCol[i]; bk_ClusterHalfCol[i]=NULL; }

  for(int i=0; i<bk_Cluster3DCol.size(); i++)
    if(!bk_Cluster3DCol[i]) { bk_Cluster3DCol[i]->Clean(); delete bk_Cluster3DCol[i]; bk_Cluster3DCol[i]=NULL; }

  Clear();

  return StatusCode::SUCCESS;
}
#endif
