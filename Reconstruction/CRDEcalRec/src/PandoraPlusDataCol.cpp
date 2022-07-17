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
  Cluster1DCol.clear();
  Cluster2DCol.clear(); 
  Cluster3DCol.clear();

  return StatusCode::SUCCESS; 
};

StatusCode PandoraPlusDataCol::Clean(){
  for(int i=0; i<bk_TrackCol.size(); i++) 
    if(!bk_TrackCol[i]) { delete bk_TrackCol[i]; bk_TrackCol[i]=NULL; }

  for(int i=0; i<bk_BarCol.size(); i++)
    if(!bk_BarCol[i]) { delete bk_BarCol[i]; bk_BarCol[i]=NULL; }

  for(int i=0; i<bk_BlockCol.size(); i++)
    if(!bk_BlockCol[i]) { bk_BlockCol[i]->Clean(); delete bk_BlockCol[i]; bk_BlockCol[i]=NULL; }

  for(int i=0; i<bk_TowerCol.size(); i++)
    if(!bk_TowerCol[i]) { bk_TowerCol[i]->Clean(); delete bk_TowerCol[i]; bk_TowerCol[i]=NULL; }

  for(int i=0; i<bk_BarClusCol.size(); i++)
    if(!bk_BarClusCol[i]) { bk_BarClusCol[i]->Clean(); delete bk_BarClusCol[i]; bk_BarClusCol[i]=NULL; }

  for(int i=0; i<bk_BarShowerCol.size(); i++)
    if(!bk_BarShowerCol[i]) { bk_BarShowerCol[i]->Clean(); delete bk_BarShowerCol[i]; bk_BarShowerCol[i]=NULL; }

  for(int i=0; i<bk_TransShowerCol.size(); i++)
    if(!bk_TransShowerCol[i]) { bk_TransShowerCol[i]->Clean(); delete bk_TransShowerCol[i]; bk_TransShowerCol[i]=NULL; }

  for(int i=0; i<bk_LongiClusCol.size(); i++)
    if(!bk_LongiClusCol[i]) { bk_LongiClusCol[i]->Clean(); delete bk_LongiClusCol[i]; bk_LongiClusCol[i]=NULL; }

  for(int i=0; i<bk_ClusterCol.size(); i++)
    if(!bk_ClusterCol[i]) { bk_ClusterCol[i]->Clean(); delete bk_ClusterCol[i]; bk_ClusterCol[i]=NULL; }

  for(int i=0; i<bk_Cluster1DCol.size(); i++)
    if(!bk_Cluster1DCol[i]) { bk_Cluster1DCol[i]->Clean(); delete bk_Cluster1DCol[i]; bk_Cluster1DCol[i]=NULL; }

  for(int i=0; i<bk_Cluster2DCol.size(); i++)
    if(!bk_Cluster2DCol[i]) { bk_Cluster2DCol[i]->Clean(); delete bk_Cluster2DCol[i]; bk_Cluster2DCol[i]=NULL; }

  for(int i=0; i<bk_Cluster3DCol.size(); i++)
    if(!bk_Cluster3DCol[i]) { bk_Cluster3DCol[i]->Clean(); delete bk_Cluster3DCol[i]; bk_Cluster3DCol[i]=NULL; }

  Clear();

  return StatusCode::SUCCESS;
}
#endif
