#ifndef CALOTOWER_C
#define CALOTOWER_C

#include "Objects/CaloTower.h"

namespace PandoraPlus{

  void CaloTower::Clear(){
    trkCol.clear();
	  blockCol.clear();
    Cluster2DCol.clear();
    longiClusXCol.clear();
    longiClusYCol.clear();
  }

  void CaloTower::Clean(){
    for(auto iter : trkCol) if(!iter) {delete iter; iter=NULL;}
    trkCol.clear();
    CleanBlock();
    CleanLongiClusters();
    for(auto iter : Cluster2DCol) if(!iter) {delete iter; iter=NULL;}
    Cluster2DCol.clear();
  }

  void CaloTower::CleanBlock() {
    for(auto iter : blockCol) if(!iter) {delete iter; iter=NULL;}
    blockCol.clear();
  }

  void CaloTower::CleanLongiClusters() {
    for(auto iter : longiClusXCol) if(!iter) {delete iter; iter=NULL;}
    for(auto iter : longiClusYCol) if(!iter) {delete iter; iter=NULL;}
    longiClusXCol.clear(); longiClusYCol.clear();
  }

  void CaloTower::Check(){
    for(int i=0; i<trkCol.size(); i++)  if(!trkCol[i]) { trkCol.erase(trkCol.begin()+i); i--; }
    for(int i=0; i<blockCol.size(); i++)  if(!blockCol[i]) { blockCol.erase(blockCol.begin()+i); i--; }
    for(int i=0; i<longiClusXCol.size(); i++)  if(!longiClusXCol[i]) { longiClusXCol.erase(longiClusXCol.begin()+i); i--; }
    for(int i=0; i<longiClusYCol.size(); i++)  if(!longiClusYCol[i]) { longiClusYCol.erase(longiClusYCol.begin()+i); i--; }
  }
  

};
#endif
