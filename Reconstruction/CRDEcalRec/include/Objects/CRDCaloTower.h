#ifndef _CRD_CALOTOWER_
#define _CRD_CALOTOWER_

#include <vector>
#include "Objects/CRDCaloBlock.h"
#include "Objects/CRDCaloHitLongiCluster.h"
#include "Objects/Track.h"

namespace CRDEcalEDM{

  class CRDCaloTower{
  public: 
    CRDCaloTower(){};
    void Clear() { blockCol.clear(); }
    inline bool operator == (const CRDCaloTower &x) const{
      return ( module==x.module && stave==x.stave && part==x.part );
    }

    std::vector<CRDEcalEDM::CRDCaloBlock> getBlocks() const { return blockCol; }
    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> getLongiClusterXCol() const { return longiClusXCol; }
    std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> getLongiClusterYCol() const { return longiClusYCol; }
    int getModule() const { return module; }
    int getStave() const { return stave; }
    int getPart() const { return part; }

    void SetTowerID(int _module, int _stave, int _part) { module=_module; stave=_stave; part=_part; }
    void AddBlock( CRDEcalEDM::CRDCaloBlock& _bl ) { blockCol.push_back(_bl); }
    void SetBlocks( std::vector<CRDEcalEDM::CRDCaloBlock>& _bls ) { blockCol=_bls; }
    void SetLongiClusters( std::vector<CRDCaloHitLongiCluster>& _clX, std::vector<CRDCaloHitLongiCluster>& _clY ) { longiClusXCol=_clX; longiClusYCol=_clY; }
    void AddLongiClusterX( CRDCaloHitLongiCluster& _clX ) { longiClusXCol.push_back(_clX); }
    void AddLongiClusterY( CRDCaloHitLongiCluster& _clY ) { longiClusYCol.push_back(_clY); }

  private: 
    int module;
    int stave;
    int part;

    std::vector<CRDEcalEDM::Track> trkCol;
    std::vector<CRDEcalEDM::CRDCaloBlock> blockCol; 
    std::vector<CRDCaloHitLongiCluster> longiClusXCol;
    std::vector<CRDCaloHitLongiCluster> longiClusYCol;

  };

};
#endif
