#ifndef CALOTOWER_H
#define CALOTOWER_H

#include <vector>
#include "Objects/CaloBlock.h"
#include "Objects/LongiCluster.h"
#include "Objects/Track.h"

namespace PandoraPlus{

  class CaloTower{
  public: 
    CaloTower(){};
    ~CaloTower() { Clean(); }

    void Clear();
    void Check();
    void CleanBlock();
    void CleanLongiClusters(); 
    void Clean();

    //inline bool operator == (const CaloTower &x) const{
    //  return ( module==x.module && stave==x.stave && part==x.part );
    //}

    std::vector<const CaloBlock*> getBlocks() const { return blockCol; }
    std::vector<const LongiCluster*> getLongiClusterXCol() const { return longiClusXCol; }
    std::vector<const LongiCluster*> getLongiClusterYCol() const { return longiClusYCol; }

    int getModule() const { return module; }
    int getStave() const { return stave; }
    int getPart() const { return part; }

    void setTowerID(int _module, int _stave, int _part) { module=_module; stave=_stave; part=_part; }
    void addBlock( const CaloBlock* _bl ) { blockCol.push_back(_bl); }
    void setBlocks( std::vector<const CaloBlock*> _bls ) { blockCol=_bls; }
    void setLongiClusters( std::vector<const LongiCluster*>& _clX, std::vector<const LongiCluster*>& _clY ) { longiClusXCol=_clX; longiClusYCol=_clY; }
    void addLongiClusterX( const LongiCluster* _clX ) { longiClusXCol.push_back(_clX); }
    void addLongiClusterY( const LongiCluster* _clY ) { longiClusYCol.push_back(_clY); }


  private: 
    int module;
    int stave;
    int part;

    std::vector<const Track *> trkCol;
    std::vector<const CaloBlock*> blockCol; 
    std::vector<const LongiCluster*> longiClusXCol;
    std::vector<const LongiCluster*> longiClusYCol;

  };

};
#endif
