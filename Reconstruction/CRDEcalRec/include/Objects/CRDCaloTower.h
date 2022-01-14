#ifndef _CRD_CALOTOWER_
#define _CRD_CALOTOWER_

#include <vector>
#include "Objects/CRDCaloBlock.h"

namespace CRDEcalEDM{

  class CRDCaloTower{
  public: 
    CRDCaloTower(){};
    void Clear() { blockCol.clear(); }
    inline bool operator == (const CRDCaloTower &x) const{
      return ( module==x.module && stave==x.stave && part==x.part );
    }

    std::vector<CRDEcalEDM::CRDCaloBlock> getBlocks() const { return blockCol; }
    int getModule() const { return module; }
    int getStave() const { return stave; }
    int getPart() const { return part; }

    void SetTowerID(int _module, int _stave, int _part) { module=_module; stave=_stave; part=_part; }
    void AddBlock( CRDEcalEDM::CRDCaloBlock _bl ) { blockCol.push_back(_bl); }
    void SetBlocks( std::vector<CRDEcalEDM::CRDCaloBlock> _bls ) { blockCol=_bls; }

  private: 
    int module;
    int stave;
    int part;
    std::vector<CRDEcalEDM::CRDCaloBlock> blockCol; 

  };

};
#endif
