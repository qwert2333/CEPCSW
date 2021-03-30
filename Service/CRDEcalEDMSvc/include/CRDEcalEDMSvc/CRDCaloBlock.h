#ifndef _CRD_CALOBLOCK_
#define _CRD_CALOBLOCK_

#include "CRDEcalEDMSvc/CRDCaloBar.h"
#include "CRDEcalEDMSvc/CRDExpEMShower.h"
namespace CRDEcalEDM {

  class CRDCaloBlock {
  public: 
    CRDCaloBlock() {};

    void Clear(){
      barXCol.clear(); 
      barYCol.clear();
      ForcedEMShower.clear();
    }
    int getNforcedEMShower() const { return ForcedEMShower.size(); }
    std::vector<CRDEcalEDM::CRDCaloBar> getBarXCol() const { return barXCol; }
    std::vector<CRDEcalEDM::CRDCaloBar> getBarYCol() const { return barYCol; }
    std::vector<CRDEcalEDM::CRDExpEMShower> getExpShower() const { return ForcedEMShower; }
    bool isForced() const { return m_isForced; }
    int getMode()  const { return m_ForcedMode; }
    int getDlayer() const { 
      if(barXCol.size()>0) return barXCol[0].getDlayer(); 
      if(barYCol.size()>0) return barYCol[0].getDlayer();
      return -1; 
    }

    void setForced(bool flag) { m_isForced=flag; }
    void setMode(int flag) { m_ForcedMode=flag; }
    void setExpShowerCol( std::vector<CRDEcalEDM::CRDExpEMShower> _showers) { ForcedEMShower = _showers; }
    void setBarCol( std::vector<CRDEcalEDM::CRDCaloBar> _colX, std::vector<CRDEcalEDM::CRDCaloBar> _colY ){ barXCol=_colX, barYCol=_colY; }
 
  private:
    std::vector<CRDEcalEDM::CRDCaloBar> barXCol;
    std::vector<CRDEcalEDM::CRDCaloBar> barYCol;

    bool m_isForced;         //false for no forced shower.
    int m_ForcedMode;       //0 for no forced, 1 for exclusive, 2 for inclusive 
    std::vector<CRDEcalEDM::CRDExpEMShower> ForcedEMShower; 

  };

};
#endif
