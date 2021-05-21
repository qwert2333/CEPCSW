#ifndef _CRD_CALOBLOCK_
#define _CRD_CALOBLOCK_

#include "CRDEcalEDMSvc/CRDCaloBar.h"
#include "CRDEcalEDMSvc/CRDShowerCandidate.h"
namespace CRDEcalEDM {

  class CRDCaloBlock {
  public: 
    CRDCaloBlock() {};

    void Clear(){
      barXCol.clear(); 
      barYCol.clear();
      CandidateCol.clear();
    }
    void ClearCandidate() { CandidateCol.clear(); }
    std::vector<CRDEcalEDM::CRDCaloBar> getBarXCol() const { return barXCol; }
    std::vector<CRDEcalEDM::CRDCaloBar> getBarYCol() const { return barYCol; }
    std::vector<CRDEcalEDM::CRDShowerCandidate> getCandidateCol() const { return CandidateCol; }
    bool isForced() const { return m_isForced; }
    int getDlayer() const { 
      if(barXCol.size()>0) return barXCol[0].getDlayer(); 
      if(barYCol.size()>0) return barYCol[0].getDlayer();
      return -1; 
    }

    void setForced(bool flag) { m_isForced=flag; }
    void setCandidateCol( std::vector<CRDEcalEDM::CRDShowerCandidate> _showers) { CandidateCol = _showers; }
    void setBarCol( std::vector<CRDEcalEDM::CRDCaloBar> _colX, std::vector<CRDEcalEDM::CRDCaloBar> _colY ){ barXCol=_colX, barYCol=_colY; }
 
  private:
    std::vector<CRDEcalEDM::CRDCaloBar> barXCol;
    std::vector<CRDEcalEDM::CRDCaloBar> barYCol;

    bool m_isForced;         //false for no forced shower.
    std::vector<CRDEcalEDM::CRDShowerCandidate> CandidateCol; 

  };

};
#endif
