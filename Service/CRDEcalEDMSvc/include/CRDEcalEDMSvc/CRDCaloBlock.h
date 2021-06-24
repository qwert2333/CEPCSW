#ifndef _CRD_CALOBLOCK_
#define _CRD_CALOBLOCK_

#include "CRDEcalEDMSvc/CRDCaloBar.h"
#include "CRDEcalEDMSvc/CRDCaloBarShower.h"
#include "CRDEcalEDMSvc/CRDShowerCandidate.h"
#include "CRDEcalEDMSvc/Track.h"
namespace CRDEcalEDM {

  class CRDCaloBlock {
  public: 
    CRDCaloBlock() {};
    CRDCaloBlock(int _module, int _stave, int _dlayer, int _part): module(_module), stave(_stave), dlayer(_dlayer), part(_part) {};

    void Clear(){
      barXCol.clear(); 
      barYCol.clear();
      barShowerXCol.clear();
      barShowerYCol.clear();
      NeuCandidateCol.clear();
      TrkCandidateCol.clear();
      trkCol.clear();
    }
    void ClearNeuCandidate() { NeuCandidateCol.clear(); }
    void ClearTrkCandidate() { TrkCandidateCol.clear(); }

    inline bool operator == (const CRDCaloBlock &x) const{
      return ( module == x.module  &&
               stave  == x.stave   &&
               part   == x.part    &&
               dlayer == x.dlayer  );
    }


    int getModule() const { return module; }
    int getStave()  const { return stave;  }
    int getDlayer() const { return dlayer; }
    int getPart()   const { return part;   }
    std::vector<CRDEcalEDM::CRDCaloBar> getBarXCol() const { return barXCol; }
    std::vector<CRDEcalEDM::CRDCaloBar> getBarYCol() const { return barYCol; }
    std::vector<CRDEcalEDM::CRDCaloBarShower> getShowerXCol() const { return barShowerXCol; }
    std::vector<CRDEcalEDM::CRDCaloBarShower> getShowerYCol() const { return barShowerYCol; }
    std::vector<CRDEcalEDM::Track> getTrkCol() const {return trkCol; }
    std::vector<CRDEcalEDM::CRDShowerCandidate> getNeuCandidateCol() const { return NeuCandidateCol; }
    std::vector<CRDEcalEDM::CRDShowerCandidate> getTrkCandidateCol() const { return TrkCandidateCol; }
    std::vector<CRDEcalEDM::CRDShowerCandidate> getAllCandidateCol() const; 
    bool isForced() const { return m_isForced; }

    bool MatchTrk( CRDEcalEDM::Track _trk ) const;
  
    void setIDInfo( int _m, int _s, int _d, int _p ){ module=_m; stave=_s; dlayer=_d; part=_p; }
    void setForced(bool flag) { m_isForced=flag; }
    void setNeuCandidateCol( std::vector<CRDEcalEDM::CRDShowerCandidate>& _showers) { NeuCandidateCol = _showers; }
    void setTrkCandidateCol( std::vector<CRDEcalEDM::CRDShowerCandidate>& _showers) { TrkCandidateCol = _showers; }
    void setBarCol( std::vector<CRDEcalEDM::CRDCaloBar>& _colX, std::vector<CRDEcalEDM::CRDCaloBar>& _colY ){ barXCol=_colX, barYCol=_colY; }
    void setTrkCol( std::vector<CRDEcalEDM::Track>& _trkvec ) { trkCol=_trkvec; } 

    void addTrk( CRDEcalEDM::Track& _trk ) { trkCol.push_back(_trk); }

    void PrintCandidates() const; 

  private:
    int module;
    int stave;
    int dlayer;
    int part;

    std::vector<CRDEcalEDM::CRDCaloBar> barXCol;
    std::vector<CRDEcalEDM::CRDCaloBar> barYCol;
    std::vector<CRDEcalEDM::CRDCaloBarShower> barShowerXCol;
    std::vector<CRDEcalEDM::CRDCaloBarShower> barShowerYCol;
    std::vector<CRDEcalEDM::Track> trkCol;

    bool m_isForced;         //false for no forced shower.
    std::vector<CRDEcalEDM::CRDShowerCandidate> NeuCandidateCol; 
    std::vector<CRDEcalEDM::CRDShowerCandidate> TrkCandidateCol;
  };

};
#endif
