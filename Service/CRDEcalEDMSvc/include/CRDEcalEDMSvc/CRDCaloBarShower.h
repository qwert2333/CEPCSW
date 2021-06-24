#ifndef _CRD_CALOBARSHOWER_
#define _CRD_CALOBARSHOWER_

#include <vector>
#include "CRDEcalEDMSvc/CRDCaloBar.h"
#include "CRDEcalEDMSvc/CRDShowerCandidate.h"
#include "TVector3.h"
namespace CRDEcalEDM {

  class CRDCaloBarShower{

  public: 
    CRDCaloBarShower() {}; 

    void Clear() { Energy=0; Bars.clear(); TrkCandidateCol.clear(); NeuCandidateCol.clear(); }
    bool isNeighbor(CRDEcalEDM::CRDCaloBar iBar); 
    bool inShower(CRDEcalEDM::CRDCaloBar iBar); 

    dd4hep::Position getPos() const;
    double getE()  const;
    double getT1() const;
    double getT2() const;
    std::vector<CRDEcalEDM::CRDCaloBar> getBars() const { return Bars; }
    CRDEcalEDM::CRDCaloBar getSeed() const { return Seed; }
    std::vector<CRDEcalEDM::CRDShowerCandidate> getTrkCandiCol() const { return TrkCandidateCol; }
    std::vector<CRDEcalEDM::CRDShowerCandidate> getNeuCandiCol() const { return NeuCandidateCol; }
    std::vector<CRDEcalEDM::CRDShowerCandidate> getAllCandiCol() const; 


    void setBars(std::vector<CRDEcalEDM::CRDCaloBar> _bars){ Bars = _bars; }
    void addBar(CRDEcalEDM::CRDCaloBar _bar) { Bars.push_back(_bar); }
    void setSeed(CRDEcalEDM::CRDCaloBar _seed ) { Seed = _seed; }
    void setPos(dd4hep::Position _pos) { pos = _pos; }
    void addTrkCandidate( CRDEcalEDM::CRDShowerCandidate _candi ) { TrkCandidateCol.push_back(_candi); }
    void addNeuCandidate( CRDEcalEDM::CRDShowerCandidate _candi ) { NeuCandidateCol.push_back(_candi); }

  private:
		double Energy;
    dd4hep::Position pos;
    std::vector<CRDEcalEDM::CRDCaloBar> Bars;
    CRDEcalEDM::CRDCaloBar Seed;
    std::vector<CRDEcalEDM::CRDShowerCandidate> TrkCandidateCol;
    std::vector<CRDEcalEDM::CRDShowerCandidate> NeuCandidateCol;
  };


};
#endif
