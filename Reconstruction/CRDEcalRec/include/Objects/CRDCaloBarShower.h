#ifndef _CRD_CALOBARSHOWER_
#define _CRD_CALOBARSHOWER_

#include <vector>
#include "Objects/CRDCaloBar.h"
#include "Objects/CRDShadowCluster.h"
#include "TVector3.h"
namespace CRDEcalEDM {

  class CRDCaloBarShower{

  public: 
    CRDCaloBarShower() {}; 

    double getlateralmomentofy() const;
    double getlateralmomentofz() const;

    void Clear() { Energy=0; Bars.clear(); TrkCandidateCol.clear(); NeuCandidateCol.clear(); }
    bool isNeighbor(CRDEcalEDM::CRDCaloBar iBar); 
    bool inShower(CRDEcalEDM::CRDCaloBar iBar); 

    dd4hep::Position getPos() const;
    double getE()  const;
    double getT1() const;
    double getT2() const;
    double getWidth() const; 
    int getDlayer() const;
    std::vector<CRDEcalEDM::CRDCaloBar> getBars() const { return Bars; }
    CRDEcalEDM::CRDCaloBar getSeed() const { return Seed; }
    std::vector<CRDEcalEDM::CRDShadowCluster> getTrkCandiCol() const { return TrkCandidateCol; }
    std::vector<CRDEcalEDM::CRDShadowCluster> getNeuCandiCol() const { return NeuCandidateCol; }
    std::vector<CRDEcalEDM::CRDShadowCluster> getAllCandiCol() const; 

    void setBars(std::vector<CRDEcalEDM::CRDCaloBar> _bars){ Bars = _bars; }
    void addBar(CRDEcalEDM::CRDCaloBar _bar) { Bars.push_back(_bar); }
    void setSeed(CRDEcalEDM::CRDCaloBar _seed ) { Seed = _seed; }
    void setPos(dd4hep::Position _pos) { pos = _pos; }
    void addTrkCandidate( CRDEcalEDM::CRDShadowCluster _candi ) { TrkCandidateCol.push_back(_candi); }
    void addNeuCandidate( CRDEcalEDM::CRDShadowCluster _candi ) { NeuCandidateCol.push_back(_candi); }

  private:
		double Energy;
    dd4hep::Position pos;
    std::vector<CRDEcalEDM::CRDCaloBar> Bars;
    CRDEcalEDM::CRDCaloBar Seed;
    std::vector<CRDEcalEDM::CRDShadowCluster> TrkCandidateCol;
    std::vector<CRDEcalEDM::CRDShadowCluster> NeuCandidateCol;
  };


};
#endif
