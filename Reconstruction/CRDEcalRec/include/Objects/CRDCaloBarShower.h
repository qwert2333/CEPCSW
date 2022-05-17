#ifndef _CRD_CALOBARSHOWER_
#define _CRD_CALOBARSHOWER_

#include <vector>
#include "Objects/CRDCaloBar.h"
#include "Objects/CRDShadowCluster.h"
namespace CRDEcalEDM {

  class CRDCaloBarShower{

  public: 
    CRDCaloBarShower() {}; 

    inline bool operator == (const CRDCaloBarShower &x) const{
      return Bars == x.Bars;
    }
    void Clear() { Energy=0; Bars.clear(); TrkShadowClusCol.clear(); NeuShadowClusCol.clear(); }
    bool isNeighbor(CRDEcalEDM::CRDCaloBar iBar); 
    bool inShower(CRDEcalEDM::CRDCaloBar iBar); 

    TVector3 getPos() const;
    double getE()  const;
    double getT1() const;
    double getT2() const;
    double getWidth() const; 
    int getModule() const { return module; }
    int getStave()  const { return stave;  }
    int getDlayer() const { return dlayer; }
    int getPart()   const { return part;   }
    int getSlayer() const { return slayer; }
    std::vector<CRDEcalEDM::CRDCaloBar> getBars() const { return Bars; }
    CRDEcalEDM::CRDCaloBar getSeed() const { return Seed; }
    std::vector<CRDEcalEDM::CRDShadowCluster> getTrkCandiCol() const { return TrkShadowClusCol; }
    std::vector<CRDEcalEDM::CRDShadowCluster> getNeuCandiCol() const { return NeuShadowClusCol; }
    std::vector<CRDEcalEDM::CRDShadowCluster> getAllCandiCol() const; 

    void setIDInfo( int _m, int _s, int _d, int _p, int _sl ){ module=_m; stave=_s; dlayer=_d; part=_p; slayer=_sl; }
    void setIDInfo(); 
    void setBars(std::vector<CRDEcalEDM::CRDCaloBar>& _bars){ Bars = _bars; }
    void addBar(CRDEcalEDM::CRDCaloBar& _bar) { Bars.push_back(_bar); }
    void setSeed(CRDEcalEDM::CRDCaloBar& _seed ) { Seed = _seed; }
    void setPos(TVector3& _pos) { pos = _pos; }
    void addTrkShadowClus( CRDEcalEDM::CRDShadowCluster& _candi ) { TrkShadowClusCol.push_back(_candi); }
    void addNeuShadowClus( CRDEcalEDM::CRDShadowCluster& _candi ) { NeuShadowClusCol.push_back(_candi); }

  private:
    int module;
    int stave;
    int part;
    int dlayer;
    int slayer;
		double Energy;
    TVector3 pos;
    std::vector<CRDEcalEDM::CRDCaloBar> Bars;
    CRDEcalEDM::CRDCaloBar Seed;
    std::vector<CRDEcalEDM::CRDShadowCluster> TrkShadowClusCol;
    std::vector<CRDEcalEDM::CRDShadowCluster> NeuShadowClusCol;
  };


};
#endif
