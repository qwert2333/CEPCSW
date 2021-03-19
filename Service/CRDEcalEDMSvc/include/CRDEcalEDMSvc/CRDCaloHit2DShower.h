#ifndef _CRD_CALOHIT_2DSHOWER_H
#define _CRD_CALOHIT_2DSHOWER_H
#include <vector>
#include "CRDEcalEDMSvc/CRDCaloBarShower.h"
#include "edm4hep/CalorimeterHit.h"
#include "TVector3.h"

namespace CRDEcalEDM{

  class CRDCaloHit2DShower{
  public: 	
    CRDCaloHit2DShower(CRDEcalEDM::CRDCaloBarShower _bars1, CRDEcalEDM::CRDCaloBarShower _bars2, std::vector<edm4hep::ConstCalorimeterHit> _hits ); 
    CRDCaloHit2DShower(){};


    void Clear(){
      barShowerX.Clear(); barShowerY.Clear(); CaloHits.clear();
    }
    dd4hep::Position getPos() const;
    double getHitsE() const; 
    double getShowerE() const { return barShowerX.getE() + barShowerY.getE();}

    int getModule() const { return module; }
    int getDlayer() const { return dlayer; }
    int getPart()   const { return part;   }
    int getBlock()  const { return block;  }
    CRDEcalEDM::CRDCaloBarShower getShowerX() const { return barShowerX; }
    CRDEcalEDM::CRDCaloBarShower getShowerY() const { return barShowerY; }
    std::vector<edm4hep::ConstCalorimeterHit> getCaloHits() const { return CaloHits; }

    void setBarShowers( CRDEcalEDM::CRDCaloBarShower _bar1, CRDEcalEDM::CRDCaloBarShower _bar2){ barShowerX=_bar1; barShowerY=_bar2; }
    void setCaloHits( std::vector<edm4hep::ConstCalorimeterHit> _hits) { CaloHits=_hits; }
    void setIDInfo( int _m, int _d, int _p, int _b ){ module=_m; dlayer=_d; part=_p; block=_b; } 

  private:
    int module;
    int dlayer;
    int part;
    int block;
    CRDEcalEDM::CRDCaloBarShower barShowerX;
    CRDEcalEDM::CRDCaloBarShower barShowerY;
    std::vector<edm4hep::ConstCalorimeterHit> CaloHits;
  };

};
#endif
