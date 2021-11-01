#ifndef _CRD_CALOHIT_2DSHOWER_H
#define _CRD_CALOHIT_2DSHOWER_H
#include <vector>
#include "Objects/CRDCaloBarShower.h"
#include "Objects/CRDShadowCluster.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Vector3f.h"
#include "TVector3.h"

namespace CRDEcalEDM{

  class CRDCaloHit2DShower{
  public: 	
    CRDCaloHit2DShower(CRDEcalEDM::CRDCaloBarShower _bars1, CRDEcalEDM::CRDCaloBarShower _bars2, std::vector<edm4hep::ConstCalorimeterHit> _hits ); 
    CRDCaloHit2DShower(){};

    double getlat() const;

    void Clear(){
      barShowerX.Clear(); barShowerY.Clear(); CaloHits.clear();
    }
    dd4hep::Position getPos() const;
    double getHitsE() const; 
    double getShowerE() const { return barShowerX.getE() + barShowerY.getE();}

    int getModule() const { return module; }
    int getStave()  const { return stave;  }
    int getDlayer() const { return dlayer; }
    int getPart()   const { return part;   }
    double getWidthX() const { return barShowerX.getWidth(); }
    double getWidthY() const { return barShowerY.getWidth(); }
    double getHitsWidth() const; 
    CRDEcalEDM::CRDCaloBarShower getShowerX() const { return barShowerX; }
    CRDEcalEDM::CRDCaloBarShower getShowerY() const { return barShowerY; }
    std::vector<edm4hep::ConstCalorimeterHit> getCaloHits() const { return CaloHits; }
    CRDEcalEDM::CRDShadowCluster getCandidate() const { return Candidate; }
    bool isTrkShower() const { return m_isTrkShower; }
    bool isNeuShower() const { return m_isNeuShower; }

    void setBarShowers( CRDEcalEDM::CRDCaloBarShower _bar1, CRDEcalEDM::CRDCaloBarShower _bar2){ barShowerX=_bar1; barShowerY=_bar2; }
    void setCaloHits( std::vector<edm4hep::ConstCalorimeterHit> _hits) { CaloHits=_hits; }
    void setIDInfo( int _m, int _s, int _d, int _p ){ module=_m; stave=_s; dlayer=_d; part=_p; } 
    void setCandidate( CRDEcalEDM::CRDShadowCluster _candi ) { Candidate = _candi; }
    void setCandidateType( int _flag );  //1 for track-type, 0 for non-track-type.

  private:
    int module;
    int stave;
    int dlayer;
    int part;
    bool m_isTrkShower;
    bool m_isNeuShower;
    CRDEcalEDM::CRDCaloBarShower barShowerX;
    CRDEcalEDM::CRDCaloBarShower barShowerY;
    std::vector<edm4hep::ConstCalorimeterHit> CaloHits;
    CRDEcalEDM::CRDShadowCluster Candidate; 
  };

};
#endif
