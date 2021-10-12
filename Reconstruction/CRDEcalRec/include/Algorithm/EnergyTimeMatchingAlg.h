#ifndef _ETMATCHING_ALG_H
#define _ETMATCHING_ALG_H

#include "PandoraPlusDataCol.h"
#include "TVector3.h"
using namespace CRDEcalEDM;

class EnergyTimeMatchingAlg{

public: 

  class Settings{
  public: 
    Settings(){};
    void SetInitialValue();
    void SetChi2Weight( double _chi2E, double _chi2T) { chi2Wi_E=_chi2E; chi2Wi_T=_chi2T; }
    void SetUseShadowCluster( bool _flag ) { fl_useShadowClus=_flag; }
    void SetUseChi2( bool _flag ) { fl_UseChi2 = _flag; }
    void SetDebug( int _debug) { Debug=_debug; }

    double chi2Wi_E;
    double chi2Wi_T;
    double th_chi2;
    double sigmaE;
    double sigmaPos = 34.89;  //sqrt(10*10/12 + pow((Tres*C/(2*nMat)),2) )
    double nMat = 2.15;
    double Emip = 0.01;
    bool fl_useShadowClus;
    bool fl_UseChi2; 
    int Debug;
  };

  EnergyTimeMatchingAlg();

  StatusCode Initialize();

  StatusCode RunAlgorithm( EnergyTimeMatchingAlg::Settings& m_settings, PandoraPlusDataCol& m_datasvc );



  StatusCode GetFullMatchedShowers( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol, std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol );

  StatusCode GetMatchedShowersL0( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol, std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol );

  StatusCode GetMatchedShowersL1( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol, std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol );

  StatusCode GetMatchedShowersL2( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol, std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol );


  StatusCode XYShowerMatchingL0( CRDEcalEDM::CRDCaloBarShower& barShowerX, CRDEcalEDM::CRDCaloBarShower& barShowerY, CRDEcalEDM::CRDCaloHitTransShower& outsh ); //1*1 case
  StatusCode XYShowerMatchingL1( CRDEcalEDM::CRDCaloBarShower& barShower1, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerNCol, std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol ); //1*N case without shadow cluster
  StatusCode XYShowerChi2Matching( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol, std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol ); //N*N case without shadow cluster
  StatusCode XYShowerChi2MatchingL1( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol, std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol ); //M*N case without shadow cluster




  Settings settings;


private: 

};
#endif
