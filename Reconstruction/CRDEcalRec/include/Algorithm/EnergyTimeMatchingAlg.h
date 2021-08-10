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
    void SetUseCandidate( int _level ) { lv_useCandidate=_level; }
    void SetUseChi2( bool _flag ) { UseChi2 = _flag; }

    double chi2Wi_E;
    double chi2Wi_T;
    double th_chi2;
    double sigmaE;
    double sigmaPos = 34.89;  //sqrt(10*10/12 + pow((Tres*C/(2*nMat)),2) )
    double nMat = 2.15;
    double Emip = 0.01;
    int lv_useCandidate; //0: Don't use. 1: Use track and neutral candidates together. 2: Use track and neutral candidates separately. 
    bool UseChi2; 
    int Debug;
  };

  EnergyTimeMatchingAlg();

  StatusCode Initialize();

  StatusCode RunAlgorithm( EnergyTimeMatchingAlg::Settings& m_settings, PandoraPlusDataCol& m_datasvc );



  StatusCode GetFullMatchedShowers( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol, std::vector<CRDEcalEDM::CRDCaloHit2DShower>& outshCol );

  StatusCode GetMatchedShowersL0( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol, std::vector<CRDEcalEDM::CRDCaloHit2DShower>& outshCol );

  StatusCode GetMatchedShowersL1( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol, std::vector<CRDEcalEDM::CRDCaloHit2DShower>& outshCol );

  StatusCode GetMatchedShowersL2( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol, std::vector<CRDEcalEDM::CRDCaloHit2DShower>& outshCol );


  StatusCode XYShowerMatchingL0( CRDEcalEDM::CRDCaloBarShower& barShowerX, CRDEcalEDM::CRDCaloBarShower& barShowerY, CRDEcalEDM::CRDCaloHit2DShower& outsh ); //1*1 case
  StatusCode XYShowerMatchingL1( CRDEcalEDM::CRDCaloBarShower& barShower1, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerNCol, std::vector<CRDEcalEDM::CRDCaloHit2DShower>& outshCol ); //1*N case without candidates
  StatusCode XYShowerChi2Matching( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol, std::vector<CRDEcalEDM::CRDCaloHit2DShower>& outshCol ); //N*N case without candidates
  StatusCode XYShowerChi2MatchingL1( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol, std::vector<CRDEcalEDM::CRDCaloHit2DShower>& outshCol ); //M*N case without candidates




  Settings settings;


private: 

};
#endif
