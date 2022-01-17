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
    void SetUseChi2( bool _flag ) { fl_UseChi2 = _flag; }
    void SetDebug( int _debug) { Debug=_debug; }

    double chi2Wi_E;
    double chi2Wi_T;
    double th_chi2;
    double sigmaE;
    double sigmaPos = 34.89;  //sqrt(10*10/12 + pow((Tres*C/(2*nMat)),2) )
    double nMat = 2.15;
    bool fl_UseChi2; 
    int Debug;
  };

  EnergyTimeMatchingAlg();

  StatusCode Initialize();
  StatusCode RunAlgorithm( EnergyTimeMatchingAlg::Settings& m_settings, PandoraPlusDataCol& m_datasvc );
  StatusCode ClearAlgorithm(); 

  StatusCode GetFullMatchedShowers( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol, std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol );

  StatusCode XYShowerMatchingL0( CRDEcalEDM::CRDCaloHitLongiCluster& m_longiClX, CRDEcalEDM::CRDCaloHitLongiCluster& m_longiClY, CRDEcalEDM::CRDCaloHit3DCluster& m_clus);
  StatusCode XYShowerMatchingL1( CRDEcalEDM::CRDCaloHitLongiCluster& m_longiCl1, std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClN, std::vector<CRDEcalEDM::CRDCaloHit3DCluster>& m_clusters );
  StatusCode XYShowerChi2Matching( std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClXCol, std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClYCol, double** chi2 );
  StatusCode XYShowerChi2MatchingL1( std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClXCol, std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClYCol, double** chi2 );


  Settings settings;


private: 

};
#endif
