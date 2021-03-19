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

    double chi2Wi_E;
    double chi2Wi_T;
    double th_chi2;
    double sigmaE;
    double sigmaPos;
    double nMat;
  };

  EnergyTimeMatchingAlg();

  StatusCode Initialize();

  StatusCode RunAlgorithm( EnergyTimeMatchingAlg::Settings& m_settings, PandoraPlusDataCol& m_datasvc );

  CRDEcalEDM::CRDCaloHit2DShower DigiHitsWithPos( CRDEcalEDM::CRDCaloBarShower& barShowerX,  CRDEcalEDM::CRDCaloBarShower& barShowerY  );

  std::vector<CRDEcalEDM::CRDCaloHit2DShower> DigiHitsWithMatching( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol );

  std::vector<CRDEcalEDM::CRDCaloHit2DShower> DigiHitsWithMatchingL2( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol );



  Settings settings;


private: 

};
#endif
