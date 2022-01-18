#ifndef _LOCALMAXFINDING_ALG_H
#define _LOCALMAXFINDING_ALG_H

#include "PandoraPlusDataCol.h"

using namespace CRDEcalEDM;

class LocalMaxFindingAlg{
public: 

  class Settings{
  public: 
    Settings() {};

    void SetInitialValue(); 
    double Eth_localMax;
    double Eth_MaxWithNeigh;
  };

  LocalMaxFindingAlg(){};
  ~LocalMaxFindingAlg(){};

  StatusCode Initialize();
  StatusCode RunAlgorithm( LocalMaxFindingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);
  StatusCode ClearAlgorithm();

  StatusCode GetLocalMax( CRDEcalEDM::CRDCaloBlock& m_block);
  StatusCode GetLocalMaxBar( std::vector<CRDEcalEDM::CRDCaloBar>& barCol, std::vector<CRDEcalEDM::CRDCaloBar>& localMaxCol );
  std::vector<CRDEcalEDM::CRDCaloBar>  getNeighbors( CRDEcalEDM::CRDCaloBar& seed, std::vector<CRDEcalEDM::CRDCaloBar>& barCol);

  Settings settings;

private: 

};
#endif
