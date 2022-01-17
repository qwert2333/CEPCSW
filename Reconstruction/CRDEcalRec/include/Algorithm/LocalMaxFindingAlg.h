#ifndef _FINDINGLOCALMAX_ALG_H
#define _FINDINGLOCALMAX_ALG_H

#include "PandoraPlusDataCol.h"

using namespace CRDEcalEDM;

class LocalMaxFindingAlg{
public: 

  class Settings{
  public: 
    void SetInitialValue(); 

  }

  LocalMaxFindingAlg(){};
  ~LocalMaxFindingAlg(){};

  StatusCode Initialize();
  StatusCode RunAlgorithm( LocalMaxFindingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);
  StatusCode ClearAlgorithm();

  StatusCode GetLocalMax( CRDEcalEDM::CRDCaloBlock& m_block);
  StatusCode GetLocalMaxBar( std::vector<CRDEcalEDM::CRDCaloBar>& barCol, std::vector<CRDEcalEDM::CRDCaloBar>& localMaxCol );


  Settings settings;

private: 

};
#endif
