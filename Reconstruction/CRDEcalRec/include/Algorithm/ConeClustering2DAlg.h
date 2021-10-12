#ifndef _CONECLUSTERING2D_ALG_H
#define _CONECLUSTERING2D_ALG_H

#include "PandoraPlusDataCol.h"

//This could be replaced by other better track finding algorithm
class ConeClustering2DAlg {

public: 
  class Settings{
  public: 
    Settings(){};
    void SetInitialValue();

    int th_stopLayer; 
    double th_ConeTheta; 
    double th_ConeR; 
  };

  ConeClustering2DAlg(){};
  ~ConeClustering2DAlg(){};

  StatusCode Initialize();
  StatusCode RunAlgorithm( ConeClustering2DAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);
  StatusCode ClearAlgorithm(); 

  StatusCode LongiConeLinking( std::map<int, std::vector<CRDEcalEDM::CRDCaloBarShower> >& orderedShower, std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& ClusterCol);

  Settings settings;

private: 


};
#endif
