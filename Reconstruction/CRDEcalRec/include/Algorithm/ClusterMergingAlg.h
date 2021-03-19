#ifndef _CLUSTERMERGING_ALG_H
#define _CLUSTERMERGING_ALG_H

#include "PandoraPlusDataCol.h"

using namespace CRDEcalEDM;

class ClusterMergingAlg{
public: 

  class Settings {
  public: 
    Settings() {};
    void SetInitialValue();

  }

  ClusterMergingAlg() {};
  StatusCode Initialize();
  StatusCode RunAlgorithm( ClusterMergingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);

  Settings settings;


};

#endif
