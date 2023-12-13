#ifndef _TRUTHCLUSTERMERGING_ALG_H
#define _TRUTHCLUSTERMERGING_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;
class TruthClusterMergingAlg: public PandoraPlus::Algorithm{
public: 

  TruthClusterMergingAlg(){};
  ~TruthClusterMergingAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new TruthClusterMergingAlg(); } 

  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();


private: 

  std::vector<const PandoraPlus::Calo3DCluster*> m_EcalClusterCol;
  std::vector<const PandoraPlus::Calo3DCluster*> m_HcalClusterCol;
  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> merged_EcalClusterCol;
  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> merged_HcalClusterCol;
  std::vector<std::shared_ptr<PandoraPlus::PFObject>> merged_CombClusterCol;

};
#endif
