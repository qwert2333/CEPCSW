#ifndef _TRUTHCLUSTERMERGING_ALG_H
#define _TRUTHCLUSTERMERGING_ALG_H

#include "Tools/Algorithm.h"

using namespace Cyber;
class TruthClusterMergingAlg: public Cyber::Algorithm{
public: 

  TruthClusterMergingAlg(){};
  ~TruthClusterMergingAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
  public: 
    Cyber::Algorithm* CreateAlgorithm() const{ return new TruthClusterMergingAlg(); } 

  };

  StatusCode ReadSettings(Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();


private: 

  std::vector<const Cyber::Calo3DCluster*> m_EcalClusterCol;
  std::vector<const Cyber::Calo3DCluster*> m_HcalClusterCol;
  std::vector<std::shared_ptr<Cyber::Calo3DCluster>> merged_EcalClusterCol;
  std::vector<std::shared_ptr<Cyber::Calo3DCluster>> merged_HcalClusterCol;
  std::vector<std::shared_ptr<Cyber::PFObject>> merged_CombClusterCol;

};
#endif
