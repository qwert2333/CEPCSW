#ifndef HCALCLUSTERING_ALG_H
#define HCALCLUSTERING_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;

class HcalClusteringAlg : public PandoraPlus::Algorithm{
public: 

  HcalClusteringAlg(){};
  ~HcalClusteringAlg(){};
 
  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public:
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new HcalClusteringAlg(); }
 
  };
 
  StatusCode ReadSettings(PandoraPlus::Settings& m_settings); 
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();
  
  template<typename T1, typename T2> StatusCode Clustering(std::vector<T1*>& m_input, std::vector<std::shared_ptr<T2>>& m_output);
  StatusCode LongiConeLinking( const std::map<int, std::vector<PandoraPlus::CaloHit*> >& orderedHit, std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster> >& ClusterCol );

private:

//   std::vector<std::shared_ptr<PandoraPlus::CaloHit>> m_hcalHits;
};
#endif