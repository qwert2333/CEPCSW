#ifndef HCALCLUSTERING_ALG_H
#define HCALCLUSTERING_ALG_H

#include "Tools/Algorithm.h"

using namespace Cyber;

class HcalClusteringAlg : public Cyber::Algorithm{
public: 

  HcalClusteringAlg(){};
  ~HcalClusteringAlg(){};
 
  class Factory : public Cyber::AlgorithmFactory
  {
  public:
    Cyber::Algorithm* CreateAlgorithm() const{ return new HcalClusteringAlg(); }
 
  };
 
  StatusCode ReadSettings(Cyber::Settings& m_settings); 
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();
  
  template<typename T1, typename T2> StatusCode Clustering(std::vector<T1*>& m_input, std::vector<std::shared_ptr<T2>>& m_output);
  StatusCode LongiConeLinking( const std::map<int, std::vector<Cyber::CaloHit*> >& orderedHit, std::vector<std::shared_ptr<Cyber::Calo3DCluster> >& ClusterCol );
  StatusCode MergeToCluster( std::vector<const Cyber::CaloHit*>& m_isohits, std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_clusters );
private:

//   std::vector<std::shared_ptr<Cyber::CaloHit>> m_hcalHits;
};
#endif
