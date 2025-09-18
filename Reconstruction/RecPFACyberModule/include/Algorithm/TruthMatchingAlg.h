#ifndef _TRUTHMATCHING_ALG_H
#define _TRUTHMATCHING_ALG_H

#include "Tools/Algorithm.h"

using namespace Cyber;
class TruthMatchingAlg: public Cyber::Algorithm{
public: 

  TruthMatchingAlg(){};
  ~TruthMatchingAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
  public: 
    Cyber::Algorithm* CreateAlgorithm() const{ return new TruthMatchingAlg(); } 

  };

  StatusCode ReadSettings(Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode TruthMatching( std::vector<const Cyber::CaloHalfCluster*>& m_ClUCol,
                            std::vector<const Cyber::CaloHalfCluster*>& m_ClVCol,
                            std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_clusters );

  StatusCode XYClusterMatchingL0( const Cyber::CaloHalfCluster* m_longiClX,
                                  const Cyber::CaloHalfCluster* m_longiClY,
                                  std::shared_ptr<Cyber::Calo3DCluster>& m_clus );

  StatusCode GetMatchedShowersL0( const Cyber::Calo1DCluster* barShowerX,
                                  const Cyber::Calo1DCluster* barShowerY,
                                  Cyber::Calo2DCluster* outsh); //1*1  

  StatusCode GetMatchedShowersL1( const Cyber::Calo1DCluster* shower1,
                                  std::vector<const Cyber::Calo1DCluster*>& showerNCol,
                                  std::vector<Cyber::Calo2DCluster*>& outshCol ); //1*N

private: 

  std::vector<Cyber::Calo3DCluster*> m_towerCol;

  std::vector<const Cyber::CaloHalfCluster*> m_HFClusUCol;
  std::vector<const Cyber::CaloHalfCluster*> m_HFClusVCol;
  std::vector<std::shared_ptr<Cyber::Calo3DCluster>> m_clusterCol;

  CyberDataCol m_bkCol;

};
#endif
