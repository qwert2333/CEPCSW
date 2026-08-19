#ifndef _TRUTHMATCHING_ALG_H
#define _TRUTHMATCHING_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;
class TruthMatchingAlg: public PandoraPlus::Algorithm{
public: 

  TruthMatchingAlg(){};
  ~TruthMatchingAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new TruthMatchingAlg(); } 

  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode TruthMatching( std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClUCol,
                            std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClVCol,
                            std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>>& m_clusters );

  StatusCode XYClusterMatchingL0( const PandoraPlus::CaloHalfCluster* m_longiClX,
                                  const PandoraPlus::CaloHalfCluster* m_longiClY,
                                  std::shared_ptr<PandoraPlus::Calo3DCluster>& m_clus );

  StatusCode GetMatchedShowersL0( const PandoraPlus::Calo1DCluster* barShowerX,
                                  const PandoraPlus::Calo1DCluster* barShowerY,
                                  PandoraPlus::Calo2DCluster* outsh); //1*1  

  StatusCode GetMatchedShowersL1( const PandoraPlus::Calo1DCluster* shower1,
                                  std::vector<const PandoraPlus::Calo1DCluster*>& showerNCol,
                                  std::vector<PandoraPlus::Calo2DCluster*>& outshCol ); //1*N

private: 

  std::vector<PandoraPlus::Calo3DCluster*> m_towerCol;

  std::vector<const PandoraPlus::CaloHalfCluster*> m_HFClusUCol;
  std::vector<const PandoraPlus::CaloHalfCluster*> m_HFClusVCol;
  std::vector<std::shared_ptr<PandoraPlus::Calo3DCluster>> m_clusterCol;

  PandoraPlusDataCol m_bkCol;

};
#endif
