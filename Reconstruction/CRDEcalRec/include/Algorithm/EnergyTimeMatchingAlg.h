#ifndef ETMATCHING_ALG_H
#define ETMATCHING_ALG_H

#include "PandoraPlusDataCol.h"
#include "Tools/Algorithm.h"

#include "TVector3.h"
using namespace PandoraPlus;

class EnergyTimeMatchingAlg: public PandoraPlus::Algorithm{

public: 

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public:
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new EnergyTimeMatchingAlg(); }

  };

  EnergyTimeMatchingAlg(){};
  ~EnergyTimeMatchingAlg(){};

  StatusCode ReadSettings( PandoraPlus::Settings& m_settings);
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datasvc );
  StatusCode ClearAlgorithm(); 


  StatusCode XYClusterMatchingL0( const PandoraPlus::CaloHalfCluster* m_longiClX, const PandoraPlus::CaloHalfCluster* m_longiClY, PandoraPlus::Calo3DCluster* m_clus);
  StatusCode XYClusterMatchingL1( const PandoraPlus::CaloHalfCluster* m_longiCl1, std::vector<const PandoraPlus::CaloHalfCluster*>& m_longiClN, std::vector<PandoraPlus::Calo3DCluster*>& m_clusters );
  StatusCode XYClusterMatchingL2( std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClUCol, std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClVCol, std::vector<PandoraPlus::Calo3DCluster*>& m_clusters );
  StatusCode XYClusterMatchingL3( std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClUCol, std::vector<const PandoraPlus::CaloHalfCluster*>& m_ClVCol, std::vector<PandoraPlus::Calo3DCluster*>& m_clusters );


  StatusCode GetFullMatchedShowers( std::vector<const PandoraPlus::Calo1DCluster*>& barShowerXCol, std::vector<const PandoraPlus::Calo1DCluster*>& barShowerYCol, std::vector<PandoraPlus::Calo2DCluster*>& outshCol );
  StatusCode GetMatchedShowersL0(const PandoraPlus::Calo1DCluster* barShowerX, const PandoraPlus::Calo1DCluster* barShowerY, PandoraPlus::Calo2DCluster* outsh); //1*1
  StatusCode GetMatchedShowersL1(const PandoraPlus::Calo1DCluster* shower1, std::vector<const PandoraPlus::Calo1DCluster*>& showerNCol, std::vector<PandoraPlus::Calo2DCluster*>& outshCol ); //1*N
  StatusCode GetMatchedShowersL2(std::vector<const PandoraPlus::Calo1DCluster*>& barShowerXCol, std::vector<const PandoraPlus::Calo1DCluster*>& barShowerYCol, std::vector<PandoraPlus::Calo2DCluster*>& outshCol ); 

  StatusCode GetMatchedShowersL3(std::vector<const PandoraPlus::Calo1DCluster*>& barShowerXCol, std::vector<const PandoraPlus::Calo1DCluster*>& barShowerYCol, std::vector<PandoraPlus::Calo2DCluster*>& outshCol );

  double** GetClusterChi2Map(std::vector<std::vector<const PandoraPlus::Calo1DCluster*>>& barShowerXCol, std::vector<std::vector<const PandoraPlus::Calo1DCluster*>>& barShowerYCol);

private: 
  std::vector<PandoraPlus::Calo3DCluster*> m_towerCol; 

  std::vector<const PandoraPlus::CaloHalfCluster*> m_HFClusUCol;
  std::vector<const PandoraPlus::CaloHalfCluster*> m_HFClusVCol;
  std::vector<PandoraPlus::Calo2DCluster*> m_transhowerCol; 
  std::vector<PandoraPlus::Calo3DCluster*> m_clusterCol; 


};
#endif
