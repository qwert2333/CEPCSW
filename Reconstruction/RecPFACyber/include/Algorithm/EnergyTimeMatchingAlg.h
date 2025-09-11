#ifndef ETMATCHING_ALG_H
#define ETMATCHING_ALG_H

#include "CyberDataCol.h"
#include "Tools/Algorithm.h"

#include "TVector3.h"
using namespace Cyber;

class EnergyTimeMatchingAlg: public Cyber::Algorithm{

public: 

  class Factory : public Cyber::AlgorithmFactory
  {
  public:
    Cyber::Algorithm* CreateAlgorithm() const{ return new EnergyTimeMatchingAlg(); }

  };

  EnergyTimeMatchingAlg(){};
  ~EnergyTimeMatchingAlg(){};

  StatusCode ReadSettings( Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datasvc );
  StatusCode ClearAlgorithm(); 

  StatusCode Matching(std::vector<const Cyber::CaloHalfCluster*> m_HFClusUCol, 
                      std::vector<const Cyber::CaloHalfCluster*> m_HFClusVCol,
                      std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_clusterCol );

  StatusCode PatternMatrixCalculation(std::vector<const Cyber::CaloHalfCluster*> m_HFClusUCol,
                                      std::vector<const Cyber::CaloHalfCluster*> m_HFClusVCol,
                                      vector<vector<int>>& matrix );

  StatusCode Chi2MatrixCalculation( std::vector<const Cyber::CaloHalfCluster*> m_HFClusUCol,
                                    std::vector<const Cyber::CaloHalfCluster*> m_HFClusVCol,
                                    vector<vector<double>>& matrix, 
                                    std::vector< std::pair<int, int> >& chi2order); 

  StatusCode ParameterMatrixCalculation(const int& M, const int& N, vector<vector<int>>& parMatrix);

  StatusCode SimplityMatrix(vector<vector<int>>& parMatrix, vector<double>& Eij, vector<vector<int>>& pattern);

  double SolveMatrix(const int& M, const int& N, vector<vector<int>>& parMatrix, vector<double>& Eij, std::vector<double>& En_clusters); 

  StatusCode leastSquares(const vector<vector<int>>& A, const vector<double>& b, vector<double>& x);

  vector<vector<double>> GetClusterChi2Map( std::vector<std::vector<const Cyber::Calo1DCluster*>>& barShowerUCol,
                                            std::vector<std::vector<const Cyber::Calo1DCluster*>>& barShowerVCol );

  StatusCode ClusterBuilding( std::vector<std::vector<double>> Ematrix, 
                              std::vector<const Cyber::CaloHalfCluster*>& m_HFClusUCol, 
                              std::vector<const Cyber::CaloHalfCluster*>& m_HFClusVCol, 
                              std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_clusterCol);

  StatusCode XYClusterMatchingL0( const Cyber::CaloHalfCluster* m_longiClU,
                                  const Cyber::CaloHalfCluster* m_longiClV,
                                  std::shared_ptr<Cyber::Calo3DCluster>& m_clus );

  StatusCode GetMatchedShowersL0( const Cyber::Calo1DCluster* barShowerU,
                                  const Cyber::Calo1DCluster* barShowerV,
                                  Cyber::Calo2DCluster* outsh );

  StatusCode GetMatchedShowersL1( const Cyber::Calo1DCluster* shower1,
                                  std::vector<const Cyber::Calo1DCluster*>& showerNCol,
                                  std::vector<Cyber::Calo2DCluster*>& outshCol );

  StatusCode ClusterReconnecting(std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_clusterCol);


private: 
  std::vector<Cyber::Calo3DCluster*> m_towerCol; 

  std::vector<const Cyber::CaloHalfCluster*> m_HFClusUCol;
  std::vector<const Cyber::CaloHalfCluster*> m_HFClusVCol;
  std::vector<const Cyber::CaloHalfCluster*> m_emptyHFClusUCol;
  std::vector<const Cyber::CaloHalfCluster*> m_emptyHFClusVCol;
  std::vector<std::shared_ptr<Cyber::Calo3DCluster>> m_clusterCol; 

  CyberDataCol m_bkCol;

};
#endif
