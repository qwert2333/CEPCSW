#ifndef HOUGHCLUSTERINGALG_H
#define HOUGHCLUSTERINGALG_H

#include "CyberDataCol.h"
#include "Tools/Algorithm.h"
#include "TVector2.h"
#include <vector>
using namespace Cyber;

class HoughClusteringAlg: public Cyber::Algorithm{

public: 

  HoughClusteringAlg () {};
  ~HoughClusteringAlg() {};

  class Factory : public Cyber::AlgorithmFactory
  {
  public:
    Cyber::Algorithm* CreateAlgorithm() const{ return new HoughClusteringAlg(); }
  };

  StatusCode ReadSettings(Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode HoughTransformation( std::vector<Cyber::HoughObject>& Hobjects );
  //StatusCode SetLineRange( int module, int slayer, double *range12, double* range34 );
  StatusCode FillHoughSpace( std::vector<Cyber::HoughObject>& Hobjects, 
                             Cyber::HoughSpace& Hspace );
  StatusCode ClusterFinding( std::vector<Cyber::HoughObject>& Hobjects, 
                             Cyber::HoughSpace& Hspace, 
                             std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& longiClusCol );
  StatusCode CleanClusters( std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& m_longiClusCol);
  

private:
	std::vector<Cyber::CaloHalfCluster*> barrel_HalfClusterV;  // Barrel ECAL, V: bars parallel to z axis
  std::vector<Cyber::CaloHalfCluster*> barrel_HalfClusterU;  // Barrel ECAL, U: bars perpendicular to z axis
  std::vector<Cyber::CaloHalfCluster*> endcap0_HalfClusterV; // Endcap ECAL at z~-2900mm, V: bars parallel to x axis
  std::vector<Cyber::CaloHalfCluster*> endcap0_HalfClusterU; // Endcap ECAL at z~-2900mm, U: bars parallel to y axis
  std::vector<Cyber::CaloHalfCluster*> endcap1_HalfClusterV; // Endcap ECAL at z~2900mm, V: bars parallel to x axis
  std::vector<Cyber::CaloHalfCluster*> endcap1_HalfClusterU; // Endcap ECAL at z~2900mm, U: bars parallel to y axis

  std::vector<const Cyber::Calo1DCluster*> m_localMaxVCol;
  std::vector<const Cyber::Calo1DCluster*> m_localMaxUCol;
  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>> m_longiClusVCol;
  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>> m_longiClusUCol;


};
#endif
