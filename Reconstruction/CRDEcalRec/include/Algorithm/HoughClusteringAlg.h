#ifndef HOUGHCLUSTERINGALG_H
#define HOUGHCLUSTERINGALG_H

#include "PandoraPlusDataCol.h"
#include "Tools/Algorithm.h"
#include "TVector2.h"
#include <vector>
using namespace PandoraPlus;

class HoughClusteringAlg: public PandoraPlus::Algorithm{

public: 

  HoughClusteringAlg () {};
  ~HoughClusteringAlg() {};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public:
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new HoughClusteringAlg(); }
  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode HoughTransformation(std::vector<PandoraPlus::HoughObject>& Hobjects);
  StatusCode SetLineRange(int module, int slayer, double *range12, double* range34);
  StatusCode FillHoughSpace(std::vector<PandoraPlus::HoughObject>& Hobjects, 
                            PandoraPlus::HoughSpace& Hspace);
  StatusCode ClusterFinding(std::vector<PandoraPlus::HoughObject>& Hobjects, 
                            PandoraPlus::HoughSpace& Hspace, 
                            std::vector<const PandoraPlus::CaloHalfCluster*>& longiClusCol);
  StatusCode CleanClusters(vector<PandoraPlus::CaloHalfCluster*>& m_longiClusCol);
  

private:
	std::vector<PandoraPlus::CaloHalfCluster*> p_HalfClusterV;
  std::vector<PandoraPlus::CaloHalfCluster*> p_HalfClusterU;

  std::vector<const PandoraPlus::Calo1DCluster*> m_localMaxVCol;
  std::vector<const PandoraPlus::Calo1DCluster*> m_localMaxUCol;
  std::vector<const PandoraPlus::CaloHalfCluster*> m_longiClusVCol;
  std::vector<const PandoraPlus::CaloHalfCluster*> m_longiClusUCol;


};
#endif
