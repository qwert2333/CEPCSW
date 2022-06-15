#ifndef HOUGHCLUSTERINGALG_H
#define HOUGHCLUSTERINGALG_H

#include "PandoraPlusDataCol.h"
#include "Tools/Algorithm.h"
#include "TVector2.h"
#include "TF1.h"
#include "TH2.h"
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
  StatusCode Initialize();
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode ConformalTransformation(std::vector<PandoraPlus::HoughObject>& m_Hobjects); 
  StatusCode HoughTransformation(std::vector<PandoraPlus::HoughObject>& m_Hobjects);
	StatusCode FillHoughSpace(std::vector<PandoraPlus::HoughObject>& m_Hobjects, PandoraPlus::HoughSpace& m_Hspace);
  StatusCode FindingHills(PandoraPlus::HoughSpace& m_Hspace);
  StatusCode Transform2Clusters(PandoraPlus::HoughSpace& m_Hspace, std::vector<PandoraPlus::HoughObject>& m_Hobjects, std::vector<const PandoraPlus::LongiCluster*>& m_longiClusCol );
  StatusCode CleanClusters( std::vector<PandoraPlus::LongiCluster*>& m_longiClusCol );
  int ExpandingPeak(TH2 &houghMap, int index_a, int index_b, PandoraPlus::HoughSpace::HoughHill& hill); 

private:
	

};
#endif
