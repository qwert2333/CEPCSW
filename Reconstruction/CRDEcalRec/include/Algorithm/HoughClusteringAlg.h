#ifndef HOUGHCLUSTERINGALG_H
#define HOUGHCLUSTERINGALG_H

#include "PandoraPlusDataCol.h"
#include "TVector2.h"
#include "TF1.h"
#include "TH2.h"

using namespace CRDEcalEDM;
class HoughClusteringAlg{

public: 

  class Settings{
  public:
    Settings(){};
    void SetInitialValue();

    //Only use first th_Layers layer.
		int th_Layers;  

    //Houghspace setting
    int Nbins_alpha;
		int Nbins_rho;
    int HoughBinDivide; 
    int th_peak;    //At least th_peak th_peak hits in one track(longiCluster).

    //Good Cluster settings
    int th_continuetrkN; //Should have at least N continue layers. 
    float th_AxisE;
    float th_intercept;

    float th_dAlpha1;
    float th_dAlpha2;
    float th_dRho1;
    float th_dRho2;
    float th_overlapE1; 
    float th_overlapE2; 

	};

	HoughClusteringAlg () {};
	~HoughClusteringAlg() {};

  StatusCode Initialize();
  StatusCode RunAlgorithm( HoughClusteringAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);
  StatusCode ClearAlgorithm();

  StatusCode ConformalTransformation(std::vector<CRDEcalEDM::CRDHoughObject>& m_Hobjects); 
  StatusCode HoughTransformation(std::vector<CRDEcalEDM::CRDHoughObject>& m_Hobjects);
	StatusCode FillHoughSpace(std::vector<CRDEcalEDM::CRDHoughObject>& m_Hobjects, CRDEcalEDM::CRDHoughSpace& m_Hspace);
  StatusCode FindingHills(CRDEcalEDM::CRDHoughSpace& m_Hspace);
  StatusCode Transform2Clusters(CRDEcalEDM::CRDHoughSpace& m_Hspace, std::vector<CRDEcalEDM::CRDHoughObject>& m_Hobjects, std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClusCol );
  StatusCode CleanClusters( std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClusCol );
  int ExpandingPeak(TH2 &houghMap, int index_a, int index_b, CRDEcalEDM::CRDHoughSpace::HoughHill& hill); 

  Settings settings;

private:
	

};
#endif
