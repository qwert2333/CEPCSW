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

    int Nbins_alpha;
		int Nbins_rho;

		int th_Layers;  //Only use first th_Layers layer. 
    int th_peak;    //At least th_peak th_peak hits in one track(longiCluster).
    int th_continuetrkN; //Should have at least N continue layers. 
    bool fl_continuetrk; //Should be continue. 
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
