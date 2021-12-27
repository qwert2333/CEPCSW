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
		int Layers; 
	};

	HoughClusteringAlg () {};
	~HoughClusteringAlg() {};

  StatusCode Initialize();
  StatusCode RunAlgorithm( BasicClusterIDAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);
  StatusCode ClearAlgorithm();

  StatusCode ConformalTransformation(std::vector<CRDEcalEDM::CRDCaloBarShower>& m_localMax, std::vector<CRDEcalEDM::CRDHoughObject>& m_Hobjects); 
  StatusCode HoughTransformation(std::vector<CRDEcalEDM::CRDHoughObject>& m_Hobjects);
	StatusCode FillHoughSpace(std::vector<CRDEcalEDM::CRDHoughObject>& m_Hobjects, CRDEcalEDM::CRDHoughSpace& m_Hspace);
  StatusCode MergingHills(CRDEcalEDM::CRDHoughSpace& m_Hspace);
  


  Settings settings;

private:
	

};
#endif
