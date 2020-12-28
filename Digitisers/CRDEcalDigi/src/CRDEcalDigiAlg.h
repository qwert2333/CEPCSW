#ifndef CRD_ECAL_DIGI_ALG_H
#define CRD_ECAL_DIGI_ALG_H
#include "CRDEcalDigiEDM.h"
#include "CRDEcalPreRecAlg.h"

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/SimCalorimeterHitConst.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/Segmentations.h> 
#include "DetInterface/IGeomSvc.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TString.h"
#include "TH3.h"
#include "TH1.h"

#define C 299.79  // unit: mm/ns
#define PI 3.141592653

class CRDEcalDigiAlg : public GaudiAlgorithm
{
 
public:
 
  CRDEcalDigiAlg(const std::string& name, ISvcLocator* svcLoc);
 
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual StatusCode initialize() ;
 
  /** Called for every event - the working horse.
   */
  virtual StatusCode execute() ; 
 
  /** Called after data processing for clean up.
   */
  virtual StatusCode finalize() ;
 

	std::vector<edm4hep::SimCalorimeterHit> MergeHits(const edm4hep::SimCalorimeterHitCollection* m_col);
	double GetBarLength(CRDEcalDigiEDM::DigiBar bar);
	dd4hep::Position GetCellPos(dd4hep::Position pos, CRDEcalDigiEDM::DigiBar bar);
	edm4hep::SimCalorimeterHit find(edm4hep::SimCalorimeterHitCollection* m_col, dd4hep::Position pos);
	edm4hep::SimCalorimeterHit find(std::vector<edm4hep::SimCalorimeterHit> m_col, unsigned long long cellid);
	unsigned long int coder(CRDEcalDigiEDM::DigiBar bar);

	std::vector<edm4hep::ConstCalorimeterHit> DigiHitsWithPos(std::vector<CRDEcalDigiEDM::DigiBar>& m_block);
	std::vector<edm4hep::ConstCalorimeterHit> DigiHitsWithTime(std::vector<CRDEcalDigiEDM::DigiBar>& m_block);
	std::vector<edm4hep::ConstCalorimeterHit> DigiHitsWithEnergy(std::vector<CRDEcalDigiEDM::DigiBar>& m_block, std::vector<CRDEcalDigiEDM::BarCollection>& barShowerX, std::vector<CRDEcalDigiEDM::BarCollection>& barShowerY);
	std::vector<edm4hep::ConstCalorimeterHit> DigiHitsWithMatching(std::vector<CRDEcalDigiEDM::BarCollection>& barShowerXCol, std::vector<CRDEcalDigiEDM::BarCollection>& barShowerYCol);
	std::vector<edm4hep::ConstCalorimeterHit> DigiHitsWithMatchingL2(std::vector<CRDEcalDigiEDM::BarCollection>& barShowerXCol, std::vector<CRDEcalDigiEDM::BarCollection>& barShowerYCol);


	void Clear();
	void ClearPreRec();

protected:

  SmartIF<IGeomSvc> m_geosvc;
  typedef std::vector<float> FloatVec;

	int _nEvt ;
	float m_length;
	TRandom3 rndm;
	TFile* m_wfile;
	TTree* t_SimBar;
	TTree* t_SimTruth;
	TTree* t_Rec;
	TTree *t_SimCont;
	TTree *t_PreRec;
	
	FloatVec m_step_x, m_step_y, m_step_z, m_step_E, m_step_T1, m_step_T2, m_stepBar_x, m_stepBar_y, m_stepBar_z;
	FloatVec m_simBar_x, m_simBar_y, m_simBar_z, m_simBar_T1, m_simBar_T2, m_simBar_Q1, m_simBar_Q2, m_simBar_dlayer, m_simBar_part, m_simBar_block, m_simBar_slayer;
	FloatVec m_simTruth_x, m_simTruth_y, m_simTruth_z, m_simTruth_E, m_simTruth_dlayer, m_simTruth_slayer;
	FloatVec m_Rec_x, m_Rec_y, m_Rec_z, m_Rec_E;

	FloatVec m_PreRec_Bar0x, m_PreRec_Bar0y, m_PreRec_Bar0z, m_PreRec_Bar0E, m_PreRec_Bar1x, m_PreRec_Bar1y, m_PreRec_Bar1z, m_PreRec_Bar1E;
	FloatVec m_PreRec_EshowerX, m_PreRec_XshowerX, m_PreRec_YshowerX, m_PreRec_ZshowerX;
	FloatVec m_PreRec_EshowerY, m_PreRec_XshowerY, m_PreRec_YshowerY, m_PreRec_ZshowerY;
	Int_t    m_PreRec_NshowerX, m_PreRec_NshowerY, m_PreRec_NclusterX, m_PreRec_NclusterY;


	dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
	dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

	Gaudi::Property<float> m_scale{ this, "Scale", 1 };

  // Input collections
  DataHandle<edm4hep::SimCalorimeterHitCollection> r_SimCaloCol{"SimCaloCol", Gaudi::DataHandle::Reader, this};
	mutable Gaudi::Property<std::string> _filename{this, "OutFileName", "testout.root", "Output file name"};
	mutable Gaudi::Property<float> _seed{this, 	"Seed", 2131, "Random Seed"};
	mutable Gaudi::Property<int>  _Debug{this, 	"Debug", 0, "Debug level"};
	mutable Gaudi::Property<float> _Eth {this, 	"EnergyThreshold", 0.001, "Energy Threshold (/GeV)"};
	mutable Gaudi::Property<float> r_cali{this, "CalibrECAL", 1, "Calibration coefficients for ECAL"};
	mutable Gaudi::Property<float> Latt{this, 	"AttenuationLength", 7000, "Crystal Attenuation Length(mm)"};
	mutable Gaudi::Property<float> Tres{this, 	"TimeResolution", 0.1, "Crystal time resolution in one side (ns)"};
	mutable Gaudi::Property<float> nMat{this, 	"MatRefractive", 2.15, "Material refractive index of crystal"};
	mutable Gaudi::Property<float> Tinit{this, 	"InitalTime", 2, "Start time (ns)"};

	mutable Gaudi::Property<float> _Qthfrac  {this, 	"ChargeThresholdFrac", 0.05, "Charge threshold fraction"};
	mutable Gaudi::Property<float> _Eth_diff {this, 	"MatchingEnergy", 0.30, "Threshold for energy matching (/100%)"};

	mutable Gaudi::Property<float> _Eth_SeedWithNeigh  {this,    "SeedWithNeighThreshold", 0.4,  ""};
	mutable Gaudi::Property<float> _Eth_SeedWithTot    {this,    "SeedWithTotThreshold",   0.15, ""};
	mutable Gaudi::Property<float> _Eth_ShowerWithTot  {this,    "ShowerWithTotThreshold", 0.05, ""};
	mutable Gaudi::Property<float> _Eth_ClusterWithTot {this,    "ClusterWithTotThreshold",0.05, ""};

	mutable Gaudi::Property<float> _chi2Wi_E {this, 	"EnergyChi2Weight", 1, "Weight for energy contribution in chi2 cal."};
	mutable Gaudi::Property<float> _chi2Wi_T {this, 	"TimeChi2Weight", 1, "Weight for time contribution in chi2 cal."};
	mutable Gaudi::Property<float> _th_chi2 {this, 	"Chi2Threshold", 0., "Threshold for chi2 matching"};

  // Output collections
  DataHandle<edm4hep::CalorimeterHitCollection>    w_DigiCaloCol{"DigiCaloCol", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::SimCalorimeterHitCollection>    w_SimCaloTruth{"SimCaloCol", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoCaloAssociationCollection>    w_CaloAssociationCol{"MCRecoCaloAssociationCollection", Gaudi::DataHandle::Writer, this};
};

#endif
