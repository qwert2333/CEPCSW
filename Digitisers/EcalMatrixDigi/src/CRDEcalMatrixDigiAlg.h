#ifndef CRD_ECAL_DIGI_ALG_H
#define CRD_ECAL_DIGI_ALG_H

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
#include "TFile.h"
#include "TString.h"
#include "TRandom3.h"
#include "TH3.h"
#include "TH1.h"

#define C 299.79  // unit: mm/ns
#define PI 3.141592653

class CRDEcalMatrixDigiAlg : public GaudiAlgorithm
{
 
public:
 
  CRDEcalMatrixDigiAlg(const std::string& name, ISvcLocator* svcLoc);
 
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
 
  void clear();

  class DigiBar {
  public:
     unsigned long long cellID = 0;
     dd4hep::Position position;
     double Q1=0;      // Q in left readout
     double Q2=0;      // Q in right readout;
     double T1=999;    // T in left readout;
     double T2=999;    // T in right readout;
  };
	class StepDigiOut {
	public: 
		double Q;
		double T;
		inline bool operator < (const StepDigiOut &x) const {
			return T<x.T ;
		}
	};

	std::vector<edm4hep::SimCalorimeterHit> MergeHits(const edm4hep::SimCalorimeterHitCollection* m_col);
	double GetBarLength(CRDEcalMatrixDigiAlg::DigiBar bar);
	dd4hep::Position GetCellPos(dd4hep::Position pos, CRDEcalMatrixDigiAlg::DigiBar bar);
	edm4hep::SimCalorimeterHit find(edm4hep::SimCalorimeterHitCollection* m_col, dd4hep::Position pos);
	edm4hep::SimCalorimeterHit find(std::vector<edm4hep::SimCalorimeterHit> m_col, unsigned long long cellid);
	std::vector<edm4hep::CalorimeterHit> CreateDigiHits(std::vector<CRDEcalMatrixDigiAlg::DigiBar> m_block);
	std::vector<edm4hep::CalorimeterHit> DigiHitsWithTime(std::vector<CRDEcalMatrixDigiAlg::DigiBar> m_block);
	void Clear();


protected:

  SmartIF<IGeomSvc> m_geosvc;
  typedef std::vector<float> FloatVec;

	TRandom3 rndm;
	int _nEvt ;
	float m_length;
	TFile* m_wfile;
	TTree* t_SimBar;
	TTree* t_SimTruth;
	TTree *t_SimCont;
	
	FloatVec m_step_x, m_step_y, m_step_z, m_step_E, m_step_T1, m_step_T2, m_stepBar_x, m_stepBar_y, m_stepBar_z;
	FloatVec m_simBar_x, m_simBar_y, m_simBar_z, m_simBar_T1, m_simBar_T2, m_simBar_Q1, m_simBar_Q2, m_simBar_dlayer, m_simBar_part, m_simBar_block, m_simBar_slayer, m_simBar_Nstep;
	FloatVec m_simTruth_x, m_simTruth_y, m_simTruth_z, m_simTruth_E, m_simTruth_dlayer, m_simTruth_slayer;

	dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
	dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

	typedef std::vector<CRDEcalMatrixDigiAlg::DigiBar> BarCol;

	Gaudi::Property<float> m_scale{ this, "Scale", 1 };

  // Input collections
  DataHandle<edm4hep::SimCalorimeterHitCollection> r_SimCaloCol{"SimCaloCol", Gaudi::DataHandle::Reader, this};
	mutable Gaudi::Property<float> _seed{this, 	"Seed", 2131, "Random Seed"};
	mutable Gaudi::Property<int>  _Debug{this, 	"Debug", 0, "Debug level"};
	mutable Gaudi::Property<float> _Eth {this, 	"EnergyThreshold", 0.001, "Energy Threshold (/GeV)"};
	mutable Gaudi::Property<float> r_cali{this, "CalibrECAL", 1, "Calibration coefficients for ECAL"};
	mutable Gaudi::Property<float> Latt{this, 	"AttenuationLength", 7000, "Crystal Attenuation Length(mm)"};
	mutable Gaudi::Property<float> Tres{this, 	"TimeResolution", 0.1, "Crystal time resolution in one side (ns)"};
	mutable Gaudi::Property<float> nMat{this, 	"MatRefractive", 2.15, "Material refractive index of crystal"};
	mutable Gaudi::Property<float> _Qthfrac  {this, 	"ChargeThresholdFrac", 0.05, "Charge threshold fraction"};
	mutable Gaudi::Property<float> _DeltaZth {this, 	"PositionThreshold", 3, "Position threshold for cross-location(/sigma)"};
	mutable Gaudi::Property<std::string> _filename{this, "OutFileName", "testout.root", "Output file name"};

  // Output collections
  DataHandle<edm4hep::CalorimeterHitCollection>    w_DigiCaloCol{"DigiCaloCol", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::SimCalorimeterHitCollection>    w_SimCaloTruth{"SimCaloCol", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoCaloAssociationCollection>    w_CaloAssociationCol{"MCRecoCaloAssociationCollection", Gaudi::DataHandle::Writer, this};
};

#endif
