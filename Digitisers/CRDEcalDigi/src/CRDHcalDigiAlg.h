#ifndef CRD_HCAL_DIGI_ALG_H 
#define CRD_HCAL_DIGI_ALG_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/MutableCaloHitContribution.h"
#include "edm4hep/MutableSimCalorimeterHit.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "HitStep.h"
#include "CaloHit.h"

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

#include <cstdlib>
#include "time.h"
#include <TTimeStamp.h> 
#include <ctime>

#define C 299.79  // unit: mm/ns
#define PI 3.141592653

class CRDHcalDigiAlg : public GaudiAlgorithm
{
 
public:
 
  CRDHcalDigiAlg(const std::string& name, ISvcLocator* svcLoc);
 
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


  StatusCode MergeHits(const edm4hep::SimCalorimeterHitCollection& m_col, std::vector<edm4hep::SimCalorimeterHit>& m_hits);  

	edm4hep::MutableSimCalorimeterHit find(const std::vector<edm4hep::MutableSimCalorimeterHit>& m_col, unsigned long long& cellid) const;

	void Clear();

protected:

  SmartIF<IGeomSvc> m_geosvc;
  //SmartIF<ICRDEcalSvc> m_edmsvc;
  typedef std::vector<float> FloatVec;

	int _nEvt ;
	TRandom3 rndm;
	TFile* m_wfile;
	TTree* t_SimCont;
	TTree* t_simHit;

	FloatVec m_step_x, m_step_y, m_step_z, m_step_E;
	FloatVec m_simHit_x, m_simHit_y, m_simHit_z, m_simHit_E, m_simHit_slice, m_simHit_layer, m_simHit_tower, m_simHit_stave, m_simHit_module;
  std::vector<unsigned long long> m_simHit_cellID;
  std::vector<int> m_simHit_steps;


	dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
	dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

	Gaudi::Property<float> m_scale{ this, "Scale", 1 };

  // Input collections
  DataHandle<edm4hep::SimCalorimeterHitCollection> r_SimCaloCol{"SimCaloCol", Gaudi::DataHandle::Reader, this};
  mutable Gaudi::Property<std::string> _readoutName{this, "ReadOutName", "CaloHitsCollection", "Readout name"};
  mutable Gaudi::Property<std::string> _filename{this, "OutFileName", "testout.root", "Output file name"};

  //Input parameters
  mutable Gaudi::Property<int>   _Nskip{this,  "SkipEvt", 0, "Skip event"};
  mutable Gaudi::Property<float> _seed{this,   "Seed", 2131, "Random Seed"};
  mutable Gaudi::Property<int>  _Debug{this,   "Debug", 0, "Debug level"};
  mutable Gaudi::Property<float> _Eth {this,   "EnergyThreshold", 0.001, "Energy Threshold (/GeV)"};
  // mutable Gaudi::Property<float> r_cali{this,  "CalibrECAL", 1, "Calibration coefficients for ECAL"};
  
  // mutable Gaudi::Property<float> _Qthfrac  {this, 	"ChargeThresholdFrac", 0.05, "Charge threshold fraction"};


  // Output collections
  DataHandle<edm4hep::CalorimeterHitCollection>    w_DigiCaloCol{"DigiCaloCol", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoCaloAssociationCollection>    w_CaloAssociationCol{"MCRecoCaloAssociationCollection", Gaudi::DataHandle::Writer, this};
};

#endif
