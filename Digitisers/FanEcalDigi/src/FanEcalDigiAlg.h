#ifndef FAN_ECAL_DIGI_ALG_H
#define FAN_ECAL_DIGI_ALG_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
//#include "edm4hep/SimCalorimeterHitConst.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "edm4hep/MCRecoCaloParticleAssociationCollection.h"
#include "edm4hep/MutableCaloHitContribution.h"
#include "edm4hep/MutableSimCalorimeterHit.h"
#include "edm4hep/CaloHitContributionCollection.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/Segmentations.h> 
#include "DetInterface/IGeomSvc.h"

#include "CaloBar.h"
//#include "CaloCrossBar.h"
//#include "CaloCrossCluster.h"
#include "CaloStep.h"
#include "CaloCluster.h"
//#include "SiPM_CaloCluster.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TString.h"
#include "TH3.h"
#include "TH1.h"
#include <TLorentzVector.h>

#define C 299.79  // unit: mm/ns
#define PI 3.141592653

class FanEcalDigiAlg : public GaudiAlgorithm
{
 
public:
 
  FanEcalDigiAlg(const std::string& name, ISvcLocator* svcLoc);
 
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
//	edm4hep::MutableSimCalorimeterHit find(edm4hep::MutableSimCalorimeterHit& m_col, dd4hep::Position& pos);
	edm4hep::MutableSimCalorimeterHit find(std::vector<edm4hep::MutableSimCalorimeterHit>& m_col, unsigned long long& cellid);
 
  StatusCode HitClassification(const edm4hep::SimCalorimeterHitCollection& m_col, std::map<std::tuple<int, double, double, double, double>, std::vector<edm4hep::SimCalorimeterHit>>& HitCollectionMap);
  StatusCode MergeSameParticleHits(const std::map<std::tuple<int, double, double, double, double>, std::vector<edm4hep::SimCalorimeterHit>>& hitsToMerge, std::map<std::tuple<int, double, double, double, double>, std::vector<edm4hep::SimCalorimeterHit>>& mergedHits);

  StatusCode NeighborClustering(std::vector<CaloBar>& m_barVec, std::vector<CaloCluster>& m_clusVec ); 
  StatusCode ToCrossCluster(vector<CaloCrossBar>& m_bars, CaloCrossCluster& m_cross_clus);
  StatusCode MergeCrossBars(vector<CaloCrossBar>& m_bars, vector<CaloCrossBar>& n_bars);
  void Clear();
  double getPhi(double x, double y);
  double getR(double x, double y, double z);
  double getTheta(double x, double y, double z);
  double getHiggsMass(vector<CaloCluster>& m_evtClus);
  StatusCode CombineInCluster(std::vector<CaloCluster>& m_clus);
  StatusCode CombineDiffCluster(CaloCluster& m_odd_clus, CaloCluster& m_evn_clus); 
  StatusCode MergeCluster(std::vector<CaloCluster>& m_oddVec, std::vector<CaloCluster>& m_evnVec, std::vector<CaloCluster>& m_evtVec);
  StatusCode CLUEClustering(std::vector<CaloBar>& m_barVec, std::vector<CaloCluster>& m_clus, string input, double dc, double deltao, double deltac, double rho);

protected:

  SmartIF<IGeomSvc> m_geosvc;
  typedef std::vector<float> FloatVec;
  typedef std::vector<int> IntVec;
  typedef std::map<const edm4hep::MCParticle, float> MCParticleToEnergyWeightMap;

	int _nEvt ;
	int recEvt;
	int rec_Evt;
	float m_length;
	TRandom3 rndm;
	TFile* m_wfile;
	TTree* t_SimStep;
	TTree* t_SimBar;
	TTree* t_ClusBar;
	TTree* t_MCPdata;
	TTree* t_MCHitInfo;
	
	FloatVec m_MChit_primary_px, m_MChit_primary_py, m_MChit_primary_pz, m_MChit_crystal, m_MChit_module, m_MChit_x, m_MChit_y, m_MChit_z, m_MChit_E, m_MChit_shower_pdg, m_MChit_primary_E, m_MChit_primary_pdg;	
	FloatVec m_step_x, m_step_y, m_step_z, m_step_E, m_step_T1, m_step_T2, m_stepBar_x, m_stepBar_y, m_stepBar_z;
	FloatVec m_simBar_x, m_simBar_y, m_simBar_z, m_simBar_T1, m_simBar_T2, m_simBar_Q1, m_simBar_Q2, m_simBar_crystal, m_simBar_module;
	IntVec m_simBar_index;
	IntVec m_MChit_index;

	//==============bars in cluster====================
	FloatVec m_clusBar_x, m_clusBar_y, m_clusBar_z, m_clusBar_NR, m_clusBar_Nphi,m_clusBar_NE, m_clusBar_Nz, m_clusBar_T1, m_clusBar_T2, m_clusBar_Q1, m_clusBar_Q2, m_clusBar_crystal, m_clusBar_module, m_InParticle_E_MC, m_InParticle_px_MC, m_InParticle_py_MC, m_InParticle_pz_MC, m_InParticle_x_MC, m_InParticle_y_MC, m_InParticle_z_MC;
        IntVec m_nclusBar, m_ncrossBar, m_clusBar_index;

	FloatVec m_oclusBar_Q1, m_oclusBar_Q2, m_oclusBar_module, m_oclusBar_crystal, m_eclusBar_Q1, m_eclusBar_Q2, m_eclusBar_module, m_eclusBar_crystal, m_evt_clusBar_Q1, m_evt_clusBar_Q2, m_evt_clusBar_module, m_evt_clusBar_crystal;
	IntVec m_oclusBar_index, m_eclusBar_index, m_evt_clusBar_index;

 	FloatVec m_MCP_pdg, m_MCP_charge, m_MCP_mass, m_MCP_endpoint_x, m_MCP_endpoint_y, m_MCP_endpoint_z, m_MCP_module, m_MCP_crystal; 
	IntVec m_MCP_index;

  int m_Ncluster;
  FloatVec m_SiPMclus_id, m_SiPMclus_E, m_SiPMclus_RStart, m_SiPMclus_aveR,  m_SiPMclus_Z, m_SiPMclus_avephi,m_SiPMclus_Ephi, m_clus_phi, m_clus_Z,m_clus_aphi,m_clus_aphi_E,m_clus_cid, m_clus_cR_start, m_clus_cR, m_clus_cphi,m_clus_cavephi, m_clus_cEphi, m_clus_cE, m_clus_cZ, m_clus_E, m_clus_phi_start,m_clus_chi2,m_clus_alpha,m_clus_beta,m_MCendpoint_x,m_MCendpoint_R,m_MCendpoint_theta,m_MCendpoint_phi,m_MCendpoint_y,m_MCendpoint_z,m_MCmomentum,m_MCmomentum_x,m_MCmomentum_y,m_MCmomentum_z;
  IntVec m_clus_Nbars,m_clus_id, m_odd_clus_index, m_evn_clus_index, m_evt_clus_index, m_cross_clus_index, m_cross_clus_bar_index;


	dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
	dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

	Gaudi::Property<float> m_scale{ this, "Scale", 1 };

  // Input collections
  DataHandle<edm4hep::MCParticleCollection> m_mcParCol{"MCParticle", Gaudi::DataHandle::Reader, this};//add
  DataHandle<edm4hep::SimCalorimeterHitCollection> r_SimCaloCol{"SimCaloCol", Gaudi::DataHandle::Reader, this};
 // DataHandle<edm4hep::CaloHitContributionCollection> r_CaloContCol{"CaloContCol", Gaudi::DataHandle::Reader, this};

  mutable Gaudi::Property<std::string> _readoutName{this, "ReadOutName", "CaloHitsCollection", "Readout name"}; 
  mutable Gaudi::Property<int>   _writeNtuple{this,  "WriteNtuple", 1, "Write ntuple"}; 
//  mutable Gaudi::Property<std::string> _readoutMC{this, "ReadOutNameMC", "MCParticle", "Readout name"};//add
  mutable Gaudi::Property<std::string> _filename{this, "OutFileName", "testout.root", "Output file name"};
  mutable Gaudi::Property<int>   _Nskip{this,  "SkipEvt", 0, "Skip event"};
  mutable Gaudi::Property<float> _seed{this,   "Seed", 2131, "Random Seed"};
  mutable Gaudi::Property<int>  _Debug{this,   "Debug", 0, "Debug level"};
  mutable Gaudi::Property<float> _Eth {this,   "EnergyThreshold", 0.001, "Energy Threshold (/GeV)"};
  mutable Gaudi::Property<float> r_cali{this,  "CalibrECAL", 1, "Calibration coefficients for ECAL"};
  mutable Gaudi::Property<float> Lbar{this, 	"CrystalBarLength", 312, "Crystal Bar Length(mm)"};
  mutable Gaudi::Property<float> Latt{this, 	"AttenuationLength", 7000, "Crystal Attenuation Length(mm)"};
  mutable Gaudi::Property<float> Tres{this, 	"TimeResolution", 0.1, "Crystal time resolution in one side (ns)"};
  mutable Gaudi::Property<float> nMat{this, 	"MatRefractive", 2.15, "Material refractive index of crystal"};
  mutable Gaudi::Property<float> Tinit{this, 	"InitalTime", 2, "Start time (ns)"};
  
  mutable Gaudi::Property<float> _Qthfrac  {this, 	"ChargeThresholdFrac", 0.05, "Charge threshold fraction"};


  //================re DIGI===========================================
    const Double_t fCryMeanLY = 850; // 850 p.e./MIP, 3x3mm SiPM
    const Double_t fMIPEnergy = 18; // 9 MeV / MIP in 1 cm BGO
    const Double_t fADCMean = 15;
    const Double_t fADCSigma = 3;
    const Double_t MIP_Thre = 0.3; // 0.3 MIP
  //==================================================================

  // Output collections
  DataHandle<edm4hep::CalorimeterHitCollection>    w_DigiCaloCol{"DigiCaloCol", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoCaloAssociationCollection>    w_CaloAssociationCol{"ECALBarrelAssoCol", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoCaloParticleAssociationCollection>    w_MCPCaloAssociationCol{"ECALBarrelParticleAssoCol", Gaudi::DataHandle::Writer, this};
};

#endif
