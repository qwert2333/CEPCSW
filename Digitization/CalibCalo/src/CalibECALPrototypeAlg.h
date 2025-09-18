#ifndef CalibECALPrototypeAlg_H
#define CalibECALPrototypeAlg_H

#include "k4FWCore/DataHandle.h"
#include <k4FWCore/MetaDataHandle.h>
#include "GaudiKernel/NTuple.h"
#include "GaudiAlg/GaudiAlgorithm.h"

#include <DDRec/DetectorData.h>
#include "DetInterface/IGeomSvc.h"

#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/RawCalorimeterHitCollection.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/Segmentations.h> 
#include "DetInterface/IGeomSvc.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TString.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"

#include <cstdlib>
#include "time.h"
#include <TTimeStamp.h> 
#include <ctime>
#include <iostream>
#include <TMath.h>
#define PI 3.141592653

struct ChannelMapEntry
{
    Int_t Board;
    Int_t Channel;
    Int_t Layer;
    Int_t Index;
    Int_t Side;
};

class CalibECALPrototypeAlg : public GaudiAlgorithm
{
 public:
  
  CalibECALPrototypeAlg(const std::string& name, ISvcLocator* svcLoc);
 
  virtual StatusCode initialize() ;
  virtual StatusCode execute() ; 
  virtual StatusCode finalize() ;


 protected:

  TFile* m_wfile;
  TTree* m_tree;
    
  SmartIF<IGeomSvc> m_geosvc;
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder_barrel;
  dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
  dd4hep::Detector* m_dd4hep;
  std::unique_ptr<podio::GenericParameters> m_metadata;

	std::vector<int> m_cellID;
	std::vector<float> m_pedestalHG, m_pedestalHGsigma, m_pedestalLG, m_pedestalLGsigma, m_MIPHG, m_MIPLG, m_GainFactor; 
	std::map<int, float> map_pedestalHG, map_pedestalHGsigma, map_pedestalLG, map_pedestalLGsigma, map_MIPHG, map_MIPLG, map_GainFactor;
	std::map<int, std::vector<int>> map_HGVec, map_LGVec;
	std::map<int, std::vector<int>> map_HGMipVec, map_LGMipVec;

  Gaudi::Property<bool>  m_dopedestal{ this, "DoPedestalCalib", 1, "Do Pedestal Calibration" };
  Gaudi::Property<bool>  m_dogain{ this, "DoGainCalib", 1, "Do HG-LG inter-Calibration" };
  Gaudi::Property<bool>  m_domip{ this, "DoMIPCalib", 1, "Do MIP Calibration" };
  Gaudi::Property<bool>  m_doana{ this, "DoDataAna", 1, "Calibrat / Analyse data " };
  Gaudi::Property<int>   m_GainLimit{this,  "ChannelHGLimit", 7800, "Write ntuple"};

  Gaudi::Property<std::string> m_metadata_path{this, "ReadinMetaData"};
  Gaudi::Property<std::string> m_channel_map{this, "ChannelMap", "NULL"}; //Temp: used as cell ID decoder
  Gaudi::Property<int> m_skipline_map{this, "MapFileHeaderLine", 8};
  Gaudi::Property<float> MIPCut_OverPedestalSigma{this, "MIPCutOverPedestalSigma", 35.};
  Gaudi::Property<int> HitNumberLimitMIP{this, "HitNumberLimitMIP", 28};
  Gaudi::Property<float> EcalMIPEnergy{this,   "EcalMIPEnergy", 8.918*2, "MIP energy deposit in 2cm BGO (MeV/MIP)"};
  Gaudi::Property<int> GainSwitchPoint{this, "GainSwitchPoint", 7900};
  Gaudi::Property<float> EcalMIPCut{this, "EcalMIPCut", 0.1};


  Gaudi::Property<int>   _writeNtuple{this,  "WriteNtuple", 1, "Write ntuple"};
  Gaudi::Property<std::string> _filename{this, "OutNtupleName", "testout.root", "Output file name"};

  // Input collections
  DataHandle<edm4hep::RawCalorimeterHitCollection> m_inputRawHGHitHdl{"RawHGCaloHit", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::RawCalorimeterHitCollection> m_inputRawLGHitHdl{"RawLGCaloHit", Gaudi::DataHandle::Reader, this};

  // Output collections
  DataHandle<edm4hep::CalorimeterHitCollection> m_outputHitHdl{"CalibCaloHit", Gaudi::DataHandle::Writer, this};

	//MetaData collections
	MetaDataHandle<std::vector<int>> m_cellIDMetaDataHdl{"CellIDVector", Gaudi::DataHandle::Writer};
	MetaDataHandle<std::vector<float>> m_pedestalHGHdl{"PedestalHGVector", Gaudi::DataHandle::Writer};
	MetaDataHandle<std::vector<float>> m_pedestalHGSigmaHdl{"PedestalHGSigmaVector", Gaudi::DataHandle::Writer};
	MetaDataHandle<std::vector<float>> m_pedestalLGHdl{"PedestalLGVector", Gaudi::DataHandle::Writer};
	MetaDataHandle<std::vector<float>> m_pedestalLGSigmaHdl{"PedestalLGSigmaVector", Gaudi::DataHandle::Writer};
	MetaDataHandle<std::vector<float>> m_MIPHGHdl{"MIPHGVector", Gaudi::DataHandle::Writer};
	MetaDataHandle<std::vector<float>> m_MIPLGHdl{"MIPLGVector", Gaudi::DataHandle::Writer};
	MetaDataHandle<std::vector<float>> m_GainFactorHdl{"GainIFVector", Gaudi::DataHandle::Writer};

  SmartIF<IRndmGenSvc> m_randSvc;
  int m_nEvt=0;
  std::map<int, ChannelMapEntry> ChannelMap;
};
#endif



