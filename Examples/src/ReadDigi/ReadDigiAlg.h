#ifndef READ_DIGI_H
#define READ_DIGI_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"

#include "DetInterface/IGeomSvc.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"

#include "UTIL/ILDConf.h"
#include "UTIL/CellIDEncoder.h"
#include <DDRec/DetectorData.h>
#include <DD4hep/Segmentations.h>
#include "TFile.h"
#include "TTree.h"

class ReadDigiAlg : public GaudiAlgorithm
{

public :

  ReadDigiAlg(const std::string& name, ISvcLocator* svcLoc);

  virtual StatusCode initialize();
  virtual StatusCode execute();
  virtual StatusCode finalize();

  StatusCode Clear(); 

private :

  DataHandle<edm4hep::MCParticleCollection>       m_MCParticleCol{"MCParticle", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::CalorimeterHitCollection>   m_HCalBarrelHitCol{"HCALBarrel", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::CalorimeterHitCollection>   m_HCalEndcapHitCol{"HCALEndcap", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::CalorimeterHitCollection>   m_ECalBarrelHitCol{"ECALBarrel", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::CalorimeterHitCollection>   m_ECalEndcapHitCol{"ECALEndcap", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::MCRecoCaloAssociationCollection>  m_MCRecoCaloAssociationCol{"MCRecoCaloAssociationCollection", Gaudi::DataHandle::Reader, this};


  mutable Gaudi::Property<std::string> _filename{this, "OutFileName", "testout.root", "Output file name"};

  typedef std::map<const edm4hep::MCParticle, float> MCParticleToEnergyWeightMap;
  typedef std::vector<float> FloatVec;
  typedef std::vector<int>   IntVec;

  int _nEvt; 
  TFile *m_wfile;
  TTree *m_wtree; 

  //MCParticle
  int N_MCP;
  FloatVec MCP_px, MCP_py, MCP_pz, MCP_E;
  IntVec MCP_pdgid, MCP_gStatus; 

  //ECal
  int Nhit_EcalB;
  IntVec EcalBHit_MCtag, EcalBHit_MCpid; 
  FloatVec EcalBHit_x, EcalBHit_y, EcalBHit_z, EcalBHit_E, EcalBHit_T;

  int Nhit_EcalE;
  IntVec EcalEHit_MCtag, EcalEHit_MCpid;
  FloatVec EcalEHit_x, EcalEHit_y, EcalEHit_z, EcalEHit_E, EcalEHit_T;

  //HCal
  int Nhit_HcalB; 
  IntVec HcalBHit_MCtag, HcalBHit_MCpid; 
  FloatVec HcalBHit_x, HcalBHit_y, HcalBHit_z, HcalBHit_E, HcalBHit_T;

  int Nhit_HcalE;
  IntVec HcalEHit_MCtag, HcalEHit_MCpid;
  FloatVec HcalEHit_x, HcalEHit_y, HcalEHit_z, HcalEHit_E, HcalEHit_T;

  

};

#endif
