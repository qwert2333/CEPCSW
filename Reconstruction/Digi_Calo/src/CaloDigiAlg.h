#ifndef Calo_DIGI_ALG_H
#define Calo_DIGI_ALG_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "edm4hep/MutableCaloHitContribution.h"
#include "edm4hep/MutableSimCalorimeterHit.h"
#include "edm4hep/MCRecoCaloParticleAssociationCollection.h"
#include "edm4hep/MCRecoParticleAssociationCollection.h"
#include "edm4hep/MCParticle.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include "DetInterface/IGeomSvc.h"




class CaloDigiAlg : public GaudiAlgorithm
{
 
public:
 
  CaloDigiAlg(const std::string& name, ISvcLocator* svcLoc);
 
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
 
protected:

  SmartIF<IGeomSvc> m_geosvc;
  typedef std::vector<float> FloatVec;

  int _nEvt ;
  float m_length;
  dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
 
  Gaudi::Property<float> m_scale{ this, "Scale", 1 };
  Gaudi::Property< std::string > name_MCParticleCol{ this, "MCParticleCollection", "MCParticle" };
  Gaudi::Property< std::vector<std::string> > name_EcalHits{ this, "ECalCaloHitCollections", {"ECALBarrel"} };
  Gaudi::Property< std::vector<std::string> > name_EcalReadout{ this, "ECalReadOutNames", {"EcalBarrelCollection"} };
  Gaudi::Property< std::vector<std::string> > name_EcalMCPAssociation{ this, "ECalMCPAssociationName", {"ECALBarrelParticleAssoCol"} };

  // Input collections
  DataHandle<edm4hep::SimCalorimeterHitCollection> r_SimCaloCol{"SimCaloCol", Gaudi::DataHandle::Reader, this};
//  DataHandle<edm4hep::SimCalorimeterHitCollection> r_MCPCol{"SimCaloCol", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::MCRecoCaloParticleAssociationCollection> r_MCRCaloPAssCol{"ECALBarrelParticleAssoCol", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::MCRecoCaloAssociationCollection>    r_CaloAssCol{"CaloAssCol", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::CalorimeterHitCollection>    r_DigiECaloCol{"DigiECaloCol", Gaudi::DataHandle::Reader, this};

  // Output collections
  DataHandle<edm4hep::CalorimeterHitCollection>    w_DigiCaloCol{"DigiCaloCol", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoCaloAssociationCollection>    w_CaloAssociationCol{"MCRecoCaloAssociationCollection", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoParticleAssociationCollection>  w_MCRecoParticleAssociationCollection {"MCRecoParticleAssociation",Gaudi::DataHandle::Writer, this};

  //For Ana
  Gaudi::Property<bool>  m_WriteAna {this, "WriteAna", false, "Write Ntuples for analysis"};
  Gaudi::Property<std::string> m_filename{this, "AnaFileName", "testout.root", "Output file name"};
};  

#endif
