#Sun Mar 10 20:54:49 2024"""Automatically generated. DO NOT EDIT please"""
import sys
if sys.version_info >= (3,):
    # Python 2 compatibility
    long = int
from GaudiKernel.DataHandle import DataHandle
from GaudiKernel.Proxy.Configurable import *

class PandoraPlusPFAlg( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : [], # list
    'ExtraOutputs' : [], # list
    'OutputLevel' : 0, # int
    'Enable' : True, # bool
    'ErrorMax' : 1, # int
    'AuditAlgorithms' : False, # bool
    'AuditInitialize' : False, # bool
    'AuditReinitialize' : False, # bool
    'AuditRestart' : False, # bool
    'AuditExecute' : False, # bool
    'AuditFinalize' : False, # bool
    'AuditStart' : False, # bool
    'AuditStop' : False, # bool
    'Timeline' : True, # bool
    'MonitorService' : 'MonitorSvc', # str
    'RegisterForContextService' : True, # bool
    'Cardinality' : 1, # int
    'NeededResources' : [  ], # list
    'Blocking' : False, # bool
    'FilterCircularDependencies' : True, # bool
    'RootInTES' : '', # str
    'ErrorsPrint' : True, # bool
    'PropertiesPrint' : False, # bool
    'TypePrint' : True, # bool
    'Context' : '', # str
    'CounterList' : [ '.*' ], # list
    'VetoObjects' : [  ], # list
    'RequireObjects' : [  ], # list
    'MCParticleCollection' : 'MCParticle', # str
    'MCRecoTrackParticleAssociationCollection' : 'MarlinTrkAssociation', # str
    'TrackCollections' : [ 'MarlinTrkTracks' ], # list
    'ECalCaloHitCollections' : [ 'ECALBarrel' ], # list
    'ECalReadOutNames' : [ 'EcalBarrelCollection' ], # list
    'ECalMCPAssociationName' : [ 'ECALBarrelParticleAssoCol' ], # list
    'HCalCaloHitCollections' : [ 'HCALBarrel' ], # list
    'HCalReadOutNames' : [ 'HcalBarrelCollection' ], # list
    'HCalMCPAssociationName' : [ 'HCALBarrelParticleAssoCol' ], # list
    'BField' : 3.00000, # float
    'Seed' : 2131.00, # float
    'Debug' : 0, # int
    'SkipEvt' : 0, # int
    'EcalType' : 'BarEcal', # str
    'UseTruthTrack' : True, # bool
    'AlgList' : [  ], # list
    'AlgParNames' : [  ], # list
    'AlgParTypes' : [  ], # list
    'AlgParValues' : [  ], # list
    'WriteAna' : False, # bool
    'AnaFileName' : 'testout.root', # str
    'RecCaloHitCollection' : DataHandle('RecECALBarrel', 'W', 'DataWrapper<edm4hep::CalorimeterHitCollection>'), # DataHandle
    'ClusterCollection' : DataHandle('PandoraClusters', 'W', 'DataWrapper<edm4hep::ClusterCollection>'), # DataHandle
    'RecoPFOCollection' : DataHandle('PandoraPFOs', 'W', 'DataWrapper<edm4hep::ReconstructedParticleCollection>'), # DataHandle
    'RecoVtxCollection' : DataHandle('PandoraPFANewStartVertices', 'W', 'DataWrapper<edm4hep::VertexCollection>'), # DataHandle
    'MCRecoPFOAssociationCollection' : DataHandle('pfoMCRecoParticleAssociation', 'W', 'DataWrapper<edm4hep::MCRecoParticleAssociationCollection>'), # DataHandle
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'ErrorMax' : """ [[deprecated]] max number of errors [Gaudi::Algorithm] """,
    'AuditAlgorithms' : """ [[deprecated]] unused [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Blocking' : """ if algorithm invokes CPU-blocking system calls (offloads computations to accelerators or quantum processors, performs disk or network I/O, is bound by resource synchronization, etc) [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'RootInTES' : """ note: overridden by parent settings [FixTESPath<Algorithm>] """,
    'ErrorsPrint' : """ print the statistics of errors/warnings/exceptions [GaudiCommon<Algorithm>] """,
    'PropertiesPrint' : """ print the properties of the component [GaudiCommon<Algorithm>] """,
    'TypePrint' : """ add the actual C++ component type into the messages [GaudiCommon<Algorithm>] """,
    'Context' : """ note: overridden by parent settings [GaudiCommon<Algorithm>] """,
    'CounterList' : """ RegEx list, of simple integer counters for CounterSummary [GaudiCommon<Algorithm>] """,
    'VetoObjects' : """ skip execute if one or more of these TES objects exist [GaudiAlgorithm] """,
    'RequireObjects' : """ execute only if one or more of these TES objects exist [GaudiAlgorithm] """,
    'MCParticleCollection' : """  [PandoraPlusPFAlg] """,
    'MCRecoTrackParticleAssociationCollection' : """  [PandoraPlusPFAlg] """,
    'TrackCollections' : """  [PandoraPlusPFAlg] """,
    'ECalCaloHitCollections' : """  [PandoraPlusPFAlg] """,
    'ECalReadOutNames' : """  [PandoraPlusPFAlg] """,
    'ECalMCPAssociationName' : """  [PandoraPlusPFAlg] """,
    'HCalCaloHitCollections' : """  [PandoraPlusPFAlg] """,
    'HCalReadOutNames' : """  [PandoraPlusPFAlg] """,
    'HCalMCPAssociationName' : """  [PandoraPlusPFAlg] """,
    'BField' : """ Magnetic field [PandoraPlusPFAlg] """,
    'Seed' : """ Random Seed [PandoraPlusPFAlg] """,
    'Debug' : """ Debug level [PandoraPlusPFAlg] """,
    'SkipEvt' : """ Skip event [PandoraPlusPFAlg] """,
    'EcalType' : """ ECAL type [PandoraPlusPFAlg] """,
    'UseTruthTrack' : """ Use truth track or reconstructed track [PandoraPlusPFAlg] """,
    'AlgList' : """  [PandoraPlusPFAlg] """,
    'AlgParNames' : """  [PandoraPlusPFAlg] """,
    'AlgParTypes' : """  [PandoraPlusPFAlg] """,
    'AlgParValues' : """  [PandoraPlusPFAlg] """,
    'WriteAna' : """ Write Ntuples for analysis [PandoraPlusPFAlg] """,
    'AnaFileName' : """ Output file name [PandoraPlusPFAlg] """,
    'RecCaloHitCollection' : """ Handle of Reconstructed CaloHit collection [unknown owner type] """,
    'ClusterCollection' : """ Handle of Cluster collection [unknown owner type] """,
    'RecoPFOCollection' : """ Handle of Reconstructed PFO collection [unknown owner type] """,
    'RecoVtxCollection' : """ Handle of Reconstructed vertex collection [unknown owner type] """,
    'MCRecoPFOAssociationCollection' : """ Handle of MC-RecoPFO association collection [unknown owner type] """,
  }
  __declaration_location__ = 'PandoraPlusPFAlg.cpp:23'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(PandoraPlusPFAlg, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'CRDEcalRec'
  def getType( self ):
      return 'PandoraPlusPFAlg'
  pass # class PandoraPlusPFAlg
