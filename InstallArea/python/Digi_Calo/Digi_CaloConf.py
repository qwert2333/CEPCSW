#Sun Mar 10 20:53:59 2024"""Automatically generated. DO NOT EDIT please"""
import sys
if sys.version_info >= (3,):
    # Python 2 compatibility
    long = int
from GaudiKernel.DataHandle import DataHandle
from GaudiKernel.Proxy.Configurable import *

class CaloDigiAlg( ConfigurableAlgorithm ) :
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
    'Scale' : 1.00000, # float
    'MCParticleCollection' : 'MCParticle', # str
    'ECalCaloHitCollections' : [ 'ECALBarrel' ], # list
    'ECalReadOutNames' : [ 'EcalBarrelCollection' ], # list
    'ECalMCPAssociationName' : [ 'ECALBarrelParticleAssoCol' ], # list
    'WriteAna' : False, # bool
    'AnaFileName' : 'testout.root', # str
    'SimCaloHitCollection' : DataHandle('SimCaloCol', 'R', 'DataWrapper<edm4hep::SimCalorimeterHitCollection>'), # DataHandle
    'ECALBarrelParticleAssoCol' : DataHandle('ECALBarrelParticleAssoCol', 'R', 'DataWrapper<edm4hep::MCRecoCaloParticleAssociationCollection>'), # DataHandle
    'CaloHitCollection' : DataHandle('DigiCaloCol', 'W', 'DataWrapper<edm4hep::CalorimeterHitCollection>'), # DataHandle
    'CaloAssociationCollection' : DataHandle('MCRecoCaloAssociationCollection', 'W', 'DataWrapper<edm4hep::MCRecoCaloAssociationCollection>'), # DataHandle
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
    'Scale' : """  [CaloDigiAlg] """,
    'MCParticleCollection' : """  [CaloDigiAlg] """,
    'ECalCaloHitCollections' : """  [CaloDigiAlg] """,
    'ECalReadOutNames' : """  [CaloDigiAlg] """,
    'ECalMCPAssociationName' : """  [CaloDigiAlg] """,
    'WriteAna' : """ Write Ntuples for analysis [CaloDigiAlg] """,
    'AnaFileName' : """ Output file name [CaloDigiAlg] """,
    'SimCaloHitCollection' : """ Handle of the Input SimCaloHit collection [unknown owner type] """,
    'ECALBarrelParticleAssoCol' : """ Handle of CaloAssociation collection [unknown owner type] """,
    'CaloHitCollection' : """ Handle of Digi CaloHit collection [unknown owner type] """,
    'CaloAssociationCollection' : """ Handle of CaloAssociation collection [unknown owner type] """,
  }
  __declaration_location__ = 'CaloDigiAlg.cpp:18'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(CaloDigiAlg, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'Digi_Calo'
  def getType( self ):
      return 'CaloDigiAlg'
  pass # class CaloDigiAlg
