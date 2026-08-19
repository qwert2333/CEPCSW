#Sun Mar 10 20:55:22 2024"""Automatically generated. DO NOT EDIT please"""
import sys
if sys.version_info >= (3,):
    # Python 2 compatibility
    long = int
from GaudiKernel.DataHandle import DataHandle
from GaudiKernel.Proxy.Configurable import *

class CylinderDigiAlg( ConfigurableAlgorithm ) :
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
    'ResRPhi' : 0.100000, # float
    'ResZ' : 2.82800, # float
    'InputSimTrackerHitCollection' : DataHandle('DriftChamberHitsCollection', 'R', 'DataWrapper<edm4hep::SimTrackerHitCollection>'), # DataHandle
    'OutputTrackerHitCollection' : DataHandle('DCTrackerHits', 'W', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'TrackerHitAssociationCollection' : DataHandle('DCTrackerHitAssociationCollection', 'W', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
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
    'ResRPhi' : """  [CylinderDigiAlg] """,
    'ResZ' : """  [CylinderDigiAlg] """,
    'InputSimTrackerHitCollection' : """ Handle of the Input SimTrackerHit collection [unknown owner type] """,
    'OutputTrackerHitCollection' : """ Handle of the output TrackerHit collection [unknown owner type] """,
    'TrackerHitAssociationCollection' : """ Handle of the Association collection between SimTrackerHit and TrackerHit [unknown owner type] """,
  }
  __declaration_location__ = 'CylinderDigiAlg.cpp:16'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(CylinderDigiAlg, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'SimpleDigi'
  def getType( self ):
      return 'CylinderDigiAlg'
  pass # class CylinderDigiAlg

class PlanarDigiAlg( ConfigurableAlgorithm ) :
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
    'ResolutionU' : [ 0.00400000 ], # list
    'ResolutionV' : [ 0.00400000 ], # list
    'IsStrip' : False, # bool
    'UsePlanarTag' : True, # bool
    'EnergyThreshold' : 0.00000, # float
    'PullCutToRetry' : 1000.00, # float
    'HeaderCol' : DataHandle('EventHeaderCol', 'R', 'DataWrapper<edm4hep::EventHeaderCollection>'), # DataHandle
    'SimTrackHitCollection' : DataHandle('VXDCollection', 'R', 'DataWrapper<edm4hep::SimTrackerHitCollection>'), # DataHandle
    'TrackerHitCollection' : DataHandle('VXDTrackerHits', 'W', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'TrackerHitAssociationCollection' : DataHandle('VXDTrackerHitRelations', 'W', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
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
    'ResolutionU' : """  [PlanarDigiAlg] """,
    'ResolutionV' : """  [PlanarDigiAlg] """,
    'IsStrip' : """  [PlanarDigiAlg] """,
    'UsePlanarTag' : """  [PlanarDigiAlg] """,
    'EnergyThreshold' : """  [PlanarDigiAlg] """,
    'PullCutToRetry' : """  [PlanarDigiAlg] """,
    'SimTrackHitCollection' : """ Handle of the Input SimTrackerHit collection [unknown owner type] """,
    'TrackerHitCollection' : """ Handle of the TrackerHit output collection [unknown owner type] """,
    'TrackerHitAssociationCollection' : """ Handle of TrackerHit SimTrackHit relation collection [unknown owner type] """,
  }
  __declaration_location__ = 'PlanarDigiAlg.cpp:38'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(PlanarDigiAlg, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'SimpleDigi'
  def getType( self ):
      return 'PlanarDigiAlg'
  pass # class PlanarDigiAlg

class TPCDigiAlg( ConfigurableAlgorithm ) :
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
    'HeaderCol' : DataHandle('EventHeader', 'R', 'DataWrapper<edm4hep::EventHeaderCollection>'), # DataHandle
    'TPCCollection' : DataHandle('TPCCollection', 'R', 'DataWrapper<edm4hep::SimTrackerHitCollection>'), # DataHandle
    'TPCSpacePointCollection' : DataHandle('TPCSpacePointCollection', 'R', 'DataWrapper<edm4hep::SimTrackerHitCollection>'), # DataHandle
    'TPCLowPtCollection' : DataHandle('TPCLowPtCollection', 'R', 'DataWrapper<edm4hep::SimTrackerHitCollection>'), # DataHandle
    'TPCTrackerHitsCol' : DataHandle('TPCTrackerHits', 'W', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'TPCTrackerHitAssCol' : DataHandle('TPCTrackerHitAss', 'W', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
    'UseRawHitsToStoreSimhitPointer' : False, # bool
    'PointResolutionPadPhi' : 0.900000, # float
    'RejectCellID0' : 1, # int
    'PointResolutionRPhi' : 0.0500000, # float
    'DiffusionCoeffRPhi' : 0.0250000, # float
    'N_eff' : 22, # int
    'PointResolutionZ' : 0.400000, # float
    'DiffusionCoeffZ' : 0.0800000, # float
    'HitSortingBinningZ' : 5.00000, # float
    'HitSortingBinningRPhi' : 2.00000, # float
    'DoubleHitResolutionZ' : 5.00000, # float
    'DoubleHitResolutionRPhi' : 2.00000, # float
    'MaxClusterSizeForMerge' : 3, # int
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
    'TPCCollection' : """ The default pad-row based SimTrackerHit collection [unknown owner type] """,
    'TPCSpacePointCollection' : """ Additional space point collection which provides additional guide hits between pad row centers. [unknown owner type] """,
    'TPCLowPtCollection' : """ The LowPt SimTrackerHit collection Produced by Mokka TPC Driver TPC0X [unknown owner type] """,
    'TPCTrackerHitsCol' : """ The output TrackerHit collection [unknown owner type] """,
    'TPCTrackerHitAssCol' : """ Handle of TrackerHit SimTrackHit relation collection [unknown owner type] """,
    'UseRawHitsToStoreSimhitPointer' : """ Store the pointer to the SimTrackerHits in RawHits (deprecated)  [unknown owner type] """,
    'PointResolutionPadPhi' : """ Pad Phi Resolution constant in TPC [unknown owner type] """,
    'RejectCellID0' : """ whether or not to use hits without proper cell ID (pad row) [unknown owner type] """,
    'PointResolutionRPhi' : """ R-Phi Resolution constant in TPC [unknown owner type] """,
    'DiffusionCoeffRPhi' : """ R-Phi Diffusion Coefficent in TPC [unknown owner type] """,
    'N_eff' : """ Number of Effective electrons per pad in TPC [unknown owner type] """,
    'PointResolutionZ' : """ TPC Z Resolution Coefficent independent of diffusion [unknown owner type] """,
    'DiffusionCoeffZ' : """ Z Diffusion Coefficent in TPC [unknown owner type] """,
    'HitSortingBinningZ' : """ Defines spatial slice in Z [unknown owner type] """,
    'HitSortingBinningRPhi' : """ Defines spatial slice in RP [unknown owner type] """,
    'DoubleHitResolutionZ' : """ Defines the minimum distance for two seperable hits in Z [unknown owner type] """,
    'DoubleHitResolutionRPhi' : """ Defines the minimum distance for two seperable hits in RPhi [unknown owner type] """,
    'MaxClusterSizeForMerge' : """ Defines the maximum number of adjacent hits which can be merged [unknown owner type] """,
  }
  __declaration_location__ = 'TPCDigiAlg.cpp:42'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(TPCDigiAlg, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'SimpleDigi'
  def getType( self ):
      return 'TPCDigiAlg'
  pass # class TPCDigiAlg
