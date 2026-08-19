#Sun Mar 10 20:54:08 2024"""Automatically generated. DO NOT EDIT please"""
import sys
if sys.version_info >= (3,):
    # Python 2 compatibility
    long = int
from GaudiKernel.DataHandle import DataHandle
from GaudiKernel.Proxy.Configurable import *

class TrackInspectAlg( ConfigurableAlgorithm ) :
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
    'RegisterForContextService' : False, # bool
    'Cardinality' : 1, # int
    'NeededResources' : [  ], # list
    'Blocking' : False, # bool
    'FilterCircularDependencies' : True, # bool
    'Weight' : 0.50000000, # float
    'useTPC' : True, # bool
    'useVXD' : True, # bool
    'useSIT' : True, # bool
    'useSET' : True, # bool
    'useFTD' : True, # bool
    'TrackCollection' : DataHandle('InputTrackCollection', 'R', 'DataWrapper<edm4hep::TrackCollection>'), # DataHandle
    'MCParticleCollection' : DataHandle('InputMCParticleCollection', 'R', 'DataWrapper<edm4hep::MCParticleCollection>'), # DataHandle
    'TPCTrackerHitRelations' : DataHandle('TPCTrackerHitRelations', 'R', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
    'VXDTrackerHitRelations' : DataHandle('VXDTrackerHitRelations', 'R', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
    'SITTrackerHitRelations' : DataHandle('SITTrackerHitRelations', 'R', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
    'SETTrackerHitRelations' : DataHandle('SETTrackerHitRelations', 'R', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
    'FTDTrackerHitRelations' : DataHandle('FTDTrackerHitRelations', 'R', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
    'useTPC' : True, # bool
    'useVXD' : True, # bool
    'useSIT' : True, # bool
    'useSET' : True, # bool
    'useFTD' : True, # bool
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
    'Weight' : """  [TrackInspectAlg] """,
    'useTPC' : """ flag whether to use TPC hits [TrackInspectAlg] """,
    'useVXD' : """ flag whether to use VXD hits [TrackInspectAlg] """,
    'useSIT' : """ flag whether to use SIT hits [TrackInspectAlg] """,
    'useSET' : """ flag whether to use SET hits [TrackInspectAlg] """,
    'useFTD' : """ flag whether to use FTD hits [TrackInspectAlg] """,
    'TrackCollection' : """ Handle of the Input Track collection [unknown owner type] """,
    'MCParticleCollection' : """ Handle of the Input MC particle collection [unknown owner type] """,
    'TPCTrackerHitRelations' : """ Handle of the TPC Tracker Hit relations [unknown owner type] """,
    'VXDTrackerHitRelations' : """ Handle of the TPC Tracker Hit relations [unknown owner type] """,
    'SITTrackerHitRelations' : """ Handle of the TPC Tracker Hit relations [unknown owner type] """,
    'SETTrackerHitRelations' : """ Handle of the TPC Tracker Hit relations [unknown owner type] """,
    'FTDTrackerHitRelations' : """ Handle of the TPC Tracker Hit relations [unknown owner type] """,
    'useTPC' : """ flag whether to use TPC hits [TrackInspectAlg] """,
    'useVXD' : """ flag whether to use VXD hits [TrackInspectAlg] """,
    'useSIT' : """ flag whether to use SIT hits [TrackInspectAlg] """,
    'useSET' : """ flag whether to use SET hits [TrackInspectAlg] """,
    'useFTD' : """ flag whether to use FTD hits [TrackInspectAlg] """,
  }
  __declaration_location__ = 'TrackInspectAlg.cpp:53'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(TrackInspectAlg, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'TrackInspect'
  def getType( self ):
      return 'TrackInspectAlg'
  pass # class TrackInspectAlg
