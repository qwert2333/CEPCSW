#Sun Mar 10 20:54:32 2024"""Automatically generated. DO NOT EDIT please"""
import sys
if sys.version_info >= (3,):
    # Python 2 compatibility
    long = int
from GaudiKernel.DataHandle import DataHandle
from GaudiKernel.Proxy.Configurable import *

class RecGenfitAlgDC( ConfigurableAlgorithm ) :
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
    'readout' : 'DriftChamberHitsCollection', # str
    'debug' : 0, # int
    'smearHit' : True, # bool
    'sigmaHit' : 0.110000, # float
    'nSigmaHit' : 5.00000, # float
    'initCovResPos' : 1.0000000, # float
    'initCovResMom' : 0.10000000, # float
    'fitterTyep' : 'DAFRef', # str
    'correctBremsstrahlung' : False, # bool
    'noMaterialEffects' : False, # bool
    'maxIteration' : 20, # int
    'resortHits' : 1, # int
    'bStart' : 100.00000, # float
    'bFinal' : 0.010000000, # float
    'dcCornerCuts' : -999.00000, # float
    'ndfCut' : 1.0000000, # float
    'chi2Cut' : 1000.0000, # float
    'debugPid' : -99, # int
    'useTruthTrack' : True, # bool
    'useTruthHit' : True, # bool
    'genfitHistRootName' : '', # str
    'showDisplay' : False, # bool
    'MCParticle' : DataHandle('MCParticle', 'R', 'DataWrapper<edm4hep::MCParticleCollection>'), # DataHandle
    'DriftChamberHitsCollection' : DataHandle('DriftChamberHitsCollection', 'R', 'DataWrapper<edm4hep::SimTrackerHitCollection>'), # DataHandle
    'DigiDCHitCollection' : DataHandle('DigiDCHitCollection', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'DCTrackCollection' : DataHandle('DCTrackCollection', 'R', 'DataWrapper<edm4hep::TrackCollection>'), # DataHandle
    'DCHitAssociationCollection' : DataHandle('DCHitAssociationCollection', 'R', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
    'DCRecParticleCollection' : DataHandle('DCRecParticleCollection', 'W', 'DataWrapper<edm4hep::ReconstructedParticleCollection>'), # DataHandle
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
    'readout' : """  [RecGenfitAlgDC] """,
    'debug' : """  [RecGenfitAlgDC] """,
    'smearHit' : """  [RecGenfitAlgDC] """,
    'sigmaHit' : """  [RecGenfitAlgDC] """,
    'nSigmaHit' : """  [RecGenfitAlgDC] """,
    'initCovResPos' : """  [RecGenfitAlgDC] """,
    'initCovResMom' : """  [RecGenfitAlgDC] """,
    'fitterTyep' : """  [RecGenfitAlgDC] """,
    'correctBremsstrahlung' : """  [RecGenfitAlgDC] """,
    'noMaterialEffects' : """  [RecGenfitAlgDC] """,
    'maxIteration' : """  [RecGenfitAlgDC] """,
    'resortHits' : """  [RecGenfitAlgDC] """,
    'bStart' : """  [RecGenfitAlgDC] """,
    'bFinal' : """  [RecGenfitAlgDC] """,
    'dcCornerCuts' : """  [RecGenfitAlgDC] """,
    'ndfCut' : """  [RecGenfitAlgDC] """,
    'chi2Cut' : """  [RecGenfitAlgDC] """,
    'debugPid' : """  [RecGenfitAlgDC] """,
    'useTruthTrack' : """  [RecGenfitAlgDC] """,
    'useTruthHit' : """  [RecGenfitAlgDC] """,
    'genfitHistRootName' : """  [RecGenfitAlgDC] """,
    'showDisplay' : """  [RecGenfitAlgDC] """,
    'MCParticle' : """ Handle of the input MCParticle collection [unknown owner type] """,
    'DriftChamberHitsCollection' : """ Handle of the input SimHit collection [unknown owner type] """,
    'DigiDCHitCollection' : """ Handle of digi DCHit collection [unknown owner type] """,
    'DCTrackCollection' : """ Handle of DC track collection [unknown owner type] """,
    'DCHitAssociationCollection' : """ Handle of simTrackerHit and TrackerHit association collection [unknown owner type] """,
    'DCRecParticleCollection' : """ Handle of drift chamber reconstructed particle collection [unknown owner type] """,
  }
  __declaration_location__ = 'RecGenfitAlgDC.cpp:37'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(RecGenfitAlgDC, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'RecGenfitAlg'
  def getType( self ):
      return 'RecGenfitAlgDC'
  pass # class RecGenfitAlgDC
