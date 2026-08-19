#Sun Mar 10 20:55:38 2024"""Automatically generated. DO NOT EDIT please"""
import sys
if sys.version_info >= (3,):
    # Python 2 compatibility
    long = int
from GaudiKernel.DataHandle import DataHandle
from GaudiKernel.Proxy.Configurable import *

class ForwardTrackingAlg( ConfigurableAlgorithm ) :
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
    'Chi2ProbCut' : 0.0050000000, # float
    'HelixFitMax' : 500.00000, # float
    'OverlappingHitsDistMax' : 3.5000000, # float
    'HitsPerTrackMin' : 4, # int
    'BestSubsetFinder' : 'SubsetHopfieldNN', # str
    'TakeBestVersionOfTrack' : True, # bool
    'MaxConnectionsAutomaton' : 100000, # int
    'MaxHitsPerSector' : 1000, # int
    'MultipleScatteringOn' : True, # bool
    'EnergyLossOn' : True, # bool
    'SmoothOn' : False, # bool
    'Criteria' : [ 'Crit2_DeltaPhi' , 'Crit2_DeltaPhi_MV' , 'Crit2_DeltaRho' , 'Crit2_DeltaTheta_MV' , 'Crit2_Distance_MV' , 'Crit2_HelixWithIP' , 'Crit2_RZRatio' , 'Crit2_StraightTrackRatio' , 'Crit3_2DAngle' , 'Crit3_2DAngleTimesR' , 'Crit3_3DAngle' , 'Crit3_3DAngleTimesR' , 'Crit3_ChangeRZRatio' , 'Crit3_IPCircleDist' , 'Crit3_IPCircleDistTimesR' , 'Crit3_NoZigZag_MV' , 'Crit3_PT' , 'Crit3_PT_MV' , 'Crit4_2DAngleChange' , 'Crit4_3DAngleChange' , 'Crit4_3DAngleChangeNormed' , 'Crit4_DistOfCircleCenters' , 'Crit4_DistToExtrapolation' , 'Crit4_NoZigZag' , 'Crit4_PhiZRatioChange' , 'Crit4_RChange' ], # list
    'CriteriaMin' : [  ], # list
    'CriteriaMax' : [  ], # list
    'DumpTime' : True, # bool
    'FTDPixelHitCollection' : DataHandle('FTDPixelTrackerHits', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'FTDSpacePointCollection' : DataHandle('FTDSpacePoints', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'FTDRawHitCollection' : DataHandle('FTDStripTrackerHits', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'ForwardTrackCollection' : DataHandle('ForwardTracks', 'W', 'DataWrapper<edm4hep::TrackCollection>'), # DataHandle
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
    'Chi2ProbCut' : """  [ForwardTrackingAlg] """,
    'HelixFitMax' : """  [ForwardTrackingAlg] """,
    'OverlappingHitsDistMax' : """  [ForwardTrackingAlg] """,
    'HitsPerTrackMin' : """  [ForwardTrackingAlg] """,
    'BestSubsetFinder' : """  [ForwardTrackingAlg] """,
    'TakeBestVersionOfTrack' : """  [ForwardTrackingAlg] """,
    'MaxConnectionsAutomaton' : """  [ForwardTrackingAlg] """,
    'MaxHitsPerSector' : """  [ForwardTrackingAlg] """,
    'MultipleScatteringOn' : """  [ForwardTrackingAlg] """,
    'EnergyLossOn' : """  [ForwardTrackingAlg] """,
    'SmoothOn' : """  [ForwardTrackingAlg] """,
    'Criteria' : """  [ForwardTrackingAlg] """,
    'CriteriaMin' : """  [ForwardTrackingAlg] """,
    'CriteriaMax' : """  [ForwardTrackingAlg] """,
    'DumpTime' : """  [ForwardTrackingAlg] """,
    'FTDPixelHitCollection' : """ Handle of the Input FTD TrackerHits collection [unknown owner type] """,
    'FTDSpacePointCollection' : """ Handle of the Input FTD SpacePoints collection [unknown owner type] """,
    'FTDRawHitCollection' : """ Handle of the FTD SpacePoints raw hit collection collection [unknown owner type] """,
    'ForwardTrackCollection' : """ Handle of the ForwarTrack output collection [unknown owner type] """,
  }
  __declaration_location__ = 'ForwardTrackingAlg.cpp:52'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(ForwardTrackingAlg, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'SiliconTracking'
  def getType( self ):
      return 'ForwardTrackingAlg'
  pass # class ForwardTrackingAlg

class SiliconTrackingAlg( ConfigurableAlgorithm ) :
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
    'LayerCombinations' : [ 8 , 6 , 5 , 8 , 6 , 4 , 8 , 6 , 3 , 8 , 6 , 2 , 8 , 5 , 3 , 8 , 5 , 2 , 8 , 4 , 3 , 8 , 4 , 2 , 6 , 5 , 3 , 6 , 5 , 2 , 6 , 4 , 3 , 6 , 4 , 2 , 6 , 3 , 1 , 6 , 3 , 0 , 6 , 2 , 1 , 6 , 2 , 0 , 5 , 3 , 1 , 5 , 3 , 0 , 5 , 2 , 1 , 5 , 2 , 0 , 4 , 3 , 1 , 4 , 3 , 0 , 4 , 2 , 1 , 4 , 2 , 0 ], # list
    'LayerCombinationsFTD' : [ 4 , 3 , 2 , 4 , 3 , 1 , 4 , 3 , 0 , 4 , 2 , 1 , 4 , 2 , 0 , 4 , 1 , 0 , 3 , 2 , 1 , 3 , 2 , 0 , 3 , 1 , 0 , 2 , 1 , 0 , 9 , 8 , 7 , 9 , 8 , 6 , 9 , 8 , 5 , 9 , 7 , 6 , 9 , 7 , 5 , 9 , 6 , 5 , 8 , 7 , 6 , 8 , 7 , 5 , 8 , 6 , 5 , 7 , 6 , 5 ], # list
    'NDivisionsInPhi' : 80, # int
    'NDivisionsInPhiFTD' : 30, # int
    'NDivisionsInTheta' : 80, # int
    'Chi2WRphiTriplet' : 1.00000, # float
    'Chi2WRphiQuartet' : 1.00000, # float
    'Chi2WRphiSeptet' : 1.00000, # float
    'Chi2WZTriplet' : 0.500000, # float
    'Chi2WZQuartet' : 0.500000, # float
    'Chi2WZSeptet' : 0.500000, # float
    'Chi2FitCut' : 120.000, # float
    'AngleCutForMerging' : 0.100000, # float
    'MinDistCutAttach' : 2.50000, # float
    'MinLayerToAttach' : -1.00000, # float
    'CutOnD0' : 100.000, # float
    'CutOnZ0' : 100.000, # float
    'CutOnPt' : 0.0500000, # float
    'MinimalHits' : 3, # int
    'NHitsChi2' : 5, # int
    'MaxHitsPerSector' : 100, # int
    'FastAttachment' : 0, # int
    'UseSIT' : True, # bool
    'InitialTrackErrorD0' : 1.00000e+06, # float
    'InitialTrackErrorPhi0' : 100.000, # float
    'InitialTrackErrorOmega' : 0.000100000, # float
    'InitialTrackErrorZ0' : 1.00000e+06, # float
    'InitialTrackErrorTanL' : 100.000, # float
    'MaxChi2PerHit' : 100.000, # float
    'CheckForDelta' : 1, # int
    'MinDistToDelta' : 0.250000, # float
    'MultipleScatteringOn' : True, # bool
    'EnergyLossOn' : True, # bool
    'SmoothOn' : True, # bool
    'HelixMaxR' : 2000.00, # float
    'DumpTime' : True, # bool
    'HeaderCol' : DataHandle('EventHeaderCol', 'R', 'DataWrapper<edm4hep::EventHeaderCollection>'), # DataHandle
    'MCParticleCollection' : DataHandle('MCParticle', 'R', 'DataWrapper<edm4hep::MCParticleCollection>'), # DataHandle
    'VTXHitCollection' : DataHandle('VXDTrackerHits', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'FTDPixelHitCollection' : DataHandle('FTDPixelTrackerHits', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'FTDSpacePointCollection' : DataHandle('FTDSpacePoints', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'SITHitCollection' : DataHandle('SITSpacePoints', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'SITRawHitCollection' : DataHandle('SITTrackerHits', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'FTDRawHitCollection' : DataHandle('FTDStripTrackerHits', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'SiTrackCollection' : DataHandle('SiTracks', 'W', 'DataWrapper<edm4hep::TrackCollection>'), # DataHandle
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
    'LayerCombinations' : """  [SiliconTrackingAlg] """,
    'LayerCombinationsFTD' : """  [SiliconTrackingAlg] """,
    'NDivisionsInPhi' : """  [SiliconTrackingAlg] """,
    'NDivisionsInPhiFTD' : """  [SiliconTrackingAlg] """,
    'NDivisionsInTheta' : """  [SiliconTrackingAlg] """,
    'Chi2WRphiTriplet' : """  [SiliconTrackingAlg] """,
    'Chi2WRphiQuartet' : """  [SiliconTrackingAlg] """,
    'Chi2WRphiSeptet' : """  [SiliconTrackingAlg] """,
    'Chi2WZTriplet' : """  [SiliconTrackingAlg] """,
    'Chi2WZQuartet' : """  [SiliconTrackingAlg] """,
    'Chi2WZSeptet' : """  [SiliconTrackingAlg] """,
    'Chi2FitCut' : """  [SiliconTrackingAlg] """,
    'AngleCutForMerging' : """  [SiliconTrackingAlg] """,
    'MinDistCutAttach' : """  [SiliconTrackingAlg] """,
    'MinLayerToAttach' : """  [SiliconTrackingAlg] """,
    'CutOnD0' : """  [SiliconTrackingAlg] """,
    'CutOnZ0' : """  [SiliconTrackingAlg] """,
    'CutOnPt' : """  [SiliconTrackingAlg] """,
    'MinimalHits' : """  [SiliconTrackingAlg] """,
    'NHitsChi2' : """  [SiliconTrackingAlg] """,
    'MaxHitsPerSector' : """  [SiliconTrackingAlg] """,
    'FastAttachment' : """  [SiliconTrackingAlg] """,
    'UseSIT' : """  [SiliconTrackingAlg] """,
    'InitialTrackErrorD0' : """  [SiliconTrackingAlg] """,
    'InitialTrackErrorPhi0' : """  [SiliconTrackingAlg] """,
    'InitialTrackErrorOmega' : """  [SiliconTrackingAlg] """,
    'InitialTrackErrorZ0' : """  [SiliconTrackingAlg] """,
    'InitialTrackErrorTanL' : """  [SiliconTrackingAlg] """,
    'MaxChi2PerHit' : """  [SiliconTrackingAlg] """,
    'CheckForDelta' : """  [SiliconTrackingAlg] """,
    'MinDistToDelta' : """  [SiliconTrackingAlg] """,
    'MultipleScatteringOn' : """  [SiliconTrackingAlg] """,
    'EnergyLossOn' : """  [SiliconTrackingAlg] """,
    'SmoothOn' : """  [SiliconTrackingAlg] """,
    'HelixMaxR' : """  [SiliconTrackingAlg] """,
    'DumpTime' : """  [SiliconTrackingAlg] """,
    'MCParticleCollection' : """ Handle of the Input MCParticle collection [unknown owner type] """,
    'VTXHitCollection' : """ Handle of the Input VTX TrackerHits collection [unknown owner type] """,
    'FTDPixelHitCollection' : """ Handle of the Input FTD TrackerHits collection [unknown owner type] """,
    'FTDSpacePointCollection' : """ Handle of the Input FTD SpacePoints collection [unknown owner type] """,
    'SITHitCollection' : """ Handle of the Input SIT TrackerHits collection [unknown owner type] """,
    'SITRawHitCollection' : """ Handle of the SIT SpacePoints raw hit collection collection [unknown owner type] """,
    'FTDRawHitCollection' : """ Handle of the FTD SpacePoints raw hit collection collection [unknown owner type] """,
    'SiTrackCollection' : """ Handle of the SiTrack output collection [unknown owner type] """,
  }
  __declaration_location__ = 'SiliconTrackingAlg.cpp:73'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(SiliconTrackingAlg, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'SiliconTracking'
  def getType( self ):
      return 'SiliconTrackingAlg'
  pass # class SiliconTrackingAlg

class SpacePointBuilderAlg( ConfigurableAlgorithm ) :
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
    'NominalVertexX' : 0.00000, # float
    'NominalVertexY' : 0.00000, # float
    'NominalVertexZ' : 0.00000, # float
    'StriplengthTolerance' : 0.100000, # float
    'MCParticleCollection' : DataHandle('MCParticle', 'R', 'DataWrapper<edm4hep::MCParticleCollection>'), # DataHandle
    'TrackerHitCollection' : DataHandle('FTDStripTrackerHits', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'TrackerHitAssociationCollection' : DataHandle('FTDStripTrackerHitsAssociation', 'R', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
    'SpacePointCollection' : DataHandle('FTDSpacePoints', 'W', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'SpacePointAssociationCollection' : DataHandle('FTDSpacePointsAssociation', 'W', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
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
    'NominalVertexX' : """  [SpacePointBuilderAlg] """,
    'NominalVertexY' : """  [SpacePointBuilderAlg] """,
    'NominalVertexZ' : """  [SpacePointBuilderAlg] """,
    'StriplengthTolerance' : """  [SpacePointBuilderAlg] """,
    'MCParticleCollection' : """ Handle of the Input MCParticle collection [unknown owner type] """,
    'TrackerHitCollection' : """ Handle of the Input TrackerHits collection [unknown owner type] """,
    'TrackerHitAssociationCollection' : """ Handle of the Input MCRecoTrackerAssociation collection [unknown owner type] """,
    'SpacePointCollection' : """ Handle of the SpacePoint output collection [unknown owner type] """,
    'SpacePointAssociationCollection' : """ Handle of the SpacePoints association output collection [unknown owner type] """,
  }
  __declaration_location__ = 'SpacePointBuilderAlg.cpp:26'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(SpacePointBuilderAlg, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'SiliconTracking'
  def getType( self ):
      return 'SpacePointBuilderAlg'
  pass # class SpacePointBuilderAlg

class TrackSubsetAlg( ConfigurableAlgorithm ) :
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
    'TrackInputCollections' : [ 'ForwardTracks' , 'SiTracks' ], # list
    'RawTrackerHitCollections' : [ 'VXDTrackerHits' , 'SITTrackerHits' , 'FTDPixelTrackerHits' , 'FTDStripTrackerHits' ], # list
    'MultipleScatteringOn' : True, # bool
    'EnergyLossOn' : True, # bool
    'SmoothOn' : False, # bool
    'InitialTrackErrorD0' : 1.00000e+06, # float
    'InitialTrackErrorPhi0' : 100.000, # float
    'InitialTrackErrorOmega' : 0.000100000, # float
    'InitialTrackErrorZ0' : 1.00000e+06, # float
    'InitialTrackErrorTanL' : 100.000, # float
    'MaxChi2PerHit' : 100.00000, # float
    'Omega' : 0.75000000, # float
    'DumpTime' : True, # bool
    'TrackSubsetCollection' : DataHandle('SubsetTracks', 'W', 'DataWrapper<edm4hep::TrackCollection>'), # DataHandle
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
    'TrackInputCollections' : """  [TrackSubsetAlg] """,
    'RawTrackerHitCollections' : """  [TrackSubsetAlg] """,
    'MultipleScatteringOn' : """  [TrackSubsetAlg] """,
    'EnergyLossOn' : """  [TrackSubsetAlg] """,
    'SmoothOn' : """  [TrackSubsetAlg] """,
    'InitialTrackErrorD0' : """  [TrackSubsetAlg] """,
    'InitialTrackErrorPhi0' : """  [TrackSubsetAlg] """,
    'InitialTrackErrorOmega' : """  [TrackSubsetAlg] """,
    'InitialTrackErrorZ0' : """  [TrackSubsetAlg] """,
    'InitialTrackErrorTanL' : """  [TrackSubsetAlg] """,
    'MaxChi2PerHit' : """  [TrackSubsetAlg] """,
    'Omega' : """  [TrackSubsetAlg] """,
    'DumpTime' : """  [TrackSubsetAlg] """,
    'TrackSubsetCollection' : """ Handle of the SiTrack output collection [unknown owner type] """,
  }
  __declaration_location__ = 'TrackSubsetAlg.cpp:22'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(TrackSubsetAlg, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'SiliconTracking'
  def getType( self ):
      return 'TrackSubsetAlg'
  pass # class TrackSubsetAlg
