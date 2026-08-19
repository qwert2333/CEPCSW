#Sun Mar 10 20:55:24 2024"""Automatically generated. DO NOT EDIT please"""
import sys
if sys.version_info >= (3,):
    # Python 2 compatibility
    long = int
from GaudiKernel.DataHandle import DataHandle
from GaudiKernel.Proxy.Configurable import *

class ClupatraAlg( ConfigurableAlgorithm ) :
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
    'DistanceCut' : 40.0000, # float
    'CosAlphaCut' : 1.00000, # float
    'NLoopForSeeding' : 4, # int
    'MinimumClusterSize' : 6, # int
    'DuplicatePadRowFraction' : 0.100000, # float
    'MaxDeltaChi2' : 35.0000, # float
    'Chi2Cut' : 100.000, # float
    'MaxStepWithoutHit' : 6, # int
    'pickUpSiHits' : False, # bool
    'CreateDebugCollections' : False, # bool
    'PadRowRange' : 15, # int
    'NumberOfZBins' : 150, # int
    'MinLayerFractionWithMultiplicity' : 0.500000, # float
    'MinLayerNumberWithMultiplicity' : 3.00000, # float
    'TrackStartsInnerDist' : 25.0000, # float
    'TrackEndsOuterCentralDist' : 25.0000, # float
    'TrackEndsOuterForwardDist' : 40.0000, # float
    'TrackIsCurlerOmega' : 0.00100000, # float
    'MultipleScatteringOn' : False, # bool
    'EnergyLossOn' : True, # bool
    'SmoothOn' : False, # bool
    'TPCHitCollection' : DataHandle('TPCTrackerHits', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'ClupatraTracks' : DataHandle('ClupatraTracks', 'W', 'DataWrapper<edm4hep::TrackCollection>'), # DataHandle
    'ClupatraTrackSegments' : DataHandle('ClupatraSegmentTracks', 'W', 'DataWrapper<edm4hep::TrackCollection>'), # DataHandle
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
    'DistanceCut' : """  [ClupatraAlg] """,
    'CosAlphaCut' : """  [ClupatraAlg] """,
    'NLoopForSeeding' : """  [ClupatraAlg] """,
    'MinimumClusterSize' : """  [ClupatraAlg] """,
    'DuplicatePadRowFraction' : """  [ClupatraAlg] """,
    'MaxDeltaChi2' : """  [ClupatraAlg] """,
    'Chi2Cut' : """  [ClupatraAlg] """,
    'MaxStepWithoutHit' : """  [ClupatraAlg] """,
    'pickUpSiHits' : """  [ClupatraAlg] """,
    'CreateDebugCollections' : """  [ClupatraAlg] """,
    'PadRowRange' : """  [ClupatraAlg] """,
    'NumberOfZBins' : """  [ClupatraAlg] """,
    'MinLayerFractionWithMultiplicity' : """  [ClupatraAlg] """,
    'MinLayerNumberWithMultiplicity' : """  [ClupatraAlg] """,
    'TrackStartsInnerDist' : """  [ClupatraAlg] """,
    'TrackEndsOuterCentralDist' : """  [ClupatraAlg] """,
    'TrackEndsOuterForwardDist' : """  [ClupatraAlg] """,
    'TrackIsCurlerOmega' : """  [ClupatraAlg] """,
    'MultipleScatteringOn' : """  [ClupatraAlg] """,
    'EnergyLossOn' : """  [ClupatraAlg] """,
    'SmoothOn' : """  [ClupatraAlg] """,
    'TPCHitCollection' : """ handler of the tpc hit input collections [unknown owner type] """,
    'ClupatraTracks' : """ handler of the collection with final TPC tracks [unknown owner type] """,
    'ClupatraTrackSegments' : """ handler of the collection that has the individual track segments [unknown owner type] """,
  }
  __declaration_location__ = 'ClupatraAlg.cpp:156'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(ClupatraAlg, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'Tracking'
  def getType( self ):
      return 'ClupatraAlg'
  pass # class ClupatraAlg

class FullLDCTrackingAlg( ConfigurableAlgorithm ) :
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
    'MultipleScatteringOn' : True, # bool
    'EnergyLossOn' : True, # bool
    'SmoothOn' : True, # bool
    'Debug' : 0, # int
    'ForceSiTPCMerging' : 0, # int
    'ForceTPCSegmentsMerging' : 0, # int
    'CutOnTPCHits' : 35, # int
    'CutOnSiHits' : 4, # int
    'AssignVTXHits' : 1, # int
    'AssignFTDHits' : 1, # int
    'AssignSITHits' : 1, # int
    'AssignTPCHits' : 1, # int
    'AssignSETHits' : 1, # int
    'AssignETDHits' : 1, # int
    'ForbidOverlapInZTPC' : 0, # int
    'ForbidOverlapInZComb' : 0, # int
    'D0CutForMerging' : 500.000, # float
    'Z0CutForMerging' : 1000.00, # float
    'OmegaCutForMerging' : 0.250000, # float
    'AngleCutForMerging' : 0.100000, # float
    'Chi2FitCut' : 100.000, # float
    'D0CutForForcedMerging' : 50.0000, # float
    'Z0CutForForcedMerging' : 200.000, # float
    'OmegaCutForForcedMerging' : 0.150000, # float
    'AngleCutForForcedMerging' : 0.0500000, # float
    'D0CutToMergeTPCSegments' : 100.000, # float
    'Z0CutToMergeTPCSegments' : 5000.00, # float
    'DeltaPCutToMergeTPCSegments' : 0.100000, # float
    'PtCutToMergeTPCSegments' : 1.20000, # float
    'cosThetaCutHighPtMerge' : 0.990000, # float
    'cosThetaCutSoftHighPtMerge' : 0.998000, # float
    'momDiffCutHighPtMerge' : 0.0100000, # float
    'momDiffCutSoftHighPtMerge' : 0.250000, # float
    'hitDistanceCutHighPtMerge' : 25.0000, # float
    'maxHitDistanceCutHighPtMerge' : 50.0000, # float
    'maxFractionOfOutliersCutHighPtMerge' : 0.950000, # float
    'NHitsExtrapolation' : 35, # int
    'VTXHitToTrackDistance' : 1.50000, # float
    'FTDHitToTrackDistance' : 2.00000, # float
    'SITHitToTrackDistance' : 2.00000, # float
    'SETHitToTrackDistance' : 2.00000, # float
    'ETDHitToTrackDistance' : 10.0000, # float
    'TPCHitToTrackDistance' : 15.0000, # float
    'CutOnTrackD0' : 500.000, # float
    'CutOnTrackZ0' : 500.000, # float
    'InitialTrackErrorD0' : 1.00000e+06, # float
    'InitialTrackErrorPhi0' : 100.000, # float
    'InitialTrackErrorOmega' : 0.000100000, # float
    'InitialTrackErrorZ0' : 1.00000e+06, # float
    'InitialTrackErrorTanL' : 100.000, # float
    'MaxChi2PerHit' : 100.000, # float
    'MinChi2ProbForSiliconTracks' : 0.00100000, # float
    'VetoMergeMomentumCut' : 2.50000, # float
    'MaxAllowedPercentageOfOutliersForTrackCombination' : 0.300000, # float
    'MaxAllowedSiHitRejectionsForTrackCombination' : 2, # int
    'DumpTime' : True, # bool
    'FTDPixelTrackerHits' : DataHandle('FTDPixelTrackerHits', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'FTDSpacePoints' : DataHandle('FTDSpacePoints', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'VTXTrackerHits' : DataHandle('VTXTrackerHits', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'SITTrackerHits' : DataHandle('SITSpacePoints', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'SETTrackerHits' : DataHandle('SETSpacePoints', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'TPCTrackerHits' : DataHandle('TPCTrackerHits', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'SITRawHits' : DataHandle('SITTrackerHits', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'SETRawHits' : DataHandle('SETTrackerHits', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'FTDRawHits' : DataHandle('FTDStripTrackerHits', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'TPCTracks' : DataHandle('ClupatraTracks', 'R', 'DataWrapper<edm4hep::TrackCollection>'), # DataHandle
    'SiTracks' : DataHandle('SiTracks', 'R', 'DataWrapper<edm4hep::TrackCollection>'), # DataHandle
    'VTXHitRelCol' : DataHandle('VTXHitRelCol', 'R', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
    'SITHitRelCol' : DataHandle('SITHitRelCol', 'R', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
    'SETHitRelCol' : DataHandle('SETHitRelCol', 'R', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
    'FTDHitRelCol' : DataHandle('FTDHitRelCol', 'R', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
    'TPCHitRelCol' : DataHandle('TPCHitRelCol', 'R', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
    'OutputTracks' : DataHandle('MarlinTrkTracks', 'W', 'DataWrapper<edm4hep::TrackCollection>'), # DataHandle
    'OutputTrackAssociationCollection' : DataHandle('MarlinTrkAssociation', 'W', 'DataWrapper<edm4hep::MCRecoTrackParticleAssociationCollection>'), # DataHandle
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
    'MultipleScatteringOn' : """  [FullLDCTrackingAlg] """,
    'EnergyLossOn' : """  [FullLDCTrackingAlg] """,
    'SmoothOn' : """  [FullLDCTrackingAlg] """,
    'Debug' : """  [FullLDCTrackingAlg] """,
    'ForceSiTPCMerging' : """  [FullLDCTrackingAlg] """,
    'ForceTPCSegmentsMerging' : """  [FullLDCTrackingAlg] """,
    'CutOnTPCHits' : """  [FullLDCTrackingAlg] """,
    'CutOnSiHits' : """  [FullLDCTrackingAlg] """,
    'AssignVTXHits' : """  [FullLDCTrackingAlg] """,
    'AssignFTDHits' : """  [FullLDCTrackingAlg] """,
    'AssignSITHits' : """  [FullLDCTrackingAlg] """,
    'AssignTPCHits' : """  [FullLDCTrackingAlg] """,
    'AssignSETHits' : """  [FullLDCTrackingAlg] """,
    'AssignETDHits' : """  [FullLDCTrackingAlg] """,
    'ForbidOverlapInZTPC' : """  [FullLDCTrackingAlg] """,
    'ForbidOverlapInZComb' : """  [FullLDCTrackingAlg] """,
    'D0CutForMerging' : """  [FullLDCTrackingAlg] """,
    'Z0CutForMerging' : """  [FullLDCTrackingAlg] """,
    'OmegaCutForMerging' : """  [FullLDCTrackingAlg] """,
    'AngleCutForMerging' : """  [FullLDCTrackingAlg] """,
    'Chi2FitCut' : """  [FullLDCTrackingAlg] """,
    'D0CutForForcedMerging' : """  [FullLDCTrackingAlg] """,
    'Z0CutForForcedMerging' : """  [FullLDCTrackingAlg] """,
    'OmegaCutForForcedMerging' : """  [FullLDCTrackingAlg] """,
    'AngleCutForForcedMerging' : """  [FullLDCTrackingAlg] """,
    'D0CutToMergeTPCSegments' : """  [FullLDCTrackingAlg] """,
    'Z0CutToMergeTPCSegments' : """  [FullLDCTrackingAlg] """,
    'DeltaPCutToMergeTPCSegments' : """  [FullLDCTrackingAlg] """,
    'PtCutToMergeTPCSegments' : """  [FullLDCTrackingAlg] """,
    'cosThetaCutHighPtMerge' : """  [FullLDCTrackingAlg] """,
    'cosThetaCutSoftHighPtMerge' : """  [FullLDCTrackingAlg] """,
    'momDiffCutHighPtMerge' : """  [FullLDCTrackingAlg] """,
    'momDiffCutSoftHighPtMerge' : """  [FullLDCTrackingAlg] """,
    'hitDistanceCutHighPtMerge' : """  [FullLDCTrackingAlg] """,
    'maxHitDistanceCutHighPtMerge' : """  [FullLDCTrackingAlg] """,
    'maxFractionOfOutliersCutHighPtMerge' : """  [FullLDCTrackingAlg] """,
    'NHitsExtrapolation' : """  [FullLDCTrackingAlg] """,
    'VTXHitToTrackDistance' : """  [FullLDCTrackingAlg] """,
    'FTDHitToTrackDistance' : """  [FullLDCTrackingAlg] """,
    'SITHitToTrackDistance' : """  [FullLDCTrackingAlg] """,
    'SETHitToTrackDistance' : """  [FullLDCTrackingAlg] """,
    'ETDHitToTrackDistance' : """  [FullLDCTrackingAlg] """,
    'TPCHitToTrackDistance' : """  [FullLDCTrackingAlg] """,
    'CutOnTrackD0' : """  [FullLDCTrackingAlg] """,
    'CutOnTrackZ0' : """  [FullLDCTrackingAlg] """,
    'InitialTrackErrorD0' : """  [FullLDCTrackingAlg] """,
    'InitialTrackErrorPhi0' : """  [FullLDCTrackingAlg] """,
    'InitialTrackErrorOmega' : """  [FullLDCTrackingAlg] """,
    'InitialTrackErrorZ0' : """  [FullLDCTrackingAlg] """,
    'InitialTrackErrorTanL' : """  [FullLDCTrackingAlg] """,
    'MaxChi2PerHit' : """  [FullLDCTrackingAlg] """,
    'MinChi2ProbForSiliconTracks' : """  [FullLDCTrackingAlg] """,
    'VetoMergeMomentumCut' : """  [FullLDCTrackingAlg] """,
    'MaxAllowedPercentageOfOutliersForTrackCombination' : """  [FullLDCTrackingAlg] """,
    'MaxAllowedSiHitRejectionsForTrackCombination' : """  [FullLDCTrackingAlg] """,
    'DumpTime' : """  [FullLDCTrackingAlg] """,
    'FTDPixelTrackerHits' : """ handler of FTD Pixel Hit Collection [unknown owner type] """,
    'FTDSpacePoints' : """ FTD FTDSpacePoint Collection [unknown owner type] """,
    'VTXTrackerHits' : """ VTX Hit Collection [unknown owner type] """,
    'SITTrackerHits' : """ SIT Hit Collection [unknown owner type] """,
    'SETTrackerHits' : """ SET Hit Collection [unknown owner type] """,
    'TPCTrackerHits' : """ TPC Hit Collection [unknown owner type] """,
    'SITRawHits' : """ SIT Raw Hit Collection of SpacePoints [unknown owner type] """,
    'SETRawHits' : """ SET Raw Hit Collection of SpacePoints [unknown owner type] """,
    'FTDRawHits' : """ FTD Raw Hit Collection of SpacePoints [unknown owner type] """,
    'TPCTracks' : """ TPC Track Collection [unknown owner type] """,
    'SiTracks' : """ Si Track Collection [unknown owner type] """,
    'VTXHitRelCol' : """ Vertex TrackerHit association name [unknown owner type] """,
    'SITHitRelCol' : """ SIT TrackerHit association name [unknown owner type] """,
    'SETHitRelCol' : """ SET TrackerHit association name [unknown owner type] """,
    'FTDHitRelCol' : """ FTD TrackerHit association name [unknown owner type] """,
    'TPCHitRelCol' : """ TPC TrackerHit association name [unknown owner type] """,
    'OutputTracks' : """ Full LDC track collection name [unknown owner type] """,
    'OutputTrackAssociationCollection' : """ Track - MCParticle association collection name [unknown owner type] """,
  }
  __declaration_location__ = 'FullLDCTrackingAlg.cpp:103'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(FullLDCTrackingAlg, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'Tracking'
  def getType( self ):
      return 'FullLDCTrackingAlg'
  pass # class FullLDCTrackingAlg

class TruthTrackerAlg( ConfigurableAlgorithm ) :
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
    'writeRecParticle' : False, # bool
    'resPT' : 0.00000, # float
    'resPz' : 0.00000, # float
    'resMomPhi' : 0.00000, # float
    'resMomTheta' : 0.00000, # float
    'resVertexX' : 0.00300000, # float
    'resVertexY' : 0.00300000, # float
    'resVertexZ' : 0.00300000, # float
    'maxDigiCut' : 1000000, # int
    'MCParticle' : DataHandle('MCParticle', 'R', 'DataWrapper<edm4hep::MCParticleCollection>'), # DataHandle
    'DigiDCHitCollection' : DataHandle('DigiDCHitCollection', 'R', 'DataWrapper<edm4hep::TrackerHitCollection>'), # DataHandle
    'DCHitAssociationCollection' : DataHandle('DCHitAssociationCollection', 'R', 'DataWrapper<edm4hep::MCRecoTrackerAssociationCollection>'), # DataHandle
    'DCTrackCollection' : DataHandle('DCTrackCollection', 'W', 'DataWrapper<edm4hep::TrackCollection>'), # DataHandle
    'SiSubsetTrackCollection' : DataHandle('SiSubsetTrackCollection', 'R', 'DataWrapper<edm4hep::TrackCollection>'), # DataHandle
    'SDTTrackCollection' : DataHandle('SDTTrackCollection', 'W', 'DataWrapper<edm4hep::TrackCollection>'), # DataHandle
    'DCRecParticleCollection' : DataHandle('DCRecParticleCollection', 'W', 'DataWrapper<edm4hep::ReconstructedParticleCollection>'), # DataHandle
    'DCRecParticleAssociationCollection' : DataHandle('DCRecMCRecoParticleAssociationCollection', 'W', 'DataWrapper<edm4hep::MCRecoParticleAssociationCollection>'), # DataHandle
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
    'readout' : """  [TruthTrackerAlg] """,
    'writeRecParticle' : """  [TruthTrackerAlg] """,
    'resPT' : """  [TruthTrackerAlg] """,
    'resPz' : """  [TruthTrackerAlg] """,
    'resMomPhi' : """  [TruthTrackerAlg] """,
    'resMomTheta' : """  [TruthTrackerAlg] """,
    'resVertexX' : """  [TruthTrackerAlg] """,
    'resVertexY' : """  [TruthTrackerAlg] """,
    'resVertexZ' : """  [TruthTrackerAlg] """,
    'maxDigiCut' : """  [TruthTrackerAlg] """,
    'MCParticle' : """ Handle of the input MCParticle collection [unknown owner type] """,
    'DigiDCHitCollection' : """ Handle of DC digi(TrackerHit) collection [unknown owner type] """,
    'DCHitAssociationCollection' : """ Handle of association collection [unknown owner type] """,
    'DCTrackCollection' : """ Handle of DC track collection [unknown owner type] """,
    'SiSubsetTrackCollection' : """ Handle of silicon subset track collection [unknown owner type] """,
    'SDTTrackCollection' : """ Handle of SDT track collection [unknown owner type] """,
    'DCRecParticleCollection' : """ Handle of drift chamber reconstructed particle collection [unknown owner type] """,
    'DCRecParticleAssociationCollection' : """ Handle of drift chamber reconstructed particle collection [unknown owner type] """,
  }
  __declaration_location__ = 'TruthTrackerAlg.cpp:22'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(TruthTrackerAlg, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'Tracking'
  def getType( self ):
      return 'TruthTrackerAlg'
  pass # class TruthTrackerAlg
