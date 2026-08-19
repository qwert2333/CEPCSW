#Sun Mar 10 20:54:38 2024"""Automatically generated. DO NOT EDIT please"""
import sys
if sys.version_info >= (3,):
    # Python 2 compatibility
    long = int
from GaudiKernel.Proxy.Configurable import *

class CalorimeterSensDetTool( ConfigurableAlgTool ) :
  __slots__ = { 
    'ExtraInputs' : [], # list
    'ExtraOutputs' : [], # list
    'OutputLevel' : 0, # int
    'MonitorService' : 'MonitorSvc', # str
    'AuditTools' : False, # bool
    'AuditInitialize' : False, # bool
    'AuditStart' : False, # bool
    'AuditStop' : False, # bool
    'AuditFinalize' : False, # bool
    'AuditReinitialize' : False, # bool
    'AuditRestart' : False, # bool
    'CalNamesMergeDisable' : [  ], # list
    'CalNamesApplyBirks' : [  ], # list
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgTool,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgTool,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [AlgTool] """,
    'MonitorService' : """ name to use for Monitor Service [AlgTool] """,
    'AuditTools' : """ [[deprecated]] unused [AlgTool] """,
    'AuditInitialize' : """ trigger auditor on initialize() [AlgTool] """,
    'AuditStart' : """ trigger auditor on start() [AlgTool] """,
    'AuditStop' : """ trigger auditor on stop() [AlgTool] """,
    'AuditFinalize' : """ trigger auditor on finalize() [AlgTool] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [AlgTool] """,
    'AuditRestart' : """ trigger auditor on restart() [AlgTool] """,
    'CalNamesMergeDisable' : """  [CalorimeterSensDetTool] """,
    'CalNamesApplyBirks' : """  [CalorimeterSensDetTool] """,
  }
  __declaration_location__ = 'CalorimeterSensDetTool.cpp:9'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(CalorimeterSensDetTool, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'DetSimSD'
  def getType( self ):
      return 'CalorimeterSensDetTool'
  pass # class CalorimeterSensDetTool

class DriftChamberSensDetTool( ConfigurableAlgTool ) :
  __slots__ = { 
    'ExtraInputs' : [], # list
    'ExtraOutputs' : [], # list
    'OutputLevel' : 0, # int
    'MonitorService' : 'MonitorSvc', # str
    'AuditTools' : False, # bool
    'AuditInitialize' : False, # bool
    'AuditStart' : False, # bool
    'AuditStop' : False, # bool
    'AuditFinalize' : False, # bool
    'AuditReinitialize' : False, # bool
    'AuditRestart' : False, # bool
    'DedxSimTool' : '', # str
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgTool,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgTool,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [AlgTool] """,
    'MonitorService' : """ name to use for Monitor Service [AlgTool] """,
    'AuditTools' : """ [[deprecated]] unused [AlgTool] """,
    'AuditInitialize' : """ trigger auditor on initialize() [AlgTool] """,
    'AuditStart' : """ trigger auditor on start() [AlgTool] """,
    'AuditStop' : """ trigger auditor on stop() [AlgTool] """,
    'AuditFinalize' : """ trigger auditor on finalize() [AlgTool] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [AlgTool] """,
    'AuditRestart' : """ trigger auditor on restart() [AlgTool] """,
    'DedxSimTool' : """  [DriftChamberSensDetTool] """,
  }
  __declaration_location__ = 'DriftChamberSensDetTool.cpp:9'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(DriftChamberSensDetTool, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'DetSimSD'
  def getType( self ):
      return 'DriftChamberSensDetTool'
  pass # class DriftChamberSensDetTool

class GenericTrackerSensDetTool( ConfigurableAlgTool ) :
  __slots__ = { 
    'ExtraInputs' : [], # list
    'ExtraOutputs' : [], # list
    'OutputLevel' : 0, # int
    'MonitorService' : 'MonitorSvc', # str
    'AuditTools' : False, # bool
    'AuditInitialize' : False, # bool
    'AuditStart' : False, # bool
    'AuditStop' : False, # bool
    'AuditFinalize' : False, # bool
    'AuditReinitialize' : False, # bool
    'AuditRestart' : False, # bool
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgTool,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgTool,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [AlgTool] """,
    'MonitorService' : """ name to use for Monitor Service [AlgTool] """,
    'AuditTools' : """ [[deprecated]] unused [AlgTool] """,
    'AuditInitialize' : """ trigger auditor on initialize() [AlgTool] """,
    'AuditStart' : """ trigger auditor on start() [AlgTool] """,
    'AuditStop' : """ trigger auditor on stop() [AlgTool] """,
    'AuditFinalize' : """ trigger auditor on finalize() [AlgTool] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [AlgTool] """,
    'AuditRestart' : """ trigger auditor on restart() [AlgTool] """,
  }
  __declaration_location__ = 'GenericTrackerSensDetTool.cpp:12'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(GenericTrackerSensDetTool, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'DetSimSD'
  def getType( self ):
      return 'GenericTrackerSensDetTool'
  pass # class GenericTrackerSensDetTool

class TimeProjectionChamberSensDetTool( ConfigurableAlgTool ) :
  __slots__ = { 
    'ExtraInputs' : [], # list
    'ExtraOutputs' : [], # list
    'OutputLevel' : 0, # int
    'MonitorService' : 'MonitorSvc', # str
    'AuditTools' : False, # bool
    'AuditInitialize' : False, # bool
    'AuditStart' : False, # bool
    'AuditStop' : False, # bool
    'AuditFinalize' : False, # bool
    'AuditReinitialize' : False, # bool
    'AuditRestart' : False, # bool
    'TypeOption' : 0, # int
    'ThresholdEnergyDeposit' : 3.2000000e-08, # float
    'LowPtCut' : 0.010000000, # float
    'LowPtMaxHitSeparation' : 0.50000000, # float
    'SameStepLimit' : True, # bool
    'WriteMCTruthForLowPtHits' : False, # bool
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgTool,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgTool,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [AlgTool] """,
    'MonitorService' : """ name to use for Monitor Service [AlgTool] """,
    'AuditTools' : """ [[deprecated]] unused [AlgTool] """,
    'AuditInitialize' : """ trigger auditor on initialize() [AlgTool] """,
    'AuditStart' : """ trigger auditor on start() [AlgTool] """,
    'AuditStop' : """ trigger auditor on stop() [AlgTool] """,
    'AuditFinalize' : """ trigger auditor on finalize() [AlgTool] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [AlgTool] """,
    'AuditRestart' : """ trigger auditor on restart() [AlgTool] """,
    'TypeOption' : """  [TimeProjectionChamberSensDetTool] """,
    'ThresholdEnergyDeposit' : """  [TimeProjectionChamberSensDetTool] """,
    'LowPtCut' : """  [TimeProjectionChamberSensDetTool] """,
    'LowPtMaxHitSeparation' : """  [TimeProjectionChamberSensDetTool] """,
    'SameStepLimit' : """  [TimeProjectionChamberSensDetTool] """,
    'WriteMCTruthForLowPtHits' : """  [TimeProjectionChamberSensDetTool] """,
  }
  __declaration_location__ = 'TimeProjectionChamberSensDetTool.cpp:11'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(TimeProjectionChamberSensDetTool, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'DetSimSD'
  def getType( self ):
      return 'TimeProjectionChamberSensDetTool'
  pass # class TimeProjectionChamberSensDetTool
