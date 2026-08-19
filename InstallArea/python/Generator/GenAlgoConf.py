#Sun Mar 10 20:53:58 2024"""Automatically generated. DO NOT EDIT please"""
import sys
if sys.version_info >= (3,):
    # Python 2 compatibility
    long = int
from GaudiKernel.DataHandle import DataHandle
from GaudiKernel.Proxy.Configurable import *

class GenAlgo( ConfigurableAlgorithm ) :
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
    'Input' : 'NULL', # str
    'FileFormat' : 'NULL', # str
    'OutputRootFile' : 'NULL', # str
    'PrintEvent' : True, # bool
    'WriteFile' : True, # bool
    'MCParticleGen' : DataHandle('MCParticleGen', 'W', 'DataWrapper<edm4hep::MCParticleCollection>'), # DataHandle
    'GenTools' : [  ], # list
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
    'Input' : """  [GenAlgo] """,
    'FileFormat' : """  [GenAlgo] """,
    'OutputRootFile' : """  [GenAlgo] """,
    'PrintEvent' : """  [GenAlgo] """,
    'WriteFile' : """  [GenAlgo] """,
    'MCParticleGen' : """ MCParticle collection (at Generator phase) [unknown owner type] """,
    'GenTools' : """ List of GenTools [unknown owner type] """,
  }
  __declaration_location__ = 'GenAlgo.cpp:18'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(GenAlgo, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'GenAlgo'
  def getType( self ):
      return 'GenAlgo'
  pass # class GenAlgo

class GenPrinter( ConfigurableAlgTool ) :
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
  __declaration_location__ = 'GenPrinter.cpp:4'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(GenPrinter, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'GenAlgo'
  def getType( self ):
      return 'GenPrinter'
  pass # class GenPrinter

class GtBeamBackgroundTool( ConfigurableAlgTool ) :
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
    'InputFileMap' : {  }, # list
    'InputFormatMap' : {  }, # list
    'InputRateMap' : {  }, # list
    'InputBeamEnergyMap' : {  }, # list
    'RotationAlongYMap' : {  }, # list
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
    'InputFileMap' : """  [GtBeamBackgroundTool] """,
    'InputFormatMap' : """  [GtBeamBackgroundTool] """,
    'InputRateMap' : """  [GtBeamBackgroundTool] """,
    'InputBeamEnergyMap' : """  [GtBeamBackgroundTool] """,
    'RotationAlongYMap' : """  [GtBeamBackgroundTool] """,
  }
  __declaration_location__ = 'GtBeamBackgroundTool.cpp:8'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(GtBeamBackgroundTool, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'GenAlgo'
  def getType( self ):
      return 'GtBeamBackgroundTool'
  pass # class GtBeamBackgroundTool

class GtGunTool( ConfigurableAlgTool ) :
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
    'Particles' : [  ], # list
    'PositionXs' : [  ], # list
    'PositionYs' : [  ], # list
    'PositionZs' : [  ], # list
    'PosZMins' : [  ], # list
    'PosZMaxs' : [  ], # list
    'PosRMins' : [  ], # list
    'PosRMaxs' : [  ], # list
    'EnergyMins' : [  ], # list
    'EnergyMaxs' : [  ], # list
    'ThetaMins' : [  ], # list
    'ThetaMaxs' : [  ], # list
    'PhiMins' : [  ], # list
    'PhiMaxs' : [  ], # list
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
    'Particles' : """  [GtGunTool] """,
    'PositionXs' : """  [GtGunTool] """,
    'PositionYs' : """  [GtGunTool] """,
    'PositionZs' : """  [GtGunTool] """,
    'PosZMins' : """  [GtGunTool] """,
    'PosZMaxs' : """  [GtGunTool] """,
    'PosRMins' : """  [GtGunTool] """,
    'PosRMaxs' : """  [GtGunTool] """,
    'EnergyMins' : """  [GtGunTool] """,
    'EnergyMaxs' : """  [GtGunTool] """,
    'ThetaMins' : """  [GtGunTool] """,
    'ThetaMaxs' : """  [GtGunTool] """,
    'PhiMins' : """  [GtGunTool] """,
    'PhiMaxs' : """  [GtGunTool] """,
  }
  __declaration_location__ = 'GtGunTool.cpp:8'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(GtGunTool, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'GenAlgo'
  def getType( self ):
      return 'GtGunTool'
  pass # class GtGunTool

class HepMCRdr( ConfigurableAlgTool ) :
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
    'Input' : '', # str
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
    'Input' : """  [HepMCRdr] """,
  }
  __declaration_location__ = 'HepMCRdr.cpp:25'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(HepMCRdr, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'GenAlgo'
  def getType( self ):
      return 'HepMCRdr'
  pass # class HepMCRdr

class SLCIORdr( ConfigurableAlgTool ) :
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
    'Input' : '', # str
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
    'Input' : """  [SLCIORdr] """,
  }
  __declaration_location__ = 'SLCIORdr.cpp:32'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(SLCIORdr, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'GenAlgo'
  def getType( self ):
      return 'SLCIORdr'
  pass # class SLCIORdr

class StdHepRdr( ConfigurableAlgTool ) :
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
    'Input' : '', # str
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
    'Input' : """  [StdHepRdr] """,
  }
  __declaration_location__ = 'StdHepRdr.cpp:27'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(StdHepRdr, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'GenAlgo'
  def getType( self ):
      return 'StdHepRdr'
  pass # class StdHepRdr
