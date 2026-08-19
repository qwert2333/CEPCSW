#Sun Mar 10 20:54:22 2024"""Automatically generated. DO NOT EDIT please"""
import sys
if sys.version_info >= (3,):
    # Python 2 compatibility
    long = int
from GaudiKernel.Proxy.Configurable import *

class BetheBlochEquationDedxSimTool( ConfigurableAlgTool ) :
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
    'material_Z' : 2.00000, # float
    'material_A' : 4.00000, # float
    'material_density' : 0.000178000, # float
    'scale' : 1.00000, # float
    'resolution' : 0.00000, # float
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
    'material_Z' : """  [BetheBlochEquationDedxSimTool] """,
    'material_A' : """  [BetheBlochEquationDedxSimTool] """,
    'material_density' : """  [BetheBlochEquationDedxSimTool] """,
    'scale' : """  [BetheBlochEquationDedxSimTool] """,
    'resolution' : """  [BetheBlochEquationDedxSimTool] """,
  }
  __declaration_location__ = 'BetheBlochEquationDedxSimTool.cpp:7'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(BetheBlochEquationDedxSimTool, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'DetSimDedx'
  def getType( self ):
      return 'BetheBlochEquationDedxSimTool'
  pass # class BetheBlochEquationDedxSimTool

class DummyDedxSimTool( ConfigurableAlgTool ) :
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
  __declaration_location__ = 'DummyDedxSimTool.cpp:5'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(DummyDedxSimTool, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'DetSimDedx'
  def getType( self ):
      return 'DummyDedxSimTool'
  pass # class DummyDedxSimTool

class GFDndxSimTool( ConfigurableAlgTool ) :
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
  __declaration_location__ = 'GFDndxSimTool.cpp:5'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(GFDndxSimTool, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'DetSimDedx'
  def getType( self ):
      return 'GFDndxSimTool'
  pass # class GFDndxSimTool

class TrackHeedSimTool( ConfigurableAlgTool ) :
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
    'readout' : 'DriftChamberHitsCollection', # str
    'gas_file' : 'He_50_isobutane_50.gas', # str
    'IonMobility_file' : 'IonMobility_He+_He.txt', # str
    'isob' : 50.0000, # float
    'he' : 50.0000, # float
    'debug' : False, # bool
    'use_max_step' : False, # bool
    'update_KE' : True, # bool
    'max_step' : 1.00000, # float
    'only_primary' : False, # bool
    'save_mc' : False, # bool
    'save_cellID' : True, # bool
    'delta_threshold' : 50.0000, # float
    'change_threshold' : 0.0500000, # float
    'BField' : -3.00000, # float
    'eps' : 1.00000e-06, # float
    'intraOpNumThreads' : 1, # int
    'interOpNumThreads' : 1, # int
    'sim_pulse' : True, # bool
    'model' : 'model_test.onnx', # str
    'batchsize' : 100, # int
    'time_scale' : 503.000, # float
    'time_shift' : 814.000, # float
    'amp_scale' : 1.15000, # float
    'amp_shift' : 0.860000, # float
    'x_scale' : 9.00000, # float
    'y_scale' : 9.00000, # float
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
    'readout' : """  [TrackHeedSimTool] """,
    'gas_file' : """  [TrackHeedSimTool] """,
    'IonMobility_file' : """  [TrackHeedSimTool] """,
    'isob' : """  [TrackHeedSimTool] """,
    'he' : """  [TrackHeedSimTool] """,
    'debug' : """  [TrackHeedSimTool] """,
    'use_max_step' : """  [TrackHeedSimTool] """,
    'update_KE' : """  [TrackHeedSimTool] """,
    'max_step' : """  [TrackHeedSimTool] """,
    'only_primary' : """  [TrackHeedSimTool] """,
    'save_mc' : """  [TrackHeedSimTool] """,
    'save_cellID' : """  [TrackHeedSimTool] """,
    'delta_threshold' : """  [TrackHeedSimTool] """,
    'change_threshold' : """  [TrackHeedSimTool] """,
    'BField' : """  [TrackHeedSimTool] """,
    'eps' : """  [TrackHeedSimTool] """,
    'intraOpNumThreads' : """  [TrackHeedSimTool] """,
    'interOpNumThreads' : """  [TrackHeedSimTool] """,
    'sim_pulse' : """  [TrackHeedSimTool] """,
    'model' : """  [TrackHeedSimTool] """,
    'batchsize' : """  [TrackHeedSimTool] """,
    'time_scale' : """  [TrackHeedSimTool] """,
    'time_shift' : """  [TrackHeedSimTool] """,
    'amp_scale' : """  [TrackHeedSimTool] """,
    'amp_shift' : """  [TrackHeedSimTool] """,
    'x_scale' : """  [TrackHeedSimTool] """,
    'y_scale' : """  [TrackHeedSimTool] """,
  }
  __declaration_location__ = 'TrackHeedSimTool.cpp:21'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(TrackHeedSimTool, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'DetSimDedx'
  def getType( self ):
      return 'TrackHeedSimTool'
  pass # class TrackHeedSimTool
