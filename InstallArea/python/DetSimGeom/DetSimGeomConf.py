#Sun Mar 10 20:54:17 2024"""Automatically generated. DO NOT EDIT please"""
import sys
if sys.version_info >= (3,):
    # Python 2 compatibility
    long = int
from GaudiKernel.Proxy.Configurable import *

class AnExampleDetElemTool( ConfigurableAlgTool ) :
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
    'X' : 30000.000, # float
    'Y' : 30000.000, # float
    'Z' : 30000.000, # float
    'detxml' : '', # str
    'SDenabled' : True, # bool
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
    'X' : """  [AnExampleDetElemTool] """,
    'Y' : """  [AnExampleDetElemTool] """,
    'Z' : """  [AnExampleDetElemTool] """,
    'detxml' : """  [AnExampleDetElemTool] """,
    'SDenabled' : """  [AnExampleDetElemTool] """,
  }
  __declaration_location__ = 'AnExampleDetElemTool.cpp:33'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(AnExampleDetElemTool, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'DetSimGeom'
  def getType( self ):
      return 'AnExampleDetElemTool'
  pass # class AnExampleDetElemTool

class WorldDetElemTool( ConfigurableAlgTool ) :
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
  __declaration_location__ = 'WorldDetElemTool.cpp:25'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(WorldDetElemTool, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'DetSimGeom'
  def getType( self ):
      return 'WorldDetElemTool'
  pass # class WorldDetElemTool
