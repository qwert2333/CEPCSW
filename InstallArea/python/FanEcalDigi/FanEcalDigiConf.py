#Wed Apr 10 22:38:24 2024"""Automatically generated. DO NOT EDIT please"""
import sys
if sys.version_info >= (3,):
    # Python 2 compatibility
    long = int
from GaudiKernel.DataHandle import DataHandle
from GaudiKernel.Proxy.Configurable import *

class FanEcalDigiAlg( ConfigurableAlgorithm ) :
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
    'ReadOutName' : 'CaloHitsCollection', # str
    'WriteNtuple' : 1, # int
    'OutFileName' : 'testout.root', # str
    'SkipEvt' : 0, # int
    'Seed' : 2131.00, # float
    'Debug' : 0, # int
    'EnergyThreshold' : 0.00100000, # float
    'CalibrECAL' : 1.00000, # float
    'CrystalBarLength' : 312.000, # float
    'AttenuationLength' : 7000.00, # float
    'TimeResolution' : 0.100000, # float
    'MatRefractive' : 2.15000, # float
    'InitalTime' : 2.00000, # float
    'ChargeThresholdFrac' : 0.0500000, # float
    'SimCaloHitCollection' : DataHandle('SimCaloCol', 'R', 'DataWrapper<edm4hep::SimCalorimeterHitCollection>'), # DataHandle
    'CaloHitCollection' : DataHandle('DigiCaloCol', 'W', 'DataWrapper<edm4hep::CalorimeterHitCollection>'), # DataHandle
    'CaloAssociationCollection' : DataHandle('ECALBarrelAssoCol', 'W', 'DataWrapper<edm4hep::MCRecoCaloAssociationCollection>'), # DataHandle
    'CaloMCPAssociationCollection' : DataHandle('ECALBarrelParticleAssoCol', 'W', 'DataWrapper<edm4hep::MCRecoCaloParticleAssociationCollection>'), # DataHandle
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
    'Scale' : """  [FanEcalDigiAlg] """,
    'ReadOutName' : """ Readout name [FanEcalDigiAlg] """,
    'WriteNtuple' : """ Write ntuple [FanEcalDigiAlg] """,
    'OutFileName' : """ Output file name [FanEcalDigiAlg] """,
    'SkipEvt' : """ Skip event [FanEcalDigiAlg] """,
    'Seed' : """ Random Seed [FanEcalDigiAlg] """,
    'Debug' : """ Debug level [FanEcalDigiAlg] """,
    'EnergyThreshold' : """ Energy Threshold (/GeV) [FanEcalDigiAlg] """,
    'CalibrECAL' : """ Calibration coefficients for ECAL [FanEcalDigiAlg] """,
    'CrystalBarLength' : """ Crystal Bar Length(mm) [FanEcalDigiAlg] """,
    'AttenuationLength' : """ Crystal Attenuation Length(mm) [FanEcalDigiAlg] """,
    'TimeResolution' : """ Crystal time resolution in one side (ns) [FanEcalDigiAlg] """,
    'MatRefractive' : """ Material refractive index of crystal [FanEcalDigiAlg] """,
    'InitalTime' : """ Start time (ns) [FanEcalDigiAlg] """,
    'ChargeThresholdFrac' : """ Charge threshold fraction [FanEcalDigiAlg] """,
    'SimCaloHitCollection' : """ Handle of the Input SimCaloHit collection [unknown owner type] """,
    'CaloHitCollection' : """ Handle of Digi CaloHit collection [unknown owner type] """,
    'CaloAssociationCollection' : """ Handle of CaloAssociation collection [unknown owner type] """,
    'CaloMCPAssociationCollection' : """ Handle of CaloAssociation collection [unknown owner type] """,
  }
  __declaration_location__ = 'FanEcalDigiAlg.cpp:35'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(FanEcalDigiAlg, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'FanEcalDigi'
  def getType( self ):
      return 'FanEcalDigiAlg'
  pass # class FanEcalDigiAlg
