#Sun Mar 10 20:53:57 2024"""Automatically generated. DO NOT EDIT please"""
import sys
if sys.version_info >= (3,):
    # Python 2 compatibility
    long = int
from GaudiKernel.Proxy.Configurable import *

class G2CDArborAlg( ConfigurableAlgorithm ) :
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
    'ECALCollections' : [ 'EcalBarrelSiliconCollection' , 'EcalEndcapSiliconCollection' , 'EcalEndcapRingCollection' ], # list
    'HcalHitCollections' : [ 'EcalBarrelSiliconPreShowerCollection' , 'EcalEndcapRingPreShowerCollection' , 'EcalEndcapSiliconPreShowerCollection' ], # list
    'HCALCollections' : [ 'HcalBarrelRegCollection' , 'HcalEndcapsCollection' , 'HcalEndcapRingsCollection' ], # list
    'ECALReadOutNames' : [ 'EcalBarrelCollection' , 'EcalEndcapsCollection' , 'EcalEndcapRingCollection' ], # list
    'HCALReadOutNames' : [ 'HcalBarrelCollection' , 'HcalEndcapsCollection' , 'HcalEndcapRingCollection' ], # list
    'DigiECALCollection' : [ 'ECALBarrel' , 'ECALEndcap' , 'ECALOther' ], # list
    'DigiHCALCollection' : [ 'HCALBarrel' , 'HCALEndcap' , 'HCALOther' ], # list
    'ChargeSpatialDistribution' : [ 0.100000 , 0.200000 , 0.400000 , 0.200000 , 0.100000 ], # list
    'CalibrECAL' : [ 40.9100 , 81.8100 ], # list
    'PositionShiftID' : [ 0.00000 , 0.00000 , 0.00000 ], # list
    'ReadLCIO' : False, # bool
    'EventReportEvery' : 100, # int
    'NumThinEcalLayer' : 20, # int
    'ECALThreshold' : 5.00000e-05, # float
    'HCALThreshold' : 0.110000, # float
    'DigiCellSize' : 10, # int
    'ShiftInX' : 0.00000, # float
    'UsingDefaultDetector' : 0, # int
    'PolyaParaA' : 0.700000, # float
    'PolyaParaB' : 0.0450000, # float
    'PolyaParaC' : 0.0300000, # float
    'ChanceOfKink' : 0.00000, # float
    'KinkHitChargeBoost' : 1.00000, # float
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
    'ECALCollections' : """ Input ECAL Hits Collection Names [G2CDArborAlg] """,
    'HcalHitCollections' : """ Hit Collection Names [G2CDArborAlg] """,
    'HCALCollections' : """ HCAL Collection Names [G2CDArborAlg] """,
    'ECALReadOutNames' : """ Name of readouts [G2CDArborAlg] """,
    'HCALReadOutNames' : """ Name of readouts [G2CDArborAlg] """,
    'DigiECALCollection' : """ Name of Digitized ECAL Hit Collections [G2CDArborAlg] """,
    'DigiHCALCollection' : """ Name of Digitized HCAL Hit Collections [G2CDArborAlg] """,
    'ChargeSpatialDistribution' : """ Spactial Distribution of MIP charge X*Y; [G2CDArborAlg] """,
    'CalibrECAL' : """ Calibration coefficients for ECAL [G2CDArborAlg] """,
    'PositionShiftID' : """ Global Position Shift For Overlay [G2CDArborAlg] """,
    'ReadLCIO' : """ Read sim file with LCIO [G2CDArborAlg] """,
    'EventReportEvery' : """ Event ID report every [G2CDArborAlg] """,
    'NumThinEcalLayer' : """ Num of thiner Ecal layers [G2CDArborAlg] """,
    'ECALThreshold' : """ Threshold for ECAL Hits in GeV [G2CDArborAlg] """,
    'HCALThreshold' : """ Threshold for HCAL Hits in GeV [G2CDArborAlg] """,
    'DigiCellSize' : """ Size of Digitized Cell (in mm) [G2CDArborAlg] """,
    'ShiftInX' : """ Shift Distance in X directoin (in mm) NP only [G2CDArborAlg] """,
    'UsingDefaultDetector' : """ Flag Parameter Setting (0 ~ self definition, 1 ~ MircoMegas, 2 ~ GRPC_PS, 3 ~ GRPC_SPS) [G2CDArborAlg] """,
    'PolyaParaA' : """ Polya: x^A*exp(-b*x) + c [G2CDArborAlg] """,
    'PolyaParaB' : """ Polya: x^a*exp(-B*x) + c [G2CDArborAlg] """,
    'PolyaParaC' : """ Polya: x^a*exp(-b*x) + C [G2CDArborAlg] """,
    'ChanceOfKink' : """ Chance of one boundary hit to create a multiple hit with boosted charge [G2CDArborAlg] """,
    'KinkHitChargeBoost' : """ Scale factor of Charge on boosted multiple hits [G2CDArborAlg] """,
  }
  __declaration_location__ = 'G2CDArborAlg.cpp:31'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(G2CDArborAlg, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'G2CDArbor'
  def getType( self ):
      return 'G2CDArborAlg'
  pass # class G2CDArborAlg
