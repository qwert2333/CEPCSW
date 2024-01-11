#!/usr/bin/env python
import os
from Gaudi.Configuration import *

NTupleSvc().Output = ["MyTuples DATAFILE='sim-rec-trackerEcal.root' OPT='NEW' TYP='ROOT'"]

from Configurables import RndmGenSvc, HepRndm__Engine_CLHEP__RanluxEngine_
rndmengine = HepRndm__Engine_CLHEP__HepJamesRandom_()
rndmengine.SetSingleton = True
rndmengine.Seeds = [1]

from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")

from Configurables import MarlinEvtSeeder
evtseeder = MarlinEvtSeeder("EventSeeder")

geometry_option = "CepC_v4.xml"

if not os.getenv("DETCEPCV4ROOT"):
    print("Can't find the geometry. Please setup envvar DETCEPCV4ROOT." )
    sys.exit(-1)

geometry_path = os.path.join(os.getenv("DETCEPCV4ROOT"), "compact", geometry_option)
if not os.path.exists(geometry_path):
    print("Can't find the compact geometry file: %s"%geometry_path)
    sys.exit(-1)

from Configurables import GeomSvc
geosvc = GeomSvc("GeomSvc")
geosvc.compact = geometry_path

from Configurables import GenAlgo
from Configurables import GtGunTool
from Configurables import StdHepRdr
from Configurables import SLCIORdr
from Configurables import HepMCRdr
from Configurables import GenPrinter

#gun = GtGunTool("GtGunTool")
#gun.Particles = ["pi-","pi-","pi-","pi-","pi-"]
#gun.PositionXs = [0.,0.,0.,0.,0.]
#gun.PositionYs = [0.,0.,0.,0.,0.]
#gun.PositionZs = [0.,0.,0.,0.,0.]
#gun.EnergyMins = [2.,2.,2.,2.,2.] # GeV
#gun.EnergyMaxs = [20.,20.,20.,20.,20.] # GeV
#gun.ThetaMins  = [15.,15.,15.,15.,15.]   # deg
#gun.ThetaMaxs  = [165.,165.,165.,165.,165.]   # deg
#gun.PhiMins    = [0.,0.,0.,0.,0.]   # deg
#gun.PhiMaxs    = [360.,360.,360.,360.,360.]   # deg

stdheprdr = StdHepRdr("StdHepRdr")
stdheprdr.Input = "/cefs/data/stdhep/CEPC240/higgs/Higgs_10M/data/E240.Pnnh_gg.e0.p0.whizard195/nnh_gg.e0.p0.00001.stdhep"

genprinter = GenPrinter("GenPrinter")

genalg = GenAlgo("GenAlgo")
#genalg.GenTools = ["GtGunTool"]
genalg.GenTools = ["StdHepRdr"]

from Configurables import DetSimSvc
detsimsvc = DetSimSvc("DetSimSvc")

from Configurables import DetSimAlg
detsimalg = DetSimAlg("DetSimAlg")
# detsimalg.VisMacs = ["vis.mac"]
detsimalg.RunCmds = [
#    "/physics_lists/factory/addOptical"
]
detsimalg.PhysicsList = "FTFP_BERT"
detsimalg.AnaElems = ["Edm4hepWriterAnaElemTool"]
detsimalg.RootDetElem = "WorldDetElemTool"

from Configurables import AnExampleDetElemTool
example_dettool = AnExampleDetElemTool("AnExampleDetElemTool")

from Configurables import TimeProjectionChamberSensDetTool
tpc_sensdettool = TimeProjectionChamberSensDetTool("TimeProjectionChamberSensDetTool")
tpc_sensdettool.TypeOption = 1


from Configurables import CalorimeterSensDetTool
from Configurables import DriftChamberSensDetTool
cal_sensdettool = CalorimeterSensDetTool("CalorimeterSensDetTool")
#cal_sensdettool.CalNamesMergeDisable = ["CaloDetector"]
cal_sensdettool.CalNamesApplyBirks = ["HcalBarrel"]


from Configurables import GearSvc
gearsvc = GearSvc("GearSvc")
#gearsvc.GearXMLFile = "Detector/DetCEPCv4/compact/FullDetGear.xml"

from Configurables import TrackSystemSvc
tracksystemsvc = TrackSystemSvc("TrackSystemSvc")

vxdhitname  = "VXDTrackerHits"
sithitname  = "SITTrackerHits"
sitspname   = "SITSpacePoints"
tpchitname  = "TPCTrackerHits"
sethitname  = "SETTrackerHits"
setspname   = "SETSpacePoints"
ftdspname   = "FTDSpacePoints"
ftdhitname = "FTDTrackerHits"
from Configurables import PlanarDigiAlg
digiVXD = PlanarDigiAlg("VXDDigi")
digiVXD.SimTrackHitCollection = "VXDCollection"
digiVXD.TrackerHitCollection = vxdhitname
digiVXD.TrackerHitAssociationCollection = "VXDTrackerHitAssociation"
digiVXD.ResolutionU = [0.0028, 0.006, 0.004, 0.004, 0.004, 0.004]
digiVXD.ResolutionV = [0.0028, 0.006, 0.004, 0.004, 0.004, 0.004]

digiSIT = PlanarDigiAlg("SITDigi")
digiSIT.IsStrip = 1
digiSIT.SimTrackHitCollection = "SITCollection"
digiSIT.TrackerHitCollection = sithitname
digiSIT.TrackerHitAssociationCollection = "SITTrackerHitAssociation"
digiSIT.ResolutionU = [0.007]
digiSIT.ResolutionV = [0.000]
#digiSIT.UsePlanarTag = True

digiSET = PlanarDigiAlg("SETDigi")
digiSET.IsStrip = 1
digiSET.SimTrackHitCollection = "SETCollection"
digiSET.TrackerHitCollection = sethitname
digiSET.TrackerHitAssociationCollection = "SETTrackerHitAssociation"
digiSET.ResolutionU = [0.007]
digiSET.ResolutionV = [0.000]
#digiSET.UsePlanarTag = True

digiFTD = PlanarDigiAlg("FTDDigi")
digiFTD.SimTrackHitCollection = "FTDCollection"
digiFTD.TrackerHitCollection = ftdhitname
digiFTD.TrackerHitAssociationCollection = "FTDTrackerHitAssociation"
digiFTD.ResolutionU = [0.003, 0.003, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007]
digiFTD.ResolutionV = [0.003, 0.003, 0,     0,     0,     0,     0,     0    ]
#digiFTD.UsePlanarTag = True
#digiFTD.OutputLevel = DEBUG

from Configurables import SpacePointBuilderAlg
spSIT = SpacePointBuilderAlg("SITBuilder")
spSIT.TrackerHitCollection = sithitname
spSIT.TrackerHitAssociationCollection = "SITTrackerHitAssociation"
spSIT.SpacePointCollection = sitspname
spSIT.SpacePointAssociationCollection = "SITSpacePointAssociation"
#spSIT.OutputLevel = DEBUG

spFTD = SpacePointBuilderAlg("FTDBuilder")
spFTD.TrackerHitCollection = ftdhitname
spFTD.TrackerHitAssociationCollection = "FTDTrackerHitAssociation"
spFTD.SpacePointCollection = ftdspname
spFTD.SpacePointAssociationCollection = "FTDSpacePointAssociation"
#spFTD.OutputLevel = DEBUG

from Configurables import TPCDigiAlg
digiTPC = TPCDigiAlg("TPCDigi")
digiTPC.TPCCollection = "TPCCollection"
digiTPC.TPCLowPtCollection = "TPCLowPtCollection"
digiTPC.TPCTrackerHitsCol = tpchitname
digiTPC.TPCTrackerHitAssCol = "TPCTrackerHitAssociation"
#digiTPC.OutputLevel = DEBUG

from Configurables import ClupatraAlg
clupatra = ClupatraAlg("Clupatra")
clupatra.TPCHitCollection = tpchitname
#clupatra.OutputLevel = DEBUG

from Configurables import SiliconTrackingAlg
tracking = SiliconTrackingAlg("SiliconTracking")
tracking.HeaderCol = "EventHeader"
tracking.VTXHitCollection = vxdhitname
tracking.SITHitCollection = sitspname
tracking.FTDPixelHitCollection = ftdhitname
tracking.FTDSpacePointCollection = ftdspname
tracking.SITRawHitCollection = sithitname
tracking.FTDRawHitCollection = ftdhitname
tracking.UseSIT = 1
tracking.SmoothOn = 0
tracking.DumpTime = False
#tracking.NDivisionsInTheta = 10
#tracking.OutputLevel = DEBUG

from Configurables import ForwardTrackingAlg
forward = ForwardTrackingAlg("ForwardTracking")
forward.FTDPixelHitCollection = ftdhitname
forward.FTDSpacePointCollection = ftdspname
forward.FTDRawHitCollection = ftdhitname
forward.Chi2ProbCut = 0.0
forward.HitsPerTrackMin = 3
forward.BestSubsetFinder = "SubsetSimple"
forward.Criteria = ["Crit2_DeltaPhi","Crit2_StraightTrackRatio","Crit3_3DAngle","Crit3_ChangeRZRatio","Crit3_IPCircleDist","Crit4_3DAngleChange","Crit4_DistToExtrapolation",
                    "Crit2_DeltaRho","Crit2_RZRatio","Crit3_PT"]
forward.CriteriaMin = [0,  0.9,  0,  0.995, 0,  0.8, 0,   20,  1.002, 0.1,      0,   0.99, 0,    0.999, 0,   0.99, 0]  
forward.CriteriaMax = [30, 1.02, 10, 1.015, 20, 1.3, 1.0, 150, 1.08,  99999999, 0.8, 1.01, 0.35, 1.001, 1.5, 1.01, 0.05] 
forward.DumpTime = False
#forward.OutputLevel = DEBUG

from Configurables import TrackSubsetAlg
subset = TrackSubsetAlg("TrackSubset")
subset.TrackInputCollections = ["ForwardTracks", "SiTracks"]
#subset.RawTrackerHitCollections = [vxdhitname, sithitname, ftdhitname, ftdspname]
subset.RawTrackerHitCollections = [vxdhitname, sithitname, ftdhitname, sitspname, ftdspname]
subset.TrackSubsetCollection = "SubsetTracks"
subset.DumpTime = False
#subset.OutputLevel = DEBUG

from Configurables import FullLDCTrackingAlg
full = FullLDCTrackingAlg("FullTracking")
full.VTXTrackerHits = vxdhitname
full.SITTrackerHits = sitspname
full.TPCTrackerHits = tpchitname
full.SETTrackerHits = setspname
full.FTDPixelTrackerHits = ftdhitname
full.FTDSpacePoints = ftdspname
full.SITRawHits     = sithitname
full.SETRawHits     = sethitname
full.FTDRawHits     = ftdhitname
full.VTXHitRelCol   = "VXDTrackerHitAssociation"
full.SITHitRelCol   = "SITTrackerHitAssociation"
full.SETHitRelCol   = "SETTrackerHitAssociation"
full.FTDHitRelCol   = "FTDTrackerHitAssociation"
full.TPCHitRelCol   = "TPCTrackerHitAssociation"
full.TPCTracks = "ClupatraTracks"
full.SiTracks  = "SubsetTracks"
full.OutputTracks  = "MarlinTrkTracks"
full.DumpTime = False
#full.SITHitToTrackDistance = 3.
#full.SETHitToTrackDistance = 5.
#full.MinChi2ProbForSiliconTracks = 0
#full.OutputLevel = DEBUG
'''
from Configurables import DumpMCParticleAlg
dumpMC = DumpMCParticleAlg("DumpMC")
dumpMC.MCParticleCollection = "MCParticle"

from Configurables import DumpTrackAlg
dumpFu = DumpTrackAlg("DumpFu")
dumpFu.TrackCollection = "MarlinTrkTracks"
#dumpFu.OutputLevel = DEBUG

dumpCl = DumpTrackAlg("DumpCl")
dumpCl.TrackCollection = "ClupatraTracks"
#dumpCl.OutputLevel = DEBUG

dumpSu = DumpTrackAlg("DumpSu")
dumpSu.TrackCollection = "SubsetTracks"
#dumpSu.OutputLevel = DEBUG

dumpSi = DumpTrackAlg("DumpSi")
dumpSi.TrackCollection = "SiTracks"
#dumpSi.OutputLevel = DEBUG

dumpFo = DumpTrackAlg("DumpFo")
dumpFo.TrackCollection = "ForwardTracks"
#dumpFo.OutputLevel = DEBUG
'''
############################################################
#from Configurables import SimHitMergeAlg
#simHitMerge = SimHitMergeAlg("SimHitMergeAlg")
#simHitMerge.InputCollections=["EcalBarrelCollection", "EcalEndcapsCollection"]
#simHitMerge.OutputCollections=["EcalBarrelCollectionMerged", "EcalEndcapsCollectionMerged"]
############################################################

from Configurables import G2CDArborAlg
caloDigi = G2CDArborAlg("G2CDArborAlg")
caloDigi.ReadLCIO = False 
#caloDigi.CalibrECAL = [48.16, 96.32]
caloDigi.CalibrECAL = [46.538, 93.0769]
caloDigi.ECALCollections = ["EcalBarrelCollection", "EcalEndcapsCollection", "EcalEndcapRingCollection"]
caloDigi.ECALReadOutNames= ["EcalBarrelCollection", "EcalEndcapsCollection", "EcalEndcapRingCollection"]
caloDigi.DigiECALCollection = ["ECALBarrel", "ECALEndcap", "ECALOther"]
caloDigi.HCALCollections = ["HcalBarrelCollection", "HcalEndcapsCollection","HcalEndcapRingCollection"]
caloDigi.HCALReadOutNames= ["HcalBarrelCollection", "HcalEndcapsCollection", "HcalEndcapRingCollection"]
caloDigi.DigiHCALCollection = ["HCALBarrel", "HCALEndcap", "HCALOther"]
caloDigi.EventReportEvery = 100
##############################################################################
from Configurables import PandoraPFAlg
from Configurables import MarlinArbor
from Configurables import BushConnect

marlinArbor = MarlinArbor("MarlinArbor")
marlinArbor.ECALCollections = ["ECALBarrel", "ECALEndcap", "ECALOther"]
marlinArbor.ECALReadOutNames = ["EcalBarrelCollection", "EcalEndcapsCollection", "EcalEndcapRingCollection"]
marlinArbor.HCALCollections = ["HCALBarrel", "HCALEndcap", "HCALOther"]
marlinArbor.HCALReadOutNames = ["HcalBarrelCollection", "HcalEndcapsCollection", "HcalEndcapRingCollection"]
bush = BushConnect("BushConnect")

from Configurables import ClusterAna
clusterAna = ClusterAna("ClusterAna")
clusterAna.TreeOutputFile="Ana_nnhgg_00001.root"

from Configurables import TotalInvMass
total_inv_mass = TotalInvMass("TotalInvMass")
total_inv_mass.TreeOutputFile="BMR_00001.root"

##############################################################################

# write PODIO file
from Configurables import PodioOutput
write = PodioOutput("write")
write.filename = "test.root"
write.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg = [genalg, detsimalg, digiVXD, digiSIT, digiSET, digiFTD, spSIT, spFTD, digiTPC, clupatra, tracking, forward, subset, full, caloDigi, marlinArbor, bush, clusterAna, total_inv_mass, write],
    #TopAlg = [genalg, detsimalg, digiVXD, digiSIT, digiSET, digiFTD, spSIT, spFTD, digiTPC, clupatra, tracking, forward, subset, full, simHitMerge, caloDigi, pandoralg, write],
    EvtSel = 'NONE',
    EvtMax = 1,
    ExtSvc = [rndmengine, dsvc, evtseeder, geosvc, gearsvc, tracksystemsvc],
    HistogramPersistency='ROOT',
    OutputLevel=INFO
)
