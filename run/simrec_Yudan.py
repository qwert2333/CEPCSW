#!/usr/bin/env python
import os
from Gaudi.Configuration import *


##############################################################################
# Random Number Svc
##############################################################################
from Configurables import RndmGenSvc, HepRndm__Engine_CLHEP__RanluxEngine_

seed = [1]

# rndmengine = HepRndm__Engine_CLHEP__RanluxEngine_() # The default engine in Gaudi
rndmengine = HepRndm__Engine_CLHEP__HepJamesRandom_("RndmGenSvc.Engine") # The default engine in Geant4
rndmengine.SetSingleton = True
rndmengine.Seeds = seed

rndmgensvc = RndmGenSvc("RndmGenSvc")
rndmgensvc.Engine = rndmengine.name()


##############################################################################

from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")

from Configurables import MarlinEvtSeeder
evtseeder = MarlinEvtSeeder("EventSeeder")

#geometry_option = "CepC_v4-onlyTrackerECAL.xml"
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

stdheprdr = StdHepRdr("StdHepRdr") 
stdheprdr.Input = "/cefs/data/stdhep/CEPC240/higgs/Higgs_10M/data/E240.Pnnh_gg.e0.p0.whizard195/nnh_gg.e0.p0.00001.stdhep"

genprinter = GenPrinter("GenPrinter")

genalg = GenAlgo("GenAlgo")
genalg.GenTools = ["StdHepRdr"]

from Configurables import DetSimSvc
detsimsvc = DetSimSvc("DetSimSvc")

from Configurables import DetSimAlg
detsimalg = DetSimAlg("DetSimAlg")
detsimalg.RandomSeeds = seed
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

from Configurables import GearSvc
gearsvc = GearSvc("GearSvc")
#gearsvc.GearXMLFile = "/cefs/higgs/fucd/rec/GearOutput_v4.xml"
gearsvc.GearXMLFile = "../Detector/DetCEPCv4/compact/FullDetGear.xml"


from Configurables import TrackSystemSvc
tracksystemsvc = TrackSystemSvc("TrackSystemSvc")

vxdhitname  = "VXDTrackerHits"
sithitname  = "SITTrackerHits"
sitspname   = "SITSpacePoints"
tpchitname  = "TPCTrackerHits"
sethitname  = "SETTrackerHits"
setspname   = "SETSpacePoints"
ftdphitname = "FTDPixelTrackerHits"
ftdshitname = "FTDStripTrackerHits"
ftdspname   = "FTDSpacePoints"
from Configurables import PlanarDigiAlg
digiVXD = PlanarDigiAlg("VXDDigi")
digiVXD.SimTrackHitCollection = "VXDCollection"
digiVXD.TrackerHitCollection = vxdhitname
digiVXD.ResolutionU = [0.0028, 0.006, 0.004, 0.004, 0.004, 0.004]
digiVXD.ResolutionV = [0.0028, 0.006, 0.004, 0.004, 0.004, 0.004]

digiSIT = PlanarDigiAlg("SITDigi")
digiSIT.IsStrip = 1
digiSIT.SimTrackHitCollection = "SITCollection"
digiSIT.TrackerHitCollection = sithitname
digiSIT.TrackerHitAssociationCollection = "SITTrackerHitAssociation"
digiSIT.ResolutionU = [0.007]
digiSIT.ResolutionV = [0.000]

digiSET = PlanarDigiAlg("SETDigi")
digiSET.IsStrip = 1
digiSET.SimTrackHitCollection = "SETCollection"
digiSET.TrackerHitCollection = sethitname
digiSET.TrackerHitAssociationCollection = "SETTrackerHitAssociation"
digiSET.ResolutionU = [0.007]
digiSET.ResolutionV = [0.000]

digiFTDP = PlanarDigiAlg("FTDPixelDigi")
digiFTDP.SimTrackHitCollection = "FTD_PIXELCollection"
digiFTDP.TrackerHitCollection = ftdphitname
digiFTDP.TrackerHitAssociationCollection = "FTDPixelTrackerHitAssociation"
digiFTDP.ResolutionU = [0.003]
digiFTDP.ResolutionV = [0.003]

digiFTDS = PlanarDigiAlg("FTDStripDigi")
digiFTDS.IsStrip = 1
digiFTDS.SimTrackHitCollection = "FTD_STRIPCollection"
digiFTDS.TrackerHitCollection = ftdshitname
digiFTDS.TrackerHitAssociationCollection = "FTDStripTrackerHitAssociation"
digiFTDS.ResolutionU = [0.007]
digiFTDS.ResolutionV = [0.000]

from Configurables import SpacePointBuilderAlg
spSIT = SpacePointBuilderAlg("SITBuilder")
spSIT.TrackerHitCollection = sithitname
spSIT.TrackerHitAssociationCollection = "SITTrackerHitAssociation"
spSIT.SpacePointCollection = sitspname
spSIT.SpacePointAssociationCollection = "SITSpacePointAssociation"
#spSIT.OutputLevel = DEBUG

spSET = SpacePointBuilderAlg("SETBuilder")
spSET.TrackerHitCollection = sethitname
spSET.TrackerHitAssociationCollection = "SETTrackerHitAssociation"
spSET.SpacePointCollection = setspname
spSET.SpacePointAssociationCollection = "SETSpacePointAssociation"
#spSET.OutputLevel = DEBUG

spFTD = SpacePointBuilderAlg("FTDBuilder")
spFTD.TrackerHitCollection = ftdshitname
spFTD.TrackerHitAssociationCollection = "FTDStripTrackerHitAssociation"
spFTD.SpacePointCollection = ftdspname
spFTD.SpacePointAssociationCollection = "FTDSpacePointAssociation"
#spFTD.OutputLevel = DEBUG

from Configurables import TPCDigiAlg
digiTPC = TPCDigiAlg("TPCDigi")
digiTPC.TPCCollection = "TPCCollection"
digiTPC.TPCLowPtCollection = "TPCLowPtCollection"
digiTPC.TPCTrackerHitsCol = tpchitname
digiTPC.TPCTrackerHitAssCol = 'TPCTrackerHitAssociation'
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
tracking.FTDPixelHitCollection = ftdphitname
tracking.FTDSpacePointCollection = ftdspname
tracking.SITRawHitCollection = sithitname
tracking.FTDRawHitCollection = ftdshitname
tracking.UseSIT = 1
tracking.SmoothOn = 0
#tracking.OutputLevel = DEBUG

from Configurables import ForwardTrackingAlg
forward = ForwardTrackingAlg("ForwardTracking")
forward.FTDPixelHitCollection = ftdphitname
forward.FTDSpacePointCollection = ftdspname
forward.FTDRawHitCollection = ftdshitname
forward.Chi2ProbCut = 0.0
forward.HitsPerTrackMin = 3
forward.BestSubsetFinder = "SubsetSimple"
forward.Criteria = ["Crit2_DeltaPhi","Crit2_StraightTrackRatio","Crit3_3DAngle","Crit3_ChangeRZRatio","Crit3_IPCircleDist","Crit4_3DAngleChange","Crit4_DistToExtrapolation",
                    "Crit2_DeltaRho","Crit2_RZRatio","Crit3_PT"]
forward.CriteriaMin = [0,  0.9,  0,  0.995, 0,  0.8, 0,   20,  1.002, 0.1,      0,   0.99, 0,    0.999, 0,   0.99, 0]  
forward.CriteriaMax = [30, 1.02, 10, 1.015, 20, 1.3, 1.0, 150, 1.08,  99999999, 0.8, 1.01, 0.35, 1.001, 1.5, 1.01, 0.05] 
#forward.OutputLevel = DEBUG

from Configurables import TrackSubsetAlg
subset = TrackSubsetAlg("TrackSubset")
subset.TrackInputCollections = ["ForwardTracks", "SiTracks"]
subset.RawTrackerHitCollections = [vxdhitname, sithitname, ftdphitname, ftdshitname, sitspname, ftdspname]
subset.TrackSubsetCollection = "SubsetTracks"

from Configurables import FullLDCTrackingAlg
full = FullLDCTrackingAlg("FullTracking")
full.VTXTrackerHits = vxdhitname
full.SITTrackerHits = sitspname
full.TPCTrackerHits = tpchitname
full.SETTrackerHits = setspname
full.FTDPixelTrackerHits = ftdphitname
full.FTDSpacePoints = ftdspname
full.SITRawHits     = sithitname
full.SETRawHits     = sethitname
full.FTDRawHits     = ftdshitname
full.TPCTracks = "ClupatraTracks"
full.SiTracks  = "SubsetTracks"
full.OutputTracks  = "MarlinTrkTracks"
#full.OutputLevel = DEBUG

from Configurables import TrackInspectAlg
trackEff = TrackInspectAlg('TrackingEfficiency')
trackEff.TrackCollection = 'ClupatraTracks'
trackEff.MCParticleCollection = 'MCParticle'
trackEff.TPCTrackerHitRelations = 'TPCTrackerHitAssociation'
trackEff.VXDTrackerHitRelations = 'VXDTrackerHitAssociation'
trackEff.SITTrackerHitRelations = 'SITTrackerHitAssociation'
trackEff.SETTrackerHitRelations = 'SETTrackerHitAssociation'
trackEff.FTDTrackerHitRelations = 'FTDTrackerHitAssociation'
trackEff.OutputLevel = DEBUG


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
from Configurables import SimHitMergeAlg
simHitMerge = SimHitMergeAlg("SimHitMergeAlg")
simHitMerge.sanity_check = False
simHitMerge.InputCollections=["EcalBarrelCollection", "EcalEndcapsCollection","EcalEndcapRingCollection", "HcalBarrelCollection", "HcalEndcapsCollection", "HcalEndcapRingCollection"]
simHitMerge.OutputCollections=["EcalBarrelCollectionMerged", "EcalEndcapsCollectionMerged", "EcalEndcapRingCollectionMerged", "HcalBarrelCollectionMerged", "HcalEndcapsCollectionMerged", "HcalEndcapRingCollectionMerged"]

############################################################

from Configurables import G2CDArborAlg
caloDigi = G2CDArborAlg("G2CDArborAlg")
#caloDigi.ReadLCIO = False 
#caloDigi.CalibrECAL = [48.16, 96.32]
caloDigi.CalibrECAL = [46.538, 93.0769]
caloDigi.ECALCollections = ["EcalBarrelCollection", "EcalEndcapsCollection","EcalEndcapRingCollection"]
caloDigi.ECALReadOutNames= ["EcalBarrelCollection", "EcalEndcapsCollection", "EcalEndcapRingCollection"]
caloDigi.DigiECALCollection = ["ECALBarrel", "ECALEndcap", "EcalEndcapRing"]
caloDigi.HCALCollections = ["HcalBarrelCollection", "HcalEndcapsCollection", "HcalEndcapRingCollection"]
caloDigi.HCALReadOutNames= ["HcalBarrelCollection", "HcalEndcapsCollection", "HcalEndcapRingCollection"]
caloDigi.DigiHCALCollection = ["HCALBarrel", "HCALEndcap", "HCALOther"]
caloDigi.EventReportEvery = 1

##############################################################################
from Configurables import MarlinArbor

marlinArbor = MarlinArbor("MarlinArbor")

marlinArbor.ECALCollections=  ["ECALBarrel"          , "ECALEndcap","EcalEndcapRing"   ]
marlinArbor.ECALReadOutNames= ["EcalBarrelCollection", "EcalEndcapsCollection","EcalEndcapRingCollection"]
marlinArbor.HCALCollections=  ["HCALBarrel", "HCALEndcap","HCALOther"]
marlinArbor.HCALReadOutNames= ["HcalBarrelCollection","HcalEndcapsCollection","HcalEndcapRingCollection"]
#marlinArbor.ReadLCIO = False
##############################################################################

from Configurables import BushConnect
bushconnect = BushConnect("BushConnect")

from Configurables import ClusterAna
clusterAna = ClusterAna("ClusterAna")
clusterAna.TreeOutputFile="Ana_nnhgg_00001.root"

from Configurables import TotalInvMass
total_inv_mass = TotalInvMass("TotalInvMass")
total_inv_mass.TreeOutputFile="BMR_00001.root"

# write PODIO file
from Configurables import PodioOutput
write = PodioOutput("write")
write.filename = "sim_rec_nnhgg_00001.root"
write.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg = [genalg, detsimalg, digiVXD, digiSIT, digiSET, digiFTDP, digiFTDS, spSIT, spSET, spFTD, digiTPC, clupatra, tracking, forward, subset, full, simHitMerge, caloDigi, marlinArbor, bushconnect, clusterAna, total_inv_mass, write],
    #TopAlg = [genalg, detsimalg, digiVXD, digiSIT, digiSET, digiFTDP, digiFTDS, spSIT, spSET, spFTD, digiTPC, clupatra, tracking, forward, subset, full, write],
    EvtSel = 'NONE',
    EvtMax = 1,
    ExtSvc = [rndmengine, rndmgensvc, dsvc, evtseeder, gearsvc, geosvc, tracksystemsvc],
    HistogramPersistency='ROOT',
    OutputLevel=INFO
)
