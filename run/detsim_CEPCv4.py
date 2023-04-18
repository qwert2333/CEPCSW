#!/usr/bin/env python

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

gun = GtGunTool("GtGunTool")
gun.Particles = ["pi-"]
gun.EnergyMins = [10]
gun.EnergyMaxs = [10]
gun.ThetaMins = [90]
gun.ThetaMaxs = [90]
gun.PhiMins = [0.]
gun.PhiMaxs = [360.]

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

from Configurables import GearSvc
gearsvc = GearSvc("GearSvc")
gearsvc.GearXMLFile = "../Detector/DetCEPCv4/compact/FullDetGear.xml"

############################################################
from Configurables import SimHitMergeAlg
simHitMerge = SimHitMergeAlg("SimHitMergeAlg")
simHitMerge.InputCollections=["EcalBarrelCollection", "EcalEndcapsCollection"]
simHitMerge.OutputCollections=["EcalBarrelCollectionMerged", "EcalEndcapsCollectionMerged"]
############################################################


# write PODIO file
from Configurables import PodioOutput
write = PodioOutput("write")
write.filename = "nnHgg.root"
write.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg = [genalg, detsimalg, simHitMerge, write],
    EvtSel = 'NONE',
    EvtMax = 10,
    ExtSvc = [rndmengine, dsvc, evtseeder, gearsvc, geosvc],
    HistogramPersistency='ROOT',
    OutputLevel=INFO
)
