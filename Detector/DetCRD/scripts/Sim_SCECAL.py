import os
import sys
import math
from Gaudi.Configuration import *

##############################################################################
# Random Number Svc
##############################################################################
from Configurables import RndmGenSvc, HepRndm__Engine_CLHEP__RanluxEngine_
seed = [${RNDM}]
# rndmengine = HepRndm__Engine_CLHEP__RanluxEngine_() # The default engine in Gaudi
rndmengine = HepRndm__Engine_CLHEP__HepJamesRandom_("RndmGenSvc.Engine") # The default engine in Geant4
rndmengine.SetSingleton = True
rndmengine.Seeds = seed

rndmgensvc = RndmGenSvc("RndmGenSvc")
rndmgensvc.Engine = rndmengine.name()

##############################################################################
# Event Data Svc
##############################################################################
from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")

##############################################################################
# Geometry Svc
##############################################################################
geometry_option = "Standalone/Standalone-EcalRotCrystal.xml"

if not os.getenv("DETCRDROOT"):
    print("Can't find the geometry. Please setup envvar DETCRDROOT." )
    sys.exit(-1)

geometry_path = os.path.join(os.getenv("DETCRDROOT"), "compact", geometry_option)
if not os.path.exists(geometry_path):
    print("Can't find the compact geometry file: %s"%geometry_path)
    sys.exit(-1)


from Configurables import GeomSvc
geosvc = GeomSvc("GeomSvc")
geosvc.compact = geometry_path

##############################################################################
# Physics Generator
##############################################################################
from Configurables import GenAlgo
from Configurables import GtGunTool
from Configurables import StdHepRdr
from Configurables import SLCIORdr
from Configurables import HepMCRdr
from Configurables import GenPrinter

def delta_phi(d, E):
    L = 1900.0
    B = 3.0
    Z = 1.0
    m = 0.1396
    deflection = math.degrees(d / L)
    magnetic = math.degrees(math.asin(0.00015 * L * Z * B / (math.sqrt(E * E - m * m))))
    return deflection + magnetic

gun = GtGunTool("GtGunTool")
gun.Particles  = ["gamma","pi+"]
gun.EnergyMins = [5.,5.] # GeV
gun.EnergyMaxs = [5.,5.] # GeV
gun.ThetaMins  = [90.,90.]  # deg
gun.ThetaMaxs  = [90.,90.]  # deg
gun.PhiMins    = [100.,100.+delta_phi(d=50,E=5)]  # deg
gun.PhiMaxs    = [100.,100.+delta_phi(d=50,E=5)]  # deg
#gun.PositionXs = [0. ] # mm
#gun.PositionYs = [0. ] # mm
#gun.PositionZs = [0. ] # mm

genprinter = GenPrinter("GenPrinter")

genalg = GenAlgo("GenAlgo")
genalg.GenTools = ["GtGunTool"]
# genalg.GenTools = ["StdHepRdr"]
# genalg.GenTools = ["StdHepRdr", "GenPrinter"]
# genalg.GenTools = ["SLCIORdr", "GenPrinter"]
# genalg.GenTools = ["HepMCRdr", "GenPrinter"]

##############################################################################
# Detector Simulation
##############################################################################

from Configurables import DetSimSvc
detsimsvc = DetSimSvc("DetSimSvc")

from Configurables import DetSimAlg
detsimalg = DetSimAlg("DetSimAlg")
detsimalg.RandomSeeds = seed
#detsimalg.VisMacs = ["vis.mac"]
#detsimalg.RunMacs = ["run.mac"]

detsimalg.RunCmds = [
   #"/run/setCutForAGivenParticle e- 1 cm",
   #"/run/setCut  1 um"
]
detsimalg.AnaElems = [
    # example_anatool.name()
    #"ExampleAnaElemTool"
    "Edm4hepWriterAnaElemTool"
]
detsimalg.RootDetElem = "WorldDetElemTool"
detsimalg.PhysicsList = "QGSP_BERT_EMV"

from Configurables import AnExampleDetElemTool
example_dettool = AnExampleDetElemTool("AnExampleDetElemTool")

from Configurables import DriftChamberSensDetTool
from Configurables import CalorimeterSensDetTool
cal_sensdettool = CalorimeterSensDetTool("CalorimeterSensDetTool")
cal_sensdettool.CalNamesMergeDisable = ["CaloDetector"]


##############################################################################
# POD I/O
##############################################################################
from Configurables import PodioOutput
out = PodioOutput("outputalg")

sample_dir = "test"
output_path = "/afs/ihep.ac.cn/users/h/huangliheng/pub/SimData/" + sample_dir + "/${NUM}"

try:
    os.makedirs(output_path)
except OSError:
    pass

out.filename = output_path + "/Sim_${RNDM}.root"
out.outputCommands = ["keep *"]

##############################################################################
# ApplicationMgr
##############################################################################

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [genalg, detsimalg, out],
                EvtSel = 'NONE',
                EvtMax = 100,
                ExtSvc = [rndmengine, dsvc, geosvc],
		#OutputLevel=DEBUG
)
