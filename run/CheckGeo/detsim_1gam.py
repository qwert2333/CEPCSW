#!/usr/bin/env python

import os
print(os.environ["DD4HEP_LIBRARY_PATH"])
import sys
# sys.exit(0)

from Gaudi.Configuration import *

##############################################################################
# Random Number Svc
##############################################################################
from Configurables import RndmGenSvc, HepRndm__Engine_CLHEP__RanluxEngine_

# rndmengine = HepRndm__Engine_CLHEP__RanluxEngine_() # The default engine in Gaudi
rndmengine = HepRndm__Engine_CLHEP__HepJamesRandom_() # The default engine in Geant4
rndmengine.SetSingleton = True
rndmengine.Seeds = [39]

# rndmgensvc = RndmGenSvc("RndmGenSvc")
# rndmgensvc.Engine = rndmengine.name()


##############################################################################
# Event Data Svc
##############################################################################
from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")


##############################################################################
# Geometry Svc
##############################################################################

#geometry_path = "/cefs/higgs/guofy/cepcsoft/CEPCSW_dev/CEPCSW/Detector/DetCRD/compact/ecalBarrel.xml"
geometry_path = "/cefs/higgs/guofy/cepcsoft/CEPCSW_dev/CEPCSW_dev02/Detector/DetCRD/compact/CRD_ECAL/CepC_v4-onlyECAL.xml"
#geometry_path = "/cefs/higgs/guofy/cepcsoft/CEPCSW_dev/CEPCSW_dev02/Detector/DetCRD/compact/CRD_ECAL/CepC_v4.xml"
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


gun = GtGunTool("GtGunTool")
gun.Particles = ["gamma"]
gun.PositionXs = [0.] # mm
gun.PositionYs = [0.] # mm
gun.PositionZs = [0.] # mm
gun.EnergyMins = [10.] # GeV
gun.EnergyMaxs = [10.] # GeV
gun.ThetaMins = [90.] # rad; 45deg
gun.ThetaMaxs = [90.] # rad; 45deg
gun.PhiMins = [0.] # rad; 0deg
gun.PhiMaxs = [0.] # rad; 360deg


# stdheprdr = StdHepRdr("StdHepRdr")
# stdheprdr.Input = "/cefs/data/stdhep/CEPC250/2fermions/E250.Pbhabha.e0.p0.whizard195/bhabha.e0.p0.00001.stdhep"

# lciordr = SLCIORdr("SLCIORdr")
# lciordr.Input = "/cefs/data/stdhep/lcio250/signal/Higgs/E250.Pbbh.whizard195/E250.Pbbh_X.e0.p0.whizard195/Pbbh_X.e0.p0.00001.slcio"

# hepmcrdr = HepMCRdr("HepMCRdr")
# hepmcrdr.Input = "example_UsingIterators.txt"

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

# from Configurables import ExampleAnaElemTool
# example_anatool = ExampleAnaElemTool("ExampleAnaElemTool")

from Configurables import DetSimAlg

detsimalg = DetSimAlg("DetSimAlg")

#detsimalg.VisMacs = ["vis.mac"]
detsimalg.RunMacs = ["run.mac"]

detsimalg.RunCmds = [
   #"/run/setCutForAGivenParticle e- 1 cm",
   #"/run/setCut  1 um"
]
detsimalg.AnaElems = [
    # example_anatool.name()
    # "ExampleAnaElemTool"
    "Edm4hepWriterAnaElemTool"
]
detsimalg.RootDetElem = "WorldDetElemTool"
detsimalg.PhysicsList = "QGSP_BERT_EMV"


from Configurables import AnExampleDetElemTool
example_dettool = AnExampleDetElemTool("AnExampleDetElemTool")


##############################################################################
# POD I/O
##############################################################################
from Configurables import PodioOutput
out = PodioOutput("outputalg")
out.filename = "CRDGlobalTest_1gamPFA.root"
out.outputCommands = ["keep *"]

##############################################################################
# ApplicationMgr
##############################################################################

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [genalg, detsimalg, out],
                EvtSel = 'NONE',
                EvtMax = 1,
                ExtSvc = [rndmengine, dsvc, geosvc],
					 #OutputLevel=DEBUG
)
