import os
import sys

from Gaudi.Configuration import *

############## GeomSvc #################
geometry_option = "/afs/ihep.ac.cn/users/h/huangliheng/script/CEPCSW-master/Detector/DetCRD/compact/Standalone/Standalone-EcalRotCrystal.xml"
#...

#if not os.getenv("DETCRDROOT"):
#    print("Can't find the geometry. Please setup envvar DETCRDROOT." )
#    sys.exit(-1)

#geometry_path = os.path.join(os.getenv("DETCRDROOT"), "compact", geometry_option)

geometry_path = geometry_option

if not os.path.exists(geometry_path):
    print("Can't find the compact geometry file: %s"%geometry_path)
    sys.exit(-1)

from Configurables import GeomSvc
geomsvc = GeomSvc("GeomSvc")
geomsvc.compact = geometry_path
#######################################

########### k4DataSvc ####################
from Configurables import k4DataSvc
podioevent = k4DataSvc("EventDataSvc", input = "/afs/ihep.ac.cn/users/h/huangliheng/script/Data/SimData/SimSCECal_gpi_21.root")##########################################


########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
#inp.collections = ["EcalBarrelCollection", "EcalEndcapRingCollection", "EcalEndcapsCollection", "MCParticle"]
inp.collections = ["EcalBarrelCollection", "MCParticleGen"]
##########################################

########## Digitalization ################
from Configurables import FanEcalDigiAlg
EcalDigi = FanEcalDigiAlg("FanEcalDigiAlg")
EcalDigi.ReadOutName = "EcalBarrelCollection"
EcalDigi.SimCaloHitCollection = "EcalBarrelCollection"
EcalDigi.CaloHitCollection = "ECALBarrel"
EcalDigi.CaloAssociationCollection = "RecoCaloAssociation_ECALBarrel"
#EcalDigi.SkipEvt = 3
EcalDigi.Seed = 2079
#Digitalization parameters
EcalDigi.CalibrECAL = 1
EcalDigi.CrystalBarLength = 338 
EcalDigi.AttenuationLength = 7e10
EcalDigi.TimeResolution = 0.5        #unit: ns
EcalDigi.EnergyThreshold = 0.001   #0.1 MeV  0.0001
EcalDigi.ChargeThresholdFrac = 0.002
EcalDigi.Debug=0
EcalDigi.OutFileName = "/afs/ihep.ac.cn/users/h/huangliheng/script/Data/DigiData/Rec_SCECAL_gpi.root"
#########################################


########################################

from Configurables import ApplicationMgr
ApplicationMgr( 
    TopAlg=[inp, EcalDigi],
    EvtSel="NONE",
    EvtMax=100,
    ExtSvc=[podioevent, geomsvc],
    #OutputLevel=DEBUG
)

