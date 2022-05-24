from Gaudi.Configuration import *

############## GeomSvc #################
geometry_option = "CRD_o1_v02/CRD_o1_v02-onlyEcalB.xml"

if not os.getenv("DETCRDROOT"):
    print("Can't find the geometry. Please setup envvar DETCRDROOT." )
    sys.exit(-1)

geometry_path = os.path.join(os.getenv("DETCRDROOT"), "compact", geometry_option)
if not os.path.exists(geometry_path):
    print("Can't find the compact geometry file: %s"%geometry_path)
    sys.exit(-1)

from Configurables import GeomSvc
geomsvc = GeomSvc("GeoSvc")
geomsvc.compact = geometry_path
#######################################

########### k4DataSvc ####################
from Configurables import k4DataSvc
podioevent = k4DataSvc("EventDataSvc")
podioevent.inputs = [
"Sim_Ecal.root"
#"/scratchfs/atlas/guofy/CEPCEcalSim/TwoParticles/ForGhost/Sim_GamGamGhost_50mm.root"
#"/scratchfs/atlas/guofy/CEPCEcalSim/TwoParticles/Diphoton/GamGam_50mm.root"
#"/scratchfs/atlas/guofy/CEPCEcalSim/TwoParticles/HadGam/GamHad_100mm.root"
]
##########################################

########### CRDEcalEdmSvc ################
from Configurables import CRDEcalSvc
Ecaldatasvc = CRDEcalSvc("CRDEcalSvc")


########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
#inp.collections = ["EcalBarrelCollection", "EcalEndcapRingCollection", "EcalEndcapsCollection", "MCParticle"]
inp.collections = ["EcalBarrelCollection", "MCParticle"]
##########################################

########## Digitalization ################
from Configurables import CRDEcalDigiAlg
EcalDigi = CRDEcalDigiAlg("CRDEcalDigiAlg")
EcalDigi.ReadOutName = "EcalBarrelCollection"
EcalDigi.SimCaloHitCollection = "EcalBarrelCollection"
EcalDigi.CaloHitCollection = "ECALBarrel"
EcalDigi.CaloAssociationCollection = "RecoCaloAssociation_ECALBarrel"
#EcalDigi.SkipEvt = 63
EcalDigi.Seed = 2079
#Digitalization parameters
EcalDigi.CalibrECAL = 1
EcalDigi.AttenuationLength = 7e10
EcalDigi.TimeResolution = 0.5        #unit: ns
EcalDigi.EnergyThreshold = 0.0001   #0.1 MeV
EcalDigi.ChargeThresholdFrac = 0.05
EcalDigi.Debug=1
EcalDigi.OutFileName = "testTree_Gam.root"
#########################################



######### Reconstruction ################
from Configurables import PandoraPlusPFAlg
PandoraPlusPFAlg = PandoraPlusPFAlg("PandoraPlusPFAlg")
PandoraPlusPFAlg.Seed = 1024
#PandoraPlusPFAlg.SkipEvt = 1
PandoraPlusPFAlg.Debug = 0
PandoraPlusPFAlg.OutFileName = "testRec_Gam.root"
########################################


##############################################################################
# POD I/O
##############################################################################
from Configurables import PodioOutput
out = PodioOutput("outputalg")
out.filename = "Digi_Gam.root"
out.outputCommands = ["keep *"]


########################################

from Configurables import ApplicationMgr
ApplicationMgr( 
    TopAlg=[inp, EcalDigi, PandoraPlusPFAlg],
    EvtSel="NONE",
    EvtMax=10,
    ExtSvc=[podioevent, geomsvc, Ecaldatasvc],
    #OutputLevel=DEBUG
)

