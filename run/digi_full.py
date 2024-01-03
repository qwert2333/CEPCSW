from Gaudi.Configuration import *
Nskip = 0
Nevt = 10
Name_suffix = 'Gamma'

############## GeomSvc #################
geometry_option = "CRD_o1_v01/CRD_o1_v01_TPC.xml"

if not os.getenv("DETCRDROOT"):
    print("Can't find the geometry. Please setup envvar DETCRDROOT." )
    sys.exit(-1)

geometry_path = os.path.join(os.getenv("DETCRDROOT"), "compact", geometry_option)
if not os.path.exists(geometry_path):
    print("Can't find the compact geometry file: %s"%geometry_path)
    sys.exit(-1)

from Configurables import GeomSvc
geomsvc = GeomSvc("GeomSvc")
geomsvc.compact = geometry_path
#######################################

########### k4DataSvc ####################
from Configurables import k4DataSvc
podioevent = k4DataSvc("EventDataSvc")
podioevent.inputs = [
"CRD_TPC_Gamma.root"
]
##########################################


########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
inp.collections = ["EcalBarrelCollection", "EcalBarrelContributionCollection", "HcalBarrelCollection", "HcalBarrelContributionCollection", "MCParticle", "MarlinTrkTracks"]
##########################################

########## Digitalization ################

##ECAL##
from Configurables import CRDEcalDigiAlg
EcalDigi = CRDEcalDigiAlg("CRDEcalDigiAlg")
EcalDigi.ReadOutName = "EcalBarrelCollection"
EcalDigi.SimCaloHitCollection = "EcalBarrelCollection"
EcalDigi.CaloHitCollection = "ECALBarrel"
EcalDigi.CaloAssociationCollection = "RecoCaloAssociation_ECALBarrel"
EcalDigi.SkipEvt = Nskip
EcalDigi.Seed = 2079
#Digitalization parameters
EcalDigi.CalibrECAL = 1
EcalDigi.AttenuationLength = 7e10
EcalDigi.TimeResolution = 0.5        #unit: ns
EcalDigi.EnergyThreshold = 0.0001   #0.1 MeV
EcalDigi.ChargeThresholdFrac = 0.05
EcalDigi.Debug=1
EcalDigi.OutFileName = "testTree_"+Name_suffix+".root"
#########################################

##HCAL##
from Configurables import G2CDArborAlg
caloDigi = G2CDArborAlg("G2CDArborAlg")
caloDigi.ReadLCIO = False
#caloDigi.CalibrECAL = [48.16, 96.32]
caloDigi.CalibrECAL = [46.538, 93.0769]
caloDigi.ECALCollections = []
caloDigi.ECALReadOutNames= []
caloDigi.DigiECALCollection = []
caloDigi.HCALCollections = ["HcalBarrelCollection"]
caloDigi.HCALReadOutNames= ["HcalBarrelCollection"]
caloDigi.DigiHCALCollection = ["HCALBarrel"]
caloDigi.EventReportEvery = 100

##############################################################################
# POD I/O
##############################################################################
from Configurables import PodioOutput
out = PodioOutput("outputalg")
out.filename = "Digi_Gamma.root"
out.outputCommands = ["keep *"]


########################################

from Configurables import ApplicationMgr
ApplicationMgr( 
    TopAlg=[inp, EcalDigi,caloDigi, out],
    EvtSel="NONE",
    EvtMax=Nevt,
    ExtSvc=[podioevent, geomsvc],
    #OutputLevel=DEBUG
)

