from Gaudi.Configuration import *
Nskip = 0
Nevt = 10

############## GeomSvc #################
geometry_option = "CRD_o1_v01/CRD_o1_v01.xml"

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
"CRD_SimGamMu_FullDet.root"
#"SimSamples/CRD_SimMu_FullDet_200evt.root"
#"CRD_Sim2Gam_FullDet.root"
#"SimSamples/CRD_Sim2Gam_FullDet_200evt.root"
]
##########################################


########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
#inp.collections = ["EcalBarrelCollection", "EcalEndcapRingCollection", "EcalEndcapsCollection", "MCParticle"]
inp.collections = ["EcalBarrelCollection", "HcalBarrelCollection", "MCParticle", "MarlinTrkTracks"]
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
EcalDigi.OutFileName = "testTree.root"
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

######### Reconstruction ################
from Configurables import PandoraPlusPFAlg
PandoraPlusPFAlg = PandoraPlusPFAlg("PandoraPlusPFAlg")
##----Global parameters----
PandoraPlusPFAlg.Seed = 1024
PandoraPlusPFAlg.BField = 3.
PandoraPlusPFAlg.Debug = 0
PandoraPlusPFAlg.SkipEvt = Nskip
PandoraPlusPFAlg.WriteAna = 1
PandoraPlusPFAlg.AnaFileName = "testRec.root"
##----Readin collections----
PandoraPlusPFAlg.MCParticleCollection = "MCParticle"
PandoraPlusPFAlg.TrackCollections = ["MarlinTrkTracks"]
PandoraPlusPFAlg.ECalCaloHitCollections = ["ECALBarrel"]
PandoraPlusPFAlg.ECalReadOutNames = ["EcalBarrelCollection"]
PandoraPlusPFAlg.HCalCaloHitCollections = ["HCALBarrel"]
PandoraPlusPFAlg.HCalReadOutNames = ["HcalBarrelCollection"]

#----Algorithms----

PandoraPlusPFAlg.AlgList = ["ExampleAlg", 
                            "GlobalClusteringAlg", 
                            "LocalMaxFindingAlg"]
PandoraPlusPFAlg.AlgParNames = [ ["Par1", "Par2"], 
                                 ["Par1"], 
                                 ["Eth_localMax", "Eth_MaxWithNeigh"] ]
PandoraPlusPFAlg.AlgParTypes = [ ["double", "double"],
                                 ["double"],
                                 ["double", "double"] ]
PandoraPlusPFAlg.AlgParValues = [ ["1.", "3.14"], 
                                  ["1."], 
                                  ["0.005", "0."] ]

########################################

##############################################################################
# POD I/O
##############################################################################
from Configurables import PodioOutput
out = PodioOutput("outputalg")
out.filename = "Digi_2Gam.root"
out.outputCommands = ["keep *"]


########################################

from Configurables import ApplicationMgr
ApplicationMgr( 
    TopAlg=[inp, EcalDigi,caloDigi, PandoraPlusPFAlg ],
    #TopAlg=[inp, EcalDigi, PandoraPlusPFAlg],
    EvtSel="NONE",
    EvtMax=Nevt,
    ExtSvc=[podioevent, geomsvc],
    #OutputLevel=DEBUG
)

