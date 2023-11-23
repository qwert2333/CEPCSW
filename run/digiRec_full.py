from Gaudi.Configuration import *
Nskip = 0
Nevt = 5
Name_suffix = 'PseudoJet_truthMatch'

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
"/cefs/higgs/guofy/CEPCSW_v214/run/PseudoJet_Gams/PseudoJet_Module6/Sim/Simu_PseudoJet_Module6_000.root"
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

######### Reconstruction ################
from Configurables import PandoraPlusPFAlg
PandoraPlusPFAlg = PandoraPlusPFAlg("PandoraPlusPFAlg")
##----Global parameters----
PandoraPlusPFAlg.Seed = 1024
PandoraPlusPFAlg.BField = 3.
PandoraPlusPFAlg.Debug = 0
PandoraPlusPFAlg.SkipEvt = Nskip
PandoraPlusPFAlg.WriteAna = 1
PandoraPlusPFAlg.AnaFileName = "testRec_"+Name_suffix+".root"
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
                            "LocalMaxFindingAlg",
                            "TrackMatchingAlg" ,
                            "HoughClusteringAlg",
                            "ConeClustering2DAlg",
                            "AxisMergingAlg",
                            "EnergySplittingAlg" ]
#                            "TruthMatchingAlg"]
PandoraPlusPFAlg.AlgParNames = [ ["Par1", "Par2"],
                                 ["Par1"],
                                 ["Eth_localMax", "Eth_MaxWithNeigh"],
                                 [""],
                                 ["th_Layers","th_peak"],
                                 ["ReadinLocalMaxName", "th_Nshowers"],
                                 ["axis_Angle", "relP_Angle", "relP_Dis"],
                                 [""] ]
#                                 [""] ]
PandoraPlusPFAlg.AlgParTypes = [ ["double", "double"],
                                 ["double"],
                                 ["double", "double"],
                                 [""],
                                 ["int", "int"], 
                                 ["string", "int"], 
                                 ["double","double", "double"],
                                 [""] ]
#                                 [""]  ]
PandoraPlusPFAlg.AlgParValues = [ ["1.", "3.14"],
                                  ["1."],
                                  ["0.005", "0."],
                                  [""],
                                  ["10","3"], 
                                  ["LeftLocalMax", "3"], 
                                  ["1.57", "1.57", "70"],  #Pi/2, Pi/2, 70.
                                  [""] ]
#                                  [""]  ]


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
    #TopAlg=[inp, EcalDigi,caloDigi],
    EvtSel="NONE",
    EvtMax=Nevt,
    ExtSvc=[podioevent, geomsvc],
    #OutputLevel=DEBUG
)

