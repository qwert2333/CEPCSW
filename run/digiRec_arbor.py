from Gaudi.Configuration import *
Nskip = 0
Nevt = 1

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
"CRDSim_Pi20GeV_FullDet.root"
]
##########################################


########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
#inp.collections = ["EcalBarrelCollection", "EcalEndcapRingCollection", "EcalEndcapsCollection", "MCParticle"]
inp.collections = ["EcalBarrelCollection", "HcalBarrelCollection", "MCParticleG4", "MarlinTrkTracks"]
#inp.collections = ["EcalBarrelCollection", "MCParticleG4"]
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
EcalDigi.OutFileName = "testTree_nnHaa.root"
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
PandoraPlusPFAlg.AnaFileName = "testRec_nnHaa.root"
##----Readin collections----
PandoraPlusPFAlg.MCParticleCollection = "MCParticleG4"
PandoraPlusPFAlg.TrackCollections = ["MarlinTrkTracks"]
PandoraPlusPFAlg.ECalCaloHitCollections = ["ECALBarrel"]
PandoraPlusPFAlg.ECalReadOutNames = ["EcalBarrelCollection"]
PandoraPlusPFAlg.HCalCaloHitCollections = ["HCALBarrel"]
PandoraPlusPFAlg.HCalReadOutNames = ["HcalBarrelCollection"]

#----Algorithms----

'''
PandoraPlusPFAlg.AlgList = ["ExampleAlg",
                            "ConeClusteringAlgHCAL" ]
PandoraPlusPFAlg.AlgParNames = [ ["Par1", "Par2"],
                                 ["ReadinHit", "OutputCluster"] ]
PandoraPlusPFAlg.AlgParTypes = [ ["double", "double"],
                                 ["string", "string"] ]
PandoraPlusPFAlg.AlgParValues = [ ["1.", "3.14"],
                                  ["HCALBarrel", "HCALCluster"] ]

'''
PandoraPlusPFAlg.AlgList = ["ExampleAlg", 
                            "GlobalClusteringAlg", 
                            "LocalMaxFindingAlg",
                            "HoughClusteringAlg",
                            "ConeClustering2DAlg",
                            "EnergySplittingAlg",
                            "EnergyTimeMatchingAlg",
                            "ConeClusteringAlgHCAL"  ]
PandoraPlusPFAlg.AlgParNames = [ ["Par1", "Par2"], 
                                 ["Par1"], 
                                 ["Eth_localMax", "Eth_MaxWithNeigh"],
                                 ["th_Layers", "th_AxisE"] ,
                                 [""],
                                 [""], 
                                 [""],
                                 ["ReadinHit", "OutputCluster"] ]
PandoraPlusPFAlg.AlgParTypes = [ ["double", "double"],
                                 ["double"],
                                 ["double", "double"],
                                 ["double", "double"],
                                 [""],
                                 [""],
                                 [""], 
                                 ["string", "string"]  ]
PandoraPlusPFAlg.AlgParValues = [ ["1.", "3.14"], 
                                  ["1."], 
                                  ["0.005", "0."],
                                  ["10", "0.5"],
                                  [""],
                                  [""],
                                  [""],
                                  ["HCALBarrel", "HCALCluster"]   ]

########################################
from Configurables import MarlinArbor
ArborRec = MarlinArbor("MarlinArbor")
ArborRec.ECALCollections = [""]
ArborRec.HCALCollections = ["HCALBarrel"]
ArborRec.ECALReadOutNames = [""]
ArborRec.HCALReadOutNames = ["HcalBarrelCollection"]


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
    TopAlg=[inp, EcalDigi,caloDigi, ArborRec, out ],
    #TopAlg=[inp, EcalDigi, PandoraPlusPFAlg],
    EvtSel="NONE",
    EvtMax=Nevt,
    ExtSvc=[podioevent, geomsvc],
    #OutputLevel=DEBUG
)

