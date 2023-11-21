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
"Digi_Gamma.root"
]
##########################################


########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
inp.collections = ["EcalBarrelCollection", 
                   "EcalBarrelContributionCollection", 
                   "ECALBarrel",
                   "HcalBarrelCollection", 
                   "HcalBarrelContributionCollection", 
                   "HCALBarrel",
                   "MCParticle", 
                   "MarlinTrkTracks",
                   "MCRecoCaloParticleAssociationCollection" ]
##########################################

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
                            "EnergySplittingAlg",
                            "EnergyTimeMatchingAlg"]
PandoraPlusPFAlg.AlgParNames = [ ["Par1", "Par2"],
                                 ["Par1"],
                                 ["Eth_localMax", "Eth_MaxWithNeigh"],
                                 [""],
                                 ["th_Layers","th_peak"],
                                 ["ReadinLocalMaxName", "th_Nshowers"],
                                 ["axis_Angle", "relP_Angle", "relP_Dis"],
                                 [""],
                                 [""] ]
PandoraPlusPFAlg.AlgParTypes = [ ["double", "double"],
                                 ["double"],
                                 ["double", "double"],
                                 [""],
                                 ["int", "int"], 
                                 ["string", "int"], 
                                 ["double","double", "double"],
                                 [""],
                                 [""]  ]
PandoraPlusPFAlg.AlgParValues = [ ["1.", "3.14"],
                                  ["1."],
                                  ["0.005", "0."],
                                  [""],
                                  ["10","3"], 
                                  ["LeftLocalMax", "3"], 
                                  ["1.57", "1.57", "70"],  #Pi/2, Pi/2, 70.
                                  [""],
                                  [""]  ]


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
    TopAlg=[inp, PandoraPlusPFAlg ],
    #TopAlg=[inp, EcalDigi,caloDigi],
    EvtSel="NONE",
    EvtMax=Nevt,
    ExtSvc=[podioevent, geomsvc],
    #OutputLevel=DEBUG
)

