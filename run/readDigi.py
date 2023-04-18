from Gaudi.Configuration import *

############## GeomSvc #################
geometry_option = "CepC_v4.xml"

if not os.getenv("DETCEPCV4ROOT"):
    print("Can't find the geometry. Please setup envvar DETCEPCV4ROOT." )
    sys.exit(-1)

geometry_path = os.path.join(os.getenv("DETCEPCV4ROOT"), "compact", geometry_option)
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
"SimSample/E240_nnHgg_50.root"
#"SimSample/GamPi_5GeV_2deg_1.root",
]
##########################################



########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
inp.collections = ["EcalBarrelCollectionMerged", "EcalEndcapRingCollection", "EcalEndcapsCollection", "HcalBarrelCollection", "HcalEndcapsCollection", "MCParticle"]
#inp.collections = ["MCParticle",
#                   "VXDTrackerHits", 
#                   "SITTrackerHits", 
#                   "SETTrackerHits", 
#                   "FTDTrackerHits", 
#                   "DCTrackerHits", 
#                   "EcalBarrelCollection",
#                   "EcalEndcapRingCollection",
#                   "EcalEndcapsCollection",
#                   "HcalBarrelCollection", 
#                   "HcalEndcapsCollection", 
#                   "MuonBarrelCollection", 
#                   "MuonEndcapsCollection" ]
##########################################

########## Digitalization ################
from Configurables import G2CDArborAlg
caloDigi = G2CDArborAlg("G2CDArborAlg")
caloDigi.ReadLCIO = False
#caloDigi.CalibrECAL = [48.16, 96.32]
caloDigi.CalibrECAL = [46.538, 93.0769]
caloDigi.ECALCollections = ["EcalBarrelCollectionMerged", "EcalEndcapsCollection"]
caloDigi.ECALReadOutNames= ["EcalBarrelCollection", "EcalEndcapsCollection"]
caloDigi.DigiECALCollection = ["ECALBarrel", "ECALEndcap"]
caloDigi.HCALCollections = ["HcalBarrelCollection", "HcalEndcapsCollection"]
caloDigi.HCALReadOutNames= ["HcalBarrelCollection", "HcalEndcapsCollection"]
caloDigi.DigiHCALCollection = ["HCALBarrel", "HCALEndcap"]
caloDigi.EventReportEvery = 1000


########## Read Info ###################
from Configurables import ReadDigiAlg
readin = ReadDigiAlg("ReadDigiAlg")
readin.OutFileName = "Digi_nnHgg.root"
#readin.OutputLevel = DEBUG
#########################################

# write PODIO file
from Configurables import PodioOutput
write = PodioOutput("write")
write.filename = "digi.root"
write.outputCommands = ["keep *"]

########################################

from Configurables import ApplicationMgr
ApplicationMgr( 
    TopAlg=[inp, caloDigi, readin],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[podioevent, geomsvc],
    #OutputLevel=DEBUG
)

