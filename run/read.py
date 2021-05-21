from Gaudi.Configuration import *

############## GeomSvc #################
geometry_path = "/cefs/higgs/guofy/cepcsoft/CEPCSW_v4/Detector/DetCRD/compact/ecalBarrel.xml"
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
"PhotonShower_10GeV.root"
]
##########################################

########### CRDEcalEdmSvc ################
from Configurables import CRDEcalEDMSvc
Ecaldatasvc = CRDEcalEDMSvc("CRDEcalEDMSvc")


########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
#inp.collections = ["EcalBarrelCollection", "EcalEndcapRingCollection", "EcalEndcapsCollection", "MCParticle"]
inp.collections = ["SimCalorimeterCol", "MCParticle"]
##########################################

########## Digitalization ################
from Configurables import CRDEcalDigiAlg
EcalDigi = CRDEcalDigiAlg("CRDEcalDigiAlg")
EcalDigi.ReadOutName = "CaloHitsCollection"
EcalDigi.SimCaloHitCollection = "SimCalorimeterCol"
EcalDigi.CaloHitCollection = "ECALBarrel"
EcalDigi.CaloAssociationCollection = "RecoCaloAssociation_ECALBarrel"
#EcalDigi.SkipEvt = 1
EcalDigi.Seed = 2079
#Digitalization parameters
EcalDigi.CalibrECAL = 1
EcalDigi.AttenuationLength = 7e10
EcalDigi.TimeResolution = 0.5
EcalDigi.EnergyThreshold = 0.0001   #0.1 MeV
EcalDigi.ChargeThresholdFrac = 0.05
EcalDigi.Debug=0
EcalDigi.OutFileName = "OutTree_testrec.root"
#########################################


######### Reconstruction ################
from Configurables import PandoraPlusPFAlg
PandoraPlusPFAlg = PandoraPlusPFAlg("PandoraPlusPFAlg")
PandoraPlusPFAlg.Seed = 1024
#PandoraPlusPFAlg.SkipEvt = 1
PandoraPlusPFAlg.Debug = 0
PandoraPlusPFAlg.OutFileName = "testMu.root"
########################################

########################################

from Configurables import ApplicationMgr
ApplicationMgr( 
    TopAlg=[inp, EcalDigi, PandoraPlusPFAlg],
    EvtSel="NONE",
    EvtMax=10,
    ExtSvc=[podioevent, geomsvc, Ecaldatasvc],
    #OutputLevel=DEBUG
)

