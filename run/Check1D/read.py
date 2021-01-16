from Gaudi.Configuration import *

############## GeomSvc #################
geometry_path = "/cefs/higgs/guofy/cepcsoft/CEPCSW_dev/CEPCSW/Detector/DetCRD/compact/ecalBlock.xml"
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
#"SimSample_2gam/CrystalBlock_2gam500MeV_3cm_m1.root",
"PhotonShower_30GeV.root",
]
##########################################

########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
#inp.collections = ["EcalBarrelCollection", "EcalEndcapRingCollection", "EcalEndcapsCollection", "MCParticle"]
inp.collections = ["SimCalorimeterCol", "MCParticle"]
##########################################

########## Digitalization ################
from Configurables import CRDEcalBlockDigiAlg
EcalDigi = CRDEcalBlockDigiAlg("CRDEcalBlockDigiAlg")
EcalDigi.SimCaloHitCollection = "SimCalorimeterCol"
EcalDigi.CaloHitCollection = "ECALBarrel"
EcalDigi.CaloAssociationCollection = "RecoCaloAssociation_ECALBarrel"
EcalDigi.Seed = 2093
#Digitalization parameters
EcalDigi.CalibrECAL = 1
EcalDigi.AttenuationLength = 7e10
EcalDigi.TimeResolution = 0.5
EcalDigi.EnergyThreshold = 0.0000   #200 keV
EcalDigi.ChargeThresholdFrac = 0.05
#Reconstruction parameters
EcalDigi.ScndMomentThreshold = 250.
EcalDigi.SeedEnergyThreshold = 0.01
EcalDigi.ShowerEnergyThreshold = 0.
EcalDigi.ClusterEnergyThreshold = 0.
EcalDigi.SeedWithNeighThreshold = 0.
EcalDigi.ShowerWithTotThreshold = 0.05
EcalDigi.ClusterWithTotThreshold = 0.05
EcalDigi.Debug=0
EcalDigi.OutFileName = "OutTree_1gam_ForFit.root"
#########################################


from Configurables import ApplicationMgr
ApplicationMgr( 
    TopAlg=[inp, EcalDigi],
    EvtSel="NONE",
    EvtMax=100,
    ExtSvc=[podioevent, geomsvc],
    #OutputLevel=DEBUG
)

