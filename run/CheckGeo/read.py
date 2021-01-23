from Gaudi.Configuration import *

############## GeomSvc #################
geometry_path = "/cefs/higgs/guofy/cepcsoft/CEPCSW_dev/CEPCSW/Detector/DetCRD/compact/ecalBarrel.xml"
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
"CRDEcalGeoTest_2gamPFA.root"
]
##########################################

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
EcalDigi.Seed = 1013
#Digitalization parameters
EcalDigi.CalibrECAL = 1
EcalDigi.AttenuationLength = 7e10
EcalDigi.TimeResolution = 0.5
EcalDigi.EnergyThreshold = 0.0002   #200 keV
EcalDigi.ChargeThresholdFrac = 0.05
#Reconstruction parameters
EcalDigi.ScndMomentThreshold = 0
EcalDigi.SeedEnergyThreshold = 0.05
EcalDigi.ShowerEnergyThreshold = 0
EcalDigi.ClusterEnergyThreshold = 0
EcalDigi.SeedWithNeighThreshold = 0.4
EcalDigi.SeedWithTotThreshold   = 0.15
EcalDigi.ShowerWithTotThreshold = 0.05
EcalDigi.ClusterWithTotThreshold = 0.05
EcalDigi.EnergyChi2Weight = 1
EcalDigi.TimeChi2Weight = 1
EcalDigi.Chi2Threshold = 0.
EcalDigi.Debug=0
EcalDigi.OutFileName = "OutTree_2gam20_case1.root"
#########################################


#from Configurables import PodioOutput
#out = PodioOutput("outputalg")
#out.filename = "CRDEcalRecv01_2gam60.root"
#out.outputCommands = ["keep *"]


from Configurables import ApplicationMgr
ApplicationMgr( 
    TopAlg=[inp, EcalDigi],
    EvtSel="NONE",
    EvtMax=3,
    ExtSvc=[podioevent, geomsvc],
    #OutputLevel=DEBUG
)

