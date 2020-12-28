from Gaudi.Configuration import *

############## GeomSvc #################
geometry_path = "/cefs/higgs/guofy/cepcsoft/CEPCSW/Detector/DetCRD/compact/ecalMatrix_cross.xml"
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
"run/gen_test.root"
#"CRDEcalMatrix_Mu.root"
]
##########################################

########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
#inp.collections = ["EcalBarrelCollection", "EcalEndcapRingCollection", "EcalEndcapsCollection", "MCParticle"]
inp.collections = ["SimCalorimeterCol", "MCParticle"]
##########################################

########## Digitalization ################
from Configurables import CRDEcalMatrixDigiAlg
EcalDigi = CRDEcalMatrixDigiAlg("CRDEcalMatrixDigiAlg")
EcalDigi.SimCaloHitCollection = "SimCalorimeterCol"
EcalDigi.CaloHitCollection = "ECALBarrel"
EcalDigi.CaloAssociationCollection = "RecoCaloAssociation_ECALBarrel"
EcalDigi.CalibrECAL = 1
EcalDigi.Seed = 1013
EcalDigi.AttenuationLength = 7e10
EcalDigi.TimeResolution = 0.
EcalDigi.EnergyThreshold = 0.003
EcalDigi.PositionThreshold = 1e5
EcalDigi.ChargeThresholdFrac = 0.05
EcalDigi.Debug=0
EcalDigi.OutFileName = "test.root"
#########################################


#from Configurables import PodioOutput
#out = PodioOutput("outputalg")
#out.filename = "CRDEcalRecv01_2gam60.root"
#out.outputCommands = ["keep *"]


from Configurables import ApplicationMgr
ApplicationMgr( 
    TopAlg=[inp, EcalDigi],
    EvtSel="NONE",
    EvtMax=2,
    ExtSvc=[podioevent, geomsvc],
    #OutputLevel=DEBUG
)

