from Gaudi.Configuration import *

############## GeomSvc #################
#geometry_path = "/cefs/higgs/guofy/cepcsoft/CEPCSW_v10/Detector/DetCRD/compact/ecalBarrel.xml"
geometry_path = os.getenv("WorkDIR")+"/Detector/DetCRD/compact/CRD_o1_v01/CRD_o1_v01.xml"
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
#"/scratchfs/atlas/guofy/CEPCEcalSim/SingleParticle/NoBfield/SinglePhoton_10GeV.root"
"../simdir/yyy_sim.root"
]
##########################################

########### CRDEcalEdmSvc ################
from Configurables import CRDEcalSvc
Ecaldatasvc = CRDEcalSvc("CRDEcalSvc")


########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
#inp.collections = ["EcalBarrelCollection", "EcalEndcapRingCollection", "EcalEndcapsCollection", "MCParticle"]
inp.collections = ["EcalBarrelCollection", "MCParticle", "MarlinTrkTracks"]
##########################################

########## Digitalization ################
from Configurables import CRDEcalDigiAlg
EcalDigi = CRDEcalDigiAlg("CRDEcalDigiAlg")
EcalDigi.ReadOutName = "EcalBarrelCollection"
EcalDigi.SimCaloHitCollection = "EcalBarrelCollection"
EcalDigi.CaloHitCollection = "ECALBarrel"
EcalDigi.CaloAssociationCollection = "RecoCaloAssociation_ECALBarrel"
#EcalDigi.SkipEvt = 21
EcalDigi.Seed = 2079
#Digitalization parameters
EcalDigi.CalibrECAL = 1
EcalDigi.AttenuationLength = 7e10
EcalDigi.TimeResolution = 0.5        #unit: ns
EcalDigi.EnergyThreshold = 0.0001   #0.1 MeV
EcalDigi.ChargeThresholdFrac = 0.05
EcalDigi.Debug=0
EcalDigi.OutFileName = "../digidir/yyy_digi.root"
#########################################


######### Reconstruction ################
from Configurables import PandoraPlusPFAlg
PandoraPlusPFAlg = PandoraPlusPFAlg("PandoraPlusPFAlg")
PandoraPlusPFAlg.Seed = 1024
#PandoraPlusPFAlg.SkipEvt = 21
PandoraPlusPFAlg.Debug = 0
PandoraPlusPFAlg.OutFileName = "../recdir/yyy_rec.root"
########################################

########################################

from Configurables import ApplicationMgr
ApplicationMgr( 
    TopAlg=[inp, EcalDigi, PandoraPlusPFAlg],
    EvtSel="NONE",
    EvtMax=yyy_events,
    ExtSvc=[podioevent, geomsvc, Ecaldatasvc],
    #OutputLevel=DEBUG
)

