from Gaudi.Configuration import *
Nskip = 0
Nevt = 10
Name_suffix = 'PseudoJet'

############## GeomSvc #################
geometry_option = "CRD_o1_v01/CRD_o1_v01_HCAL.xml"

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
#"CRD_HCal_Gamma10GeV.root"
"PseudoJet_3Gam2Pi/Sim/Simu_PseudoJet_3Gam2Pi_000.root"
]
##########################################


########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
inp.collections = ["EcalBarrelCollection", "EcalBarrelContributionCollection", "HcalBarrelCollection", "HcalBarrelContributionCollection", "MCParticle", "MarlinTrkTracks", "MarlinTrkAssociation"]
##########################################

########## Digitalization ################

##ECAL##
from Configurables import CRDEcalDigiAlg
EcalDigi = CRDEcalDigiAlg("CRDEcalDigiAlg")
EcalDigi.ReadOutName = "EcalBarrelCollection"
EcalDigi.SimCaloHitCollection = "EcalBarrelCollection"
EcalDigi.CaloHitCollection = "ECALBarrel"
EcalDigi.CaloAssociationCollection = "ECALBarrelAssoCol"
EcalDigi.CaloMCPAssociationCollection = "ECALBarrelParticleAssoCol"
EcalDigi.SkipEvt = Nskip
EcalDigi.Seed = 2079
#Digitalization parameters
EcalDigi.CalibrECAL = 1
EcalDigi.AttenuationLength = 7e10
EcalDigi.TimeResolution = 0.5        #unit: ns
EcalDigi.EnergyThreshold = 0.0001   #0.1 MeV
EcalDigi.ChargeThresholdFrac = 0.05
EcalDigi.Debug=1
EcalDigi.OutFileName = "DigiEcal_"+Name_suffix+".root"
#########################################

##HCAL##
from Configurables import CRDHcalDigiAlg
HcalDigi = CRDHcalDigiAlg("CRDHcalDigiAlg")
HcalDigi.ReadOutName = "HcalBarrelCollection"
HcalDigi.SimCaloHitCollection = "HcalBarrelCollection"
HcalDigi.CaloHitCollection = "HCALBarrel"
HcalDigi.CaloAssociationCollection = "HCALBarrelAssoCol"
HcalDigi.CaloMCPAssociationCollection = "HCALBarrelParticleAssoCol"
HcalDigi.SkipEvt = Nskip
HcalDigi.Seed = 2079
#Digitalization parameters
HcalDigi.MIPResponse = 0.0005  # 0.5 MeV / MIP
HcalDigi.MIPThreshold = 0.5    # Unit: MIP
HcalDigi.Debug=0
HcalDigi.OutFileName = "DigiHcal_"+Name_suffix+".root"

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
PandoraPlusPFAlg.ECalMCPAssociationName = ["ECALBarrelParticleAssoCol"]
PandoraPlusPFAlg.HCalCaloHitCollections = ["HCALBarrel"]
PandoraPlusPFAlg.HCalReadOutNames = ["HcalBarrelCollection"]
PandoraPlusPFAlg.HCalMCPAssociationName = ["HCALBarrelParticleAssoCol"]

#----Algorithms----

PandoraPlusPFAlg.AlgList = ["TruthClusteringAlg",
                            "TruthMatchingAlg",
                            "TruthClusterMergingAlg" ]
PandoraPlusPFAlg.AlgParNames = [ [""] ,
                                 [""] ,
                                 [""] ]
PandoraPlusPFAlg.AlgParTypes = [ [""] ,
                                 [""] ,
                                 [""] ]
PandoraPlusPFAlg.AlgParValues = [ [""] ,
                                  [""] ,
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
    TopAlg=[inp, EcalDigi, HcalDigi, PandoraPlusPFAlg],
    #TopAlg=[inp, EcalDigi,caloDigi],
    EvtSel="NONE",
    EvtMax=Nevt,
    ExtSvc=[podioevent, geomsvc],
    #OutputLevel=DEBUG
)

