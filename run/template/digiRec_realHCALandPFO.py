from Gaudi.Configuration import *
Nskip = 0
Nevt = NEVT

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
"SIMFILE"
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
EcalDigi.WriteNtuple = 0
EcalDigi.OutFileName = "ECALDIGIFILE"
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
HcalDigi.WriteNtuple = 0
HcalDigi.OutFileName = "HCALDIGIFILE"

######### Reconstruction ################
from Configurables import PandoraPlusPFAlg
PandoraPlusPFAlg = PandoraPlusPFAlg("PandoraPlusPFAlg")
##----Global parameters----
PandoraPlusPFAlg.Seed = 1024
PandoraPlusPFAlg.BField = 3.
PandoraPlusPFAlg.Debug = 0
PandoraPlusPFAlg.SkipEvt = Nskip
PandoraPlusPFAlg.WriteAna = 1
PandoraPlusPFAlg.AnaFileName = "RECFILE"
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

PandoraPlusPFAlg.AlgList = ["GlobalClusteringAlg",      #1
                            "LocalMaxFindingAlg",       #2
                            "TruthTrackMatchingAlg",    #7
                            "TruthPatternRecAlg",       #8
                            "TruthEnergySplittingAlg",  #10
                            "TruthMatchingAlg",         #12
                            "HcalClusteringAlg",        #13
                            "PFOCreatingAlg" ]          #14
PandoraPlusPFAlg.AlgParNames = [ ["InputECALBars","OutputECAL1DClusters","OutputECALHalfClusters"],#1
                                 ["OutputLocalMaxName"],#2
                                 ["ReadinLocalMaxName","OutputLongiClusName"],#7
                                 ["ReadinLocalMaxName", "OutputLongiClusName", "DoAxisMerging", "ReadinAxisName", "OutputMergedAxisName"],#8
                                 ["ReadinAxisName", "OutputClusName", "OutputTowerName"],#10
                                 ["ReadinHFClusterName", "ReadinTowerName", "OutputClusterName"],#12
                                 ["InputHCALHits", "OutputHCALClusters"],#13
                                 ["ReadinECALClusters","ReadinHCALClusters","OutputCombPFO"] ]#14
PandoraPlusPFAlg.AlgParTypes = [ ["string","string","string"],#1
                                 ["string"],#2
                                 ["string","string"],#7
                                 ["string","string","bool","string","string"],#8
                                 ["string","string","string"],#10
                                 ["string","string","string"],#12
                                 ["string","string"],#13
                                 ["string","string","string"] ]#14
PandoraPlusPFAlg.AlgParValues = [ ["BarCol","Cluster1DCol","HalfClusterCol"],#1
                                  ["AllLocalMax"],#2
                                  ["AllLocalMax","TruthTrackAxis"],#7
                                  ["AllLocalMax","TruthAxis","1","TruthTrackAxis","TruthMergedAxis"],#8
                                  ["TruthMergedAxis","TruthESCluster","TruthESTower"],#10
                                  ["TruthESCluster","TruthESTower","EcalCluster"],#12
                                  ["HCALBarrel","HCALCluster"],#13
                                  ["EcalCluster","HCALCluster","CombPFO"] ]#14


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

