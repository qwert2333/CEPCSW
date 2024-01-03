from Gaudi.Configuration import *
Nskip = 0
Nevt = 1
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

PandoraPlusPFAlg.AlgList = ["GlobalClusteringAlg",      #1
                            "LocalMaxFindingAlg",       #2
                            "TrackMatchingAlg",         #3
                            "HoughClusteringAlg",       #4
                            "ConeClustering2DAlg",      #5
                            "AxisMergingAlg",           #6
                            "TruthTrackMatchingAlg",    #7
                            "TruthPatternRecAlg",       #8
                            "EnergySplittingAlg",       #9
                            "TruthEnergySplittingAlg",  #10
                            "EnergyTimeMatchingAlg",    #11
                            "TruthMatchingAlg",         #12
                            "HcalClusteringAlg",        #13
                            "PFOCreatingAlg",           #14
                            "TruthClusteringAlg",       #15
                            "TruthClusterMergingAlg" ]  #16  
PandoraPlusPFAlg.AlgParNames = [ ["InputECALBars","OutputECAL1DClusters","OutputECALHalfClusters"],#1
                                 ["OutputLocalMaxName"],#2
                                 ["ReadinLocalMaxName","OutputLongiClusName"],#3
                                 ["ReadinLocalMaxName","LeftLocalMaxName","OutputLongiClusName"],#4
                                 ["ReadinLocalMaxName", "OutputLongiClusName"],#5
                                 ["OutputAxisName"],#6
                                 ["ReadinLocalMaxName","OutputLongiClusName"],#7
                                 ["ReadinLocalMaxName", "OutputLongiClusName", "DoAxisMerging", "ReadinAxisName", "OutputMergedAxisName"],#8
                                 ["ReadinAxisName", "OutputClusName", "OutputTowerName"],#9
                                 ["ReadinAxisName", "OutputClusName", "OutputTowerName"],#10
                                 ["ReadinHFClusterName", "ReadinTowerName","OutputClusterName"],#11
                                 ["ReadinHFClusterName", "ReadinTowerName", "OutputClusterName"],#12
                                 ["InputHCALHits", "OutputHCALClusters"],#13
                                 ["ReadinECALClusters","ReadinHCALClusters","OutputCombPFO"],#14
                                 ["DoECALClustering","DoHCALClustering","InputHCALHits","OutputHCALClusters"],#15
                                 ["DoECALMerge","DoHCALMerge","DoECALHCALConnection","ReadinECALClusters","ReadinHCALClusters","OutputECALCluster","OutputHCALCluster","OutputCombPFO"] ]#16
PandoraPlusPFAlg.AlgParTypes = [ ["string","string","string"],#1
                                 ["string"],#2
                                 ["string","string"],#3
                                 ["string","string","string"],#4
                                 ["string","string"],#5
                                 ["string"],#6
                                 ["string","string"],#7
                                 ["string","string","bool","string","string"],#8
                                 ["string","string","string"],#9
                                 ["string","string","string"],#10
                                 ["string","string","string"],#11
                                 ["string","string","string"],#12
                                 ["string","string"],#13
                                 ["string","string","string"],#14
                                 ["bool","bool","string","string"],#15
                                 ["bool","bool","bool","string","string","string","string","string"] ]#16
PandoraPlusPFAlg.AlgParValues = [ ["BarCol","Cluster1DCol","HalfClusterCol"],#1
                                  ["AllLocalMax"],#2
                                  ["AllLocalMax","TrackAxis"],#3
                                  ["AllLocalMax","LeftLocalMax","HoughAxis"],#4
                                  ["LeftLocalMax","ConeAxis"],#5
                                  ["MergedAxis"],#6
                                  ["AllLocalMax","TruthTrackAxis"],#7
                                  ["AllLocalMax","TruthAxis","1","TruthTrackAxis","TruthMergedAxis"],#8
                                  ["MergedAxis","ESHalfCluster","ESTower"],#9
                                  ["TruthMergedAxis","TruthESCluster","TruthESTower"],#10
                                  ["ESHalfCluster","ESTower","EcalCluster"],#11
                                  ["TruthESCluster","TruthESTower","TruthEcalCluster"],#12
                                  ["HCALBarrel","HCALCluster"],#13
                                  ["EcalCluster","HCALCluster","outputPFO"],#14
                                  ["0","1","HCALBarrel","TruthHcalCluster"],#15
                                  ["1","1","1","TruthEcalCluster","TruthHcalCluster","TruthMergedEcalCluster","TruthMergedHcalCluster","TruthCombPFO"]  ]#16


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

