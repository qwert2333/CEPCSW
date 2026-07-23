import os, sys
from Gaudi.Configuration import *

############## GeomSvc #################
geometry_option = "TDR_o1_v01/TDR_o1_v01.xml"

if not os.getenv("DETCRDROOT"):
    print("Can't find the geometry. Please setup envvar DETCRDROOT." )
    sys.exit(-1)

geometry_path = os.path.join(os.getenv("DETCRDROOT"), "compact", geometry_option)
if not os.path.exists(geometry_path):
    print("Can't find the compact geometry file: %s"%geometry_path)
    sys.exit(-1)

from Configurables import DetGeomSvc
geomsvc = DetGeomSvc("GeomSvc")
geomsvc.compact = geometry_path
#######################################

########### k4DataSvc ####################
from Configurables import k4DataSvc
podioevent = k4DataSvc("EventDataSvc", input="Tracking_TDR_o1_v01.root")
##########################################

########## CEPCSWData ################# 
cepcswdatatop ="/cvmfs/cepcsw.ihep.ac.cn/prototype/releases/data/latest"
#######################################

########## CrystalEcalEnergyCorrectionSvc ########
from Configurables import CrystalEcalEnergyCorrectionSvc
crystalecalcorr = CrystalEcalEnergyCorrectionSvc("CrystalEcalEnergyCorrectionSvc")
crystalecalcorr.CorrectionFile = os.path.join(cepcswdatatop, "CEPCSWData/offline-data/Service/CrystalEcalSvc/data/CrackRegionEnergyCorrection.root")
##################################################

########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
inp.collections = [ 
                    "ECALBarrel",
                    "ECALBarrelParticleAssoCol",
#                    "ECALEndcaps",
#                    "ECALEndcapsParticleAssoCol",
                    "HCALBarrel",
                    "HCALBarrelParticleAssoCol",
#                    "HCALEndcaps",
#                    "HCALEndcapsParticleAssoCol",
                    "MCParticle", 
                    "CompleteTracks", 
                    "CompleteTracksParticleAssociation",
                    "RecTofCollection",
                    "DndxTracks"]
##########################################

######### Reconstruction ################
from Configurables import CyberPFAlg
CyberPFAlg = CyberPFAlg("CyberPFAlg")
##----Global parameters----
CyberPFAlg.Seed = 1024
CyberPFAlg.BField = 3.
CyberPFAlg.Debug = 0
CyberPFAlg.SkipEvt = 0
CyberPFAlg.WriteAna = 1
CyberPFAlg.AnaFileName = "RecAnaTuple_TDR_o1_v01.root"
CyberPFAlg.UseMCPTrack = 0
CyberPFAlg.UseTruthMatchTrack = 0
CyberPFAlg.DoCleanTrack = 1
CyberPFAlg.TrackIDFile = "/cvmfs/cepcsw.ihep.ac.cn/prototype/releases/data/latest/CEPCSWData/offline-data/Reconstruction/CyberPFA_trackID/TrkID_BDT_BDTG.weights.xml"
CyberPFAlg.TrackIDMethod = "BDTG"
CyberPFAlg.EcalChargedCalib = 1.26
CyberPFAlg.HcalChargedCalib = 4.0
CyberPFAlg.EcalNeutralCalib = 1.0
CyberPFAlg.HcalNeutralCalib = 4.0
##----Readin collections----
CyberPFAlg.MCParticleCollection = "MCParticle"
CyberPFAlg.TrackCollections = ["CompleteTracks"]
CyberPFAlg.MCRecoTrackParticleAssociationCollection = "CompleteTracksParticleAssociation"
#CyberPFAlg.ECalCaloHitCollections = ["ECALBarrel","ECALEndcaps"]
#CyberPFAlg.ECalReadOutNames = ["EcalBarrelCollection","EcalEndcapsCollection"]
#CyberPFAlg.ECalMCPAssociationName = ["ECALBarrelParticleAssoCol", "ECALEndcapsParticleAssoCol"]
#CyberPFAlg.HCalCaloHitCollections = ["HCALBarrel", "HCALEndcaps"]
#CyberPFAlg.HCalReadOutNames = ["HcalBarrelCollection", "HcalEndcapsCollection"]
#CyberPFAlg.HCalMCPAssociationName = ["HCALBarrelParticleAssoCol", "HCALEndcapsParticleAssoCol"]
CyberPFAlg.ECalCaloHitCollections = ["ECALBarrel"]
CyberPFAlg.ECalReadOutNames = ["EcalBarrelCollection"]
CyberPFAlg.ECalMCPAssociationName = ["ECALBarrelParticleAssoCol"]
CyberPFAlg.HCalCaloHitCollections = ["HCALBarrel"]
CyberPFAlg.HCalReadOutNames = ["HcalBarrelCollection"]
CyberPFAlg.HCalMCPAssociationName = ["HCALBarrelParticleAssoCol"]

##--- Output collections ---
CyberPFAlg.OutputPFO = "outputPFO";
CyberPFAlg.RecoPFOCollection = "CyberPFO"

#----Algorithms----
CyberPFAlg.AlgList = ["GlobalClusteringAlg",      #1
                            "LocalMaxFindingAlg",       #2
                            "TrackMatchingAlg",         #3
                            "HoughClusteringAlg",       #4
                            "ConeClustering2DAlg",      #5
                            "AxisMergingAlg",           #6
                            "EnergySplittingAlg",       #9
                            "EnergyTimeMatchingAlg",    #11
                            "HcalClusteringAlg",        #12
                            "TruthClusteringAlg",       #15
                            "TrackClusterConnectingAlg",  #16
                            "PFOReclusteringAlg" ]  #17
CyberPFAlg.AlgParNames = [ ["InputECALBars","OutputECAL1DClusters","OutputECALHalfClusters"],#1
                                 ["OutputLocalMaxName"],#2
                                 ["ReadinLocalMaxName","OutputLongiClusName"],#3
                                 ["ReadinLocalMaxName","LeftLocalMaxName","OutputLongiClusName"],#4
                                 ["ReadinLocalMaxName", "OutputLongiClusName"], #5
                                 ["OutputAxisName"], #6
                                 ["ReadinAxisName", "OutputClusName", "OutputTowerName"],  #9
                                 ["ReadinHFClusterName", "ReadinTowerName","OutputClusterName"], #11
                                 ["InputHCALHits", "OutputHCALClusters"], #12
                                 ["DoECALClustering","DoHCALClustering","InputHCALHits","OutputHCALClusters"], #15
                                 ["ReadinECALClusterName", "ReadinHCALClusterName", "OutputCombPFO"],  #16
                                 ["ECALChargedCalib", "HCALChargedCalib", "ECALNeutralCalib", "HCALNeutralCalib"] ]#17
CyberPFAlg.AlgParTypes = [ ["string","string","string"],#1
                                 ["string"],#2
                                 ["string","string"],#3
                                 ["string","string","string"],#4
                                 ["string","string"], #5
                                 ["string"], #6
                                 ["string","string","string"],  #9
                                 ["string","string","string"], #11
                                 ["string", "string"], #12
                                 ["bool","bool","string","string"], #15
                                 ["string","string","string"],  #16
                                 ["double","double", "double","double"] ]#17
CyberPFAlg.AlgParValues = [ ["BarCol","Cluster1DCol","HalfClusterCol"],#1
                                  ["AllLocalMax"],#2
                                  ["AllLocalMax","TrackAxis"],#3
                                  ["AllLocalMax","LeftLocalMax","HoughAxis"],#4
                                  ["LeftLocalMax","ConeAxis"], #5
                                  ["MergedAxis"], #6
                                  ["MergedAxis","ESHalfCluster","ESTower"],  #9
                                  ["ESHalfCluster","ESTower","EcalCluster"], #11
                                  ["HCALBarrel", "SimpleHCALCluster"], #12
                                  ["0","1","HCALBarrel","HCALCluster"], #15
                                  ["EcalCluster", "SimpleHCALCluster", "outputPFO"],  #16
                                  ["1.26","4.", "1.", "4."]  ]#17


from Configurables import FinalPIDAlg
pid = FinalPIDAlg("FinalPIDAlg")
pid.OutputPFOName = "CyberPFOPID"

##############################################################################
# POD I/O
##############################################################################
from Configurables import PodioOutput
out = PodioOutput("outputalg")
out.filename = "Rec_TDR_o1_v01.root"
out.outputCommands = ["keep *"]


########################################

from Configurables import ApplicationMgr
ApplicationMgr( 
    TopAlg=[inp, CyberPFAlg, pid, out ],
    EvtSel="NONE",
    EvtMax=10,
    ExtSvc=[podioevent, geomsvc],
    #OutputLevel=DEBUG
)
