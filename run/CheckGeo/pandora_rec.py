from Gaudi.Configuration import *

############## GeomSvc #################
#geometry_path = "/cefs/higgs/guofy/cepcsoft/CEPCSW_dev/CEPCSW/Detector/DetCRD/compact/CRD_ECAL/CepC_v4-onlyECAL.xml"
geometry_path = "/cefs/higgs/guofy/cepcsoft/CEPCSW_dev/CEPCSW/Detector/DetCRD/compact/CRD_ECAL/CepC_v4.xml"
if not os.path.exists(geometry_path):
    print("Can't find the compact geometry file: %s"%geometry_path)
    sys.exit(-1)

from Configurables import GeomSvc
geomsvc = GeomSvc("GeoSvc")
geomsvc.compact = geometry_path
#######################################


##############GEAR Svc#################
from Configurables import GearSvc
gearSvc  = GearSvc("GearSvc")
gearSvc.GearXMLFile = "/cefs/higgs/guofy/cepcsoft/CEPCSW_dev/CEPCSW/Detector/DetCEPCv4/compact/FullDetGear.xml"


########### k4DataSvc ####################
from Configurables import k4DataSvc
podioevent = k4DataSvc("EventDataSvc")
podioevent.inputs = [
"CRDGlobalTest_eePFA.root"
]
##########################################

########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
inp.collections = [
	"VXDCollection",
	"SITCollection",
	"TPCCollection",
	"SETCollection",
	"EcalBarrelCollection", 
	"EcalEndcapRingCollection", 
	"EcalEndcapsCollection", 
	"HcalBarrelCollection",
	"MCParticle"]
#inp.collections = ["SimCalorimeterCol", "MCParticle"]
##########################################

########## Digitalization ################
from Configurables import CRDEcalDigiAlg
EcalDigi = CRDEcalDigiAlg("CRDEcalDigiAlg")
EcalDigi.SimCaloHitCollection = "EcalBarrelCollection"
EcalDigi.CaloHitCollection = "ECALBarrel"
EcalDigi.CaloAssociationCollection = "MCRecoCaloAssociationCollection"
EcalDigi.CalibrECAL = 1
EcalDigi.Seed = 1013
#Digitalization parameters
EcalDigi.AttenuationLength = 7e10
EcalDigi.TimeResolution = 0.3
EcalDigi.EnergyThreshold = 0.003
EcalDigi.ChargeThresholdFrac = 0.05
#Reconstruction parameters
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

##############################################################################
from Configurables import NTupleSvc
ntsvc = NTupleSvc("NTupleSvc")
ntsvc.Output = ["MyTuples DATAFILE='detsim_Pan_ana.root' OPT='NEW' TYP='ROOT'"]
##############################################################################

##############################################################################
# Pandora 
##############################################################################
from Configurables import PandoraPFAlg

pandoralg = PandoraPFAlg("PandoraPFAlg")
pandoralg.debug              = True
pandoralg.use_dd4hep_geo     = True
pandoralg.use_dd4hep_decoder = False
pandoralg.use_preshower      = False
pandoralg.WriteAna           = True
pandoralg.collections = [
        "MCParticle:MCParticle",
        "CalorimeterHit:ECALBarrel",
        "CalorimeterHit:ECALEndcap",
        "CalorimeterHit:ECALOther" ,
        "CalorimeterHit:HCALBarrel",
        "CalorimeterHit:HCALEndcap",
        "CalorimeterHit:HCALOther" ,
        "CalorimeterHit:MUON",
        "CalorimeterHit:LCAL",
        "CalorimeterHit:LHCAL",
        "CalorimeterHit:BCAL",
        "Vertex:KinkVertices",
        "Vertex:ProngVertices",
        "Vertex:SplitVertices",
        "Vertex:V0Vertices",
        "Track:MarlinTrkTracks",
        "MCRecoCaloAssociation:MCRecoCaloAssociationCollection"
        ]
pandoralg.WriteClusterCollection               = "PandoraClusters"
pandoralg.WriteReconstructedParticleCollection = "PandoraPFOs"
pandoralg.WriteVertexCollection                = "PandoraPFANewStartVertices"

pandoralg.PandoraSettingsDefault_xml = "/cefs/higgs/guofy/cepcsoft/CEPCSW_dev/CEPCSW/Reconstruction/PFA/Pandora/PandoraSettingsDefault.xml"
#### Do not chage the collection name, only add or remove ###############
pandoralg.TrackCollections      =  ["MarlinTrkTracks"]
pandoralg.ECalCaloHitCollections=  ["ECALBarrel", "ECALEndcap", "ECALOther"]
pandoralg.ECalReadOutNames      =  ["EcalBarrelCollection", "EcalEndcapsCollection", "ECALOther"]
pandoralg.HCalCaloHitCollections=  ["HCALBarrel", "HCALEndcap", "HCALOther"]
pandoralg.HCalReadOutNames      =  ["HcalBarrelCollection", "HcalEndcapsCollection", "HCALOther"]
pandoralg.LCalCaloHitCollections=  ["LCAL"]
pandoralg.LCalReadOutNames      =  ["LcalCollection"]
pandoralg.LHCalCaloHitCollections= ["LHCAL"]
pandoralg.LHCalReadOutNames      = ["LHcalCollection"]
pandoralg.MuonCaloHitCollections=  ["MUON"]
pandoralg.MuonCalReadOutNames    = ["MuoncalCollection"]
pandoralg.MCParticleCollections =  ["MCParticle"]
pandoralg.RelCaloHitCollections =  ["MCRecoCaloAssociationCollection"]
pandoralg.RelTrackCollections   =  ["MarlinTrkTracksMCTruthLink"]
pandoralg.KinkVertexCollections =  ["KinkVertices"]
pandoralg.ProngVertexCollections=  ["ProngVertices"]
pandoralg.SplitVertexCollections=  ["SplitVertices"]
pandoralg.V0VertexCollections   =  ["V0Vertices"]

pandoralg.ECalToMipCalibration  = 112 #1000MeV/8.918
pandoralg.HCalToMipCalibration  = 34.8
pandoralg.ECalMipThreshold      = 0.225# 8.918*0.225=2.00655
pandoralg.HCalMipThreshold      = 0.3
pandoralg.ECalToEMGeVCalibration= 1.# BGO, to be tuned
pandoralg.HCalToEMGeVCalibration= 1.
pandoralg.ECalToHadGeVCalibrationBarrel= 1.
pandoralg.ECalToHadGeVCalibrationEndCap= 1.
pandoralg.HCalToHadGeVCalibration= 1.
pandoralg.MuonToMipCalibration= 10.0
pandoralg.DigitalMuonHits= 0
pandoralg.MaxHCalHitHadronicEnergy   = 1.0
pandoralg.UseOldTrackStateCalculation= 0
pandoralg.AbsorberRadLengthECal= 0.08945 #BG0: 1/11.18 mm 
pandoralg.AbsorberIntLengthECal= 0.00448 #BG0: 1/223.2 mm 
pandoralg.AbsorberRadLengthHCal= 0.0569
pandoralg.AbsorberIntLengthHCal= 0.006
pandoralg.AbsorberRadLengthOther= 0.0569
pandoralg.AbsorberIntLengthOther= 0.006


# write PODIO file
from Configurables import PodioOutput
write = PodioOutput("write")
write.filename = "test.root"
write.outputCommands = ["keep *"]



from Configurables import ApplicationMgr
ApplicationMgr( 
    TopAlg=[inp, EcalDigi, pandoralg],
    EvtSel="NONE",
    EvtMax=1,
    ExtSvc=[podioevent, geomsvc, gearSvc], 
    HistogramPersistency = "ROOT",
    #OutputLevel=DEBUG
)

