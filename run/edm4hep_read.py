#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc",
                 # input="test.root"
                 input="/cefs/higgs/zyang/cepcsoft/CEPCSW/yang/pion/simdir/sim_pi0_5GeV_theta50_phiall.root"
)

from Configurables import PodioInput
podioinput = PodioInput("PodioReader", collections=[
    "MCParticleG4",
    "EcalBarrelCollection"
    ])

from Configurables import Edm4hepReadAlg
alg = Edm4hepReadAlg("Edm4hepReadAlg")
#alg.HeaderCol.Path = "EventHeader"
alg.MCParticleCol.Path = "MCParticleG4"
alg.SimCalorimeterHitCol.Path = "EcalBarrelCollection"

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [podioinput, alg],
                EvtSel = 'NONE',
                EvtMax = 1000,
                ExtSvc = [dsvc]
                #OutputLevel=DEBUG
)
