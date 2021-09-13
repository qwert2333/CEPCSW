# Modules
gaudi_add_library(CRDEcalRecLib
                  SOURCES src/Objects/CRDArborNode.cc
                          src/Objects/CRDArborTree.cc
                          src/Objects/CRDCaloBar.cc
                          src/Objects/CRDCaloBarCluster.cc
                          src/Objects/CRDCaloBarShower.cc
                          src/Objects/CRDCaloBlock.cc
                          src/Objects/CRDCaloHit2DShower.cc
                          src/Objects/CRDCaloHit3DCluster.cc
                          src/Objects/Track.cc
                          src/Objects/TrackFitInEcal.cpp
                  LINK    k4FWCore::k4FWCore
                          ROOT::TMVA
                          GearSvc
                          DetInterface
                          Gaudi::GaudiKernel
                          Gaudi::GaudiAlgLib
                          ${CLHEP_LIBRARIES}
                          ${GEAR_LIBRARIES}
                          ${GSL_LIBRARIES}
                          ${LCIO_LIBRARIES}
                          ${ROOT_LIBRARIES}
                          ${TMVA_LIBRARIES}
                          EDM4HEP::edm4hep EDM4HEP::edm4hepDict

)

gaudi_add_module(CRDEcalRec
                 SOURCES src/PandoraPlusPFAlg.cpp
                         src/EcalClusterReconstruction.cpp
                         src/PandoraPlusDataCol.cpp
                         src/Algorithm/ArborClusteringAlg.cpp
                         src/Algorithm/ArborTreeMergingAlg.cpp
                         src/Algorithm/BasicClusterIDAlg.cpp
                         src/Algorithm/ShadowMakingAlg.cpp
                         src/Algorithm/ClusterMergingAlg.cpp
                         src/Algorithm/ConeClusteringAlg.cpp
                         src/Algorithm/EnergySplittingAlg.cpp
                         src/Algorithm/EnergyTimeMatchingAlg.cpp
                         src/Tools/PFOCreator.cpp
                 LINK    CRDEcalRecLib
                         CRDEcalSvcLib
)
target_include_directories(CRDEcalRec PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS CRDEcalRecLib CRDEcalRec
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

