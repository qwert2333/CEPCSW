# Modules
gaudi_add_library(CrystalCaloRecLib
                  SOURCES src/Objects/CaloUnit.cc
                          src/Objects/CaloHit.cc
                          src/Objects/Calo1DCluster.cc
                          src/Objects/Calo2DCluster.cc
                          src/Objects/Calo3DCluster.cc
                          src/Objects/CaloHalfCluster.cc
                          src/Objects/Track.cc
                          src/Objects/HoughSpace.cc
                          src/Objects/HoughObject.cc
                          src/Objects/PFObject.cc
                  LINK    k4FWCore::k4FWCore
                          ROOT::TMVA
                          GearSvc
                          CrystalEcalSvcLib
                          DetInterface
                          Gaudi::GaudiKernel
                          Gaudi::GaudiAlgLib
                          ${CLHEP_LIBRARIES}
                          ${GEAR_LIBRARIES}
                          ${GSL_LIBRARIES}
                          ${LCIO_LIBRARIES}
                          ${ROOT_LIBRARIES}
                          ${TMVA_LIBRARIES}
                          EDM4CEPC::edm4cepc EDM4CEPC::edm4cepcDict
                          EDM4HEP::edm4hep EDM4HEP::edm4hepDict
                          ${DD4hep_COMPONENT_LIBRARIES}
                          DetIdentifier

)

gaudi_add_module(CrystalCaloRec
                 SOURCES src/CyberPFAlg.cpp
                         src/CyberDataCol.cpp
                         src/Tools/MCParticleCreator.cpp
                         src/Tools/TrackCreator.cpp
                         src/Tools/CaloHitsCreator.cpp
                         src/Tools/OutputCreator.cpp
                         src/Tools/TrackFitInEcal.cpp
                         src/Algorithm/ExampleAlg.cpp
                         src/Algorithm/TrackExtrapolatingAlg.cpp
                         src/Algorithm/GlobalClusteringAlg.cpp
                         src/Algorithm/LocalMaxFindingAlg.cpp
                         src/Algorithm/TrackMatchingAlg.cpp
                         src/Algorithm/HoughClusteringAlg.cpp
                         src/Algorithm/ConeClustering2DAlg.cpp
                         src/Algorithm/AxisMergingAlg.cpp
                         src/Algorithm/EnergySplittingAlg.cpp
                         src/Algorithm/EnergyTimeMatchingAlg.cpp
                         src/Algorithm/HcalClusteringAlg.cpp
                         src/Algorithm/ConeClusteringAlg.cpp
                         src/Algorithm/PFOCreatingAlg.cpp
                         src/Algorithm/TrackClusterConnectingAlg.cpp
                         src/Algorithm/PFOReclusteringAlg.cpp
                         src/Algorithm/TruthTrackMatchingAlg.cpp
                         src/Algorithm/TruthMatchingAlg.cpp
                         src/Algorithm/TruthPatternRecAlg.cpp
                         src/Algorithm/TruthEnergySplittingAlg.cpp
                         src/Algorithm/TruthClusteringAlg.cpp
                         src/Algorithm/TruthClusterMergingAlg.cpp
                 LINK    CrystalCaloRecLib
)
target_include_directories(CrystalCaloRec PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS CrystalCaloRecLib CrystalCaloRec
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

