# Modules
gaudi_add_library(CRDEcalRecLib
                  SOURCES src/Objects/CaloUnit.cc
                          src/Objects/Calo1DCluster.cc
                          src/Objects/Calo2DCluster.cc
                          src/Objects/Calo3DCluster.cc
                          src/Objects/CaloBarShower.cc
                          src/Objects/CaloBarCluster.cc
                          src/Objects/CaloBlock.cc
                          src/Objects/CaloTower.cc
                          src/Objects/TransShower.cc
                          src/Objects/CaloCluster.cc .
                          src/Objects/Track.cc
                          src/Objects/LongiCluster.cc
                  LINK    k4FWCore::k4FWCore
                          GearSvc
                          DetInterface
                          Gaudi::GaudiKernel
                          Gaudi::GaudiAlgLib
                          ${CLHEP_LIBRARIES}
                          ${GEAR_LIBRARIES}
                          ${GSL_LIBRARIES}
                          ${LCIO_LIBRARIES}
                          ${ROOT_LIBRARIES}
                          EDM4HEP::edm4hep EDM4HEP::edm4hepDict

)

gaudi_add_module(CRDEcalRec
                 SOURCES src/PandoraPlusPFAlg.cpp
                         src/PandoraPlusDataCol.cpp
                         src/Tools/MCParticleCreator.cpp
                         src/Tools/TrackCreator.cpp
                         src/Tools/CaloHitsCreator.cpp
                         src/Tools/OutputCreator.cpp
                         src/Tools/TrackFitInEcal.cpp
                         src/Algorithm/ExampleAlg.cpp
                         src/Algorithm/TrackExtrapolatingAlg.cpp
                         src/Algorithm/GlobalClusteringAlg.cpp
                         src/Algorithm/LocalMaxFindingAlg.cpp
                         src/Algorithm/HoughClusteringAlg.cpp
                         src/Algorithm/EnergySplittingAlg.cpp
                         src/Algorithm/EnergyTimeMatchingAlg.cpp
                         src/Algorithm/ConeClusteringAlg.cpp
                 LINK    CRDEcalRecLib
)
target_include_directories(CRDEcalRec PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS CRDEcalRecLib CRDEcalRec
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

