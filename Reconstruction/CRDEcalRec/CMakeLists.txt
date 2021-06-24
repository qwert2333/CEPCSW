# Modules
gaudi_add_module(CRDEcalRec
                 SOURCES src/PandoraPlusPFAlg.cpp
                         src/EcalClusterReconstruction.cpp
                         src/PandoraPlusDataCol.cpp
                         src/Algorithm/EnergySplittingAlg.cpp
                         src/Algorithm/EnergyTimeMatchingAlg.cpp
                         src/Algorithm/ConeClusteringAlg.cpp
                         src/Algorithm/CandidateMakingAlg.cpp
                         src/Algorithm/ClusterMergingAlg.cpp
                         src/Tools/PFOCreator.cpp
                 LINK k4FWCore::k4FWCore
                      GearSvc
                      DetInterface
                      CRDEcalEDMSvcLib
                      Gaudi::GaudiKernel
                      Gaudi::GaudiAlgLib 
                      ${CLHEP_LIBRARIES}
                      ${GEAR_LIBRARIES} 
                      ${GSL_LIBRARIES} 
                      ${LCIO_LIBRARIES}
                      ${ROOT_LIBRARIES}
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)
target_include_directories(CRDEcalRec PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS CRDEcalRec
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

