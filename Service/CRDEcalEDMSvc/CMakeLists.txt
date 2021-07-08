# gaudi_subdir(CRDEcalEDMSvc v0r0)

gaudi_add_library(CRDEcalEDMSvcLib
                  SOURCES src/CRDCaloBar.cc
                          src/CRDCaloBarCluster.cc
                          src/CRDCaloBarShower.cc
                          src/CRDCaloBlock.cc
                          src/CRDCaloHit2DShower.cc
                          src/CRDCaloHit3DShower.cc
                          src/Track.cc
                          src/TrackFitInEcal.cpp
                          src/CRDArborNode.cc
                          src/CRDArborTree.cc
                  LINK DataHelperLib 
                       GearSvc
                       Gaudi::GaudiKernel 
                       ${ROOT_LIBRARIES}
                       ${CLHEP_LIBRARIES} 
                       ${GEAR_LIBRARIES}
                       ${LCIO_LIBRARIES}
                       EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)

gaudi_add_module(CRDEcalEDMSvc 
                 SOURCES src/CRDEcalEDMSvc.cpp
                 LINK CRDEcalEDMSvcLib 
                 Gaudi::GaudiKernel
                 ${GEAR_LIBRARIES}
)

install(TARGETS CRDEcalEDMSvcLib CRDEcalEDMSvc
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

