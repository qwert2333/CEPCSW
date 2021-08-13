# gaudi_subdir(CRDEcalSvc v0r0)

gaudi_add_header_only_library(CRDEcalSvcLib)

gaudi_add_module(CRDEcalSvc 
                 SOURCES src/CRDEcalSvc.cpp
                 LINK Gaudi::GaudiKernel
                      CRDEcalSvcLib
                      DataHelperLib
                      CRDEcalRecLib
                      GearSvc
                      Gaudi::GaudiKernel
                      ${ROOT_LIBRARIES}
                      ${CLHEP_LIBRARIES}
                      ${GEAR_LIBRARIES}
                      ${LCIO_LIBRARIES}
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)

#target_include_directories(CRDEcalSvc PUBLIC
#  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
#  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS  CRDEcalSvcLib CRDEcalSvc
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

