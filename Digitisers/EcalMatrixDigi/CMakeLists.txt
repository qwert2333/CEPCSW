# Modules
gaudi_add_module(EcalMatrixDigi
                 SOURCES src/CRDEcalMatrixDigiAlg.cpp
                 LINK k4FWCore::k4FWCore
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
install(TARGETS EcalMatrixDigi
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

