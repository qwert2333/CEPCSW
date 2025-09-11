find_package(xgboost)

gaudi_add_header_only_library(FinalPIDSvc)

gaudi_add_library(FinalPIDSvcLib
                SOURCES src/WorkingPoint.cc
                        src/MuonExtrapolator.cc
                LINK GearSvc
                    Gaudi::GaudiKernel
                    DataHelperLib
                    DetInterface
                    DetSegmentation
                    TrackSystemSvcLib
                    ${DD4hep_COMPONENT_LIBRARIES}
                    ${ROOT_LIBRARIES}
                    ${GSL_LIBRARIES}
                    ${GEAR_LIBRARIES}
                    ${LCIO_LIBRARIES}
                    EDM4CEPC::edm4cepc EDM4CEPC::edm4cepcDict
                    EDM4HEP::edm4hep EDM4HEP::edm4hepDict
                    k4FWCore::k4FWCore
                    xgboost::xgboost
)
gaudi_add_module(FinalPIDSvcPlugins
                SOURCES src/FinalPIDSvc.cpp
                LINK FinalPIDSvc
                    FinalPIDSvcLib
                    GearSvc
                    Gaudi::GaudiKernel
                    DataHelperLib
                    DetInterface
                    DetSegmentation
                    TrackSystemSvcLib
                    ${DD4hep_COMPONENT_LIBRARIES}
                    ${ROOT_LIBRARIES}
                    ${GSL_LIBRARIES}
                    ${GEAR_LIBRARIES}
                    ${LCIO_LIBRARIES}
                    EDM4CEPC::edm4cepc EDM4CEPC::edm4cepcDict
                    EDM4HEP::edm4hep EDM4HEP::edm4hepDict
                    k4FWCore::k4FWCore
                    xgboost::xgboost
)

install(TARGETS FinalPIDSvc FinalPIDSvcPlugins
    EXPORT CEPCSWTargets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
    COMPONENT dev)