########################################

gaudi_add_module(GenAlgo
                 SOURCES src/IGenTool.cpp 
                         src/GenAlgo.cpp 
                         src/GenEvent.cpp 
                         src/GenPrinter.cpp
                         # ------- Readers -------
                         src/GenReader.cpp 
                         src/StdHepRdr.cpp 
                         src/LCAscHepRdr.cc
                         src/HepevtRdr.cpp
                         src/SLCIORdr.cpp
                         src/HepMCRdr.cpp
                         # ------- Signals -------
                         src/GtGunTool.cpp
                         src/GtPythiaTool.cpp
                         # ------- Beam Background -------
                         src/GtBeamBackgroundTool.cpp
                         src/BeamBackgroundFileParserV0.cpp
                         src/BeamBackgroundFileParserV1.cpp
                         src/BeamBackgroundFileParserV2.cpp
                         src/GuineaPigPairsFileParser.cpp
												 # ------- Data ----------
												 src/DataProcessor.cpp
                 LINK ${ROOT_LIBRARIES}
                      k4FWCore::k4FWCore 
                      Gaudi::GaudiAlgLib
                      Gaudi::GaudiKernel
                      ${HEPMC_LIBRARIES}
                      ${CLHEP_LIBRARIES}
                      ${LCIO_LIBRARIES}
                      ${PYTHIA8_LIBRARY}
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
                      ROOT::EG
)

target_include_directories(GenAlgo
  PUBLIC ${LCIO_INCLUDE_DIRS}
  PUBLIC ${HEPMC_INCLUDE_DIRS}
)

install(TARGETS GenAlgo
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
