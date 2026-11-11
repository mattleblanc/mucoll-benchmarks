cmake_minimum_required(VERSION 3.16)

project(CalibHits LANGUAGES CXX)

find_package(DD4hep REQUIRED COMPONENTS DDCore DDG4)
find_package(Geant4 REQUIRED)

dd4hep_set_compiler_flags()
dd4hep_configure_output()

dd4hep_add_plugin(CalibHits
  SOURCES src/CalibHitCalorimeterSDAction.cpp src/EnergyBalanceCheckAction.cpp
  USES    DD4hep::DDCore DD4hep::DDG4 ${Geant4_LIBRARIES}
)

install(TARGETS CalibHits LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
