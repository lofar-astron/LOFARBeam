cmake_minimum_required(VERSION 2.8.8)

project(lofarbeam)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures)
include_directories(${CASACORE_INCLUDE_DIR})

add_compile_options(-std=c++11 -Wall -DNDEBUG)

add_library(stationresponse SHARED
  AntennaField.cc
  AntennaFieldHBA.cc
  AntennaFieldLBA.cc
  AntennaModelHBA.cc
  AntennaModelLBA.cc
  DualDipoleAntenna.cc
  ElementResponse.cc
  ITRFDirection.cc
  Dir2ITRF.cc
  LofarMetaDataUtil.cc
  MathUtil.cc
  Station.cc
  TileAntenna.cc
  Types.cc)

add_library(elementresponse SHARED
  ElementResponse.cc)


target_link_libraries(stationresponse ${CASACORE_LIBRARIES})
target_link_libraries(elementresponse ${CASACORE_LIBRARIES})

install (TARGETS stationresponse DESTINATION lib)
install (TARGETS elementresponse DESTINATION lib)

install (FILES
  AntennaField.h
  Constants.h
  LofarMetaDataUtil.h
  Station.h
  Types.h
  ITRFDirection.h
  Dir2ITRF.h
DESTINATION include/StationResponse)

install (FILES
  ElementResponse.h
DESTINATION include/ElementResponse)
