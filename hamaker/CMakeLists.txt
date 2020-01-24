# directory for config.h
include_directories(${CMAKE_BINARY_DIR})

# build libhamaker.so
add_library(hamaker SHARED
  HamakerElementResponse.cc
  HamakerCoeff.cc)


# install libhamaker.so
install(TARGETS hamaker DESTINATION lib)

# install coefficients
message("install hamaker in: " ${CMAKE_INSTALL_DATA_DIR})
install(FILES "${CMAKE_SOURCE_DIR}/hamaker/HamakerHBACoeff.h5" DESTINATION ${CMAKE_INSTALL_DATA_DIR})
install(FILES "${CMAKE_SOURCE_DIR}/hamaker/HamakerLBACoeff.h5" DESTINATION ${CMAKE_INSTALL_DATA_DIR})