cmake_minimum_required(VERSION 3.0.0)

project(lofarbeam)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures)
include_directories(${CASACORE_INCLUDE_DIR})

add_compile_options(-std=c++11 -Wall -DNDEBUG -Wl,--no-undefined)

add_library(stationresponse SHARED
  AntennaField.cc
  AntennaFieldHBA.cc
  AntennaFieldLBA.cc
  AntennaModelHBA.cc
  AntennaModelLBA.cc
  DualDipoleAntenna.cc
  ElementResponse.cc
  ITRFConverter.cc
  ITRFDirection.cc
  LofarMetaDataUtil.cc
  MathUtil.cc
  Station.cc
  TileAntenna.cc
  Types.cc)

add_library(elementresponse SHARED
  ElementResponse.cc)

set_target_properties(stationresponse PROPERTIES VERSION 3)
set_target_properties(elementresponse PROPERTIES VERSION 3)

# The elementresponse lib does not require casacore, so
# linking casacore only to stationresponse is enough.
target_link_libraries(stationresponse ${CASACORE_LIBRARIES})

install (TARGETS stationresponse DESTINATION lib)
install (TARGETS elementresponse DESTINATION lib)

install (FILES
  AntennaField.h
  Constants.h
  LofarMetaDataUtil.h
  Station.h
  Types.h
  ITRFConverter.h
  ITRFDirection.h
DESTINATION include/StationResponse)

install (FILES
  ElementResponse.h
DESTINATION include/ElementResponse)

find_package(Doxygen)

if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating developer API documentation with Doxygen" VERBATIM)
elseif(DOXYGEN_FOUND)

endif(DOXYGEN_FOUND)

find_package (Python COMPONENTS Development)
message("Compiling Python ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR} support. If this is not the intended target use 'DPYTHON_EXECUTABLE' to control.")
# bind against boost-python for compilation of the python API
# depending on major release version of the python libraries
if (${PYTHON_VERSION_MAJOR} VERSION_LESS "3")
  find_package (Boost REQUIRED COMPONENTS python numpy)
  find_package (Casacore REQUIRED COMPONENTS python)
  message(STATUS "Boost python library: ${Boost_PYTHON_LIBRARY}")
  message(STATUS "Boost numpy library: ${Boost_NUMPY_LIBRARY}")
  message(STATUS "Python library ${PYTHON_LIBRARIES}")
  set(BOOST_PY_FOUND TRUE)
else (${PYTHON_VERSION_MAJOR} VERSION_LESS "3")
  find_package (Boost REQUIRED COMPONENTS python3 numpy3)
  find_package (Casacore REQUIRED COMPONENTS python3)
  message(STATUS "Boost python library: ${Boost_PYTHON3_LIBRARY}")
  message(STATUS "Boost numpy library: ${Boost_NUMPY3_LIBRARY}")
  message(STATUS "Python library ${PYTHON_LIBRARIES}")
  set(BOOST_PY_FOUND TRUE)
endif (${PYTHON_VERSION_MAJOR} VERSION_LESS "3")


if(BOOST_PY_FOUND)
    add_library(_stationresponse MODULE pystationresponse.cc)
    set_target_properties(_stationresponse PROPERTIES PREFIX "")
    target_include_directories(_stationresponse PRIVATE ${PYTHON_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
    if (${PYTHON_VERSION_MAJOR} VERSION_LESS "3")
      target_link_libraries(_stationresponse stationresponse ${CASA_PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARY} ${Boost_NUMPY_LIBRARY} ${PYTHON_LIBRARIES})
    else (${PYTHON_VERSION_MAJOR} VERSION_LESS "3")
      target_link_libraries(_stationresponse stationresponse ${CASA_PYTHON3_LIBRARY} ${Boost_PYTHON3_LIBRARY} ${Boost_NUMPY3_LIBRARY} ${PYTHON_LIBRARIES})
    endif (${PYTHON_VERSION_MAJOR} VERSION_LESS "3")
    install(TARGETS _stationresponse DESTINATION ${PYTHON_INSTALL_DIR}/lofar/stationresponse)
    install(FILES __init__.py DESTINATION ${PYTHON_INSTALL_DIR}/lofar/stationresponse)
endif(BOOST_PY_FOUND)

