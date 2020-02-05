# directory for config.h
include_directories(${CMAKE_BINARY_DIR})

# build liboskar.so
add_library(oskar SHARED
    OSKARElementResponse.cc
    OskarSphericalWaveCoeff.cc
    oskar_evaluate_dipole_pattern.cc
    oskar_evaluate_spherical_wave_sum.cc)

# install liboskar.so
install(TARGETS oskar DESTINATION lib)
