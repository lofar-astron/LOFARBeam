# build liboskar.so
add_library(oskar SHARED
    OSKARElementResponse.cc
    oskar_evaluate_dipole_pattern.cc
    oskar_evaluate_spherical_wave_sum.cc)

# install liboskar.so
install(TARGETS oskar DESTINATION lib)