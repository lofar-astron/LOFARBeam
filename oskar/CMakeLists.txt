add_library(oskar SHARED
    OSKARElementResponse.cc
    oskar_evaluate_dipole_pattern.cc)

install(TARGETS oskar DESTINATION lib)