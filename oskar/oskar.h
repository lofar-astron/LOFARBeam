#include "oskar_vector_types.h"
#include "oskar_helper.h"


void oskar_evaluate_dipole_pattern_double(
    const int num_points,
    const double* theta,
    const double* phi,
    const double freq_hz,
    const double dipole_length_m,
    std::complex<double>* pattern);
