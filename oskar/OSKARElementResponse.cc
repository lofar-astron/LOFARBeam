#include "OSKARElementResponse.h"

#include "oskar_vector_types.h"
#include "oskar_helper.h"

#include <iostream>

void oskar_evaluate_dipole_pattern_double(
    const int num_points,
    const double* theta,
    const double* phi,
    const double freq_hz,
    const double dipole_length_m,
    std::complex<double>* pattern);

void OSKARElementResponseDipole::element_response(
    double freq,
    double theta,
    double phi,
    std::complex<double> (&response)[2][2]) const
{
    double dipole_length_m = 1; // TODO
    std::complex<double>* response_ptr = (std::complex<double> *) response;

    double phi_x = phi;
    double phi_y = phi + M_PI/2;
    oskar_evaluate_dipole_pattern_double(1, &theta, &phi_x, freq, dipole_length_m, response_ptr);
    oskar_evaluate_dipole_pattern_double(1, &theta, &phi_y, freq, dipole_length_m, response_ptr + 2);
}
