#include "OSKARElementResponse.h"
#include "oskar.h"

#include <iostream>

namespace LOFAR {
namespace StationResponse{

void OSKARElementResponseDipole::response(
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

} // namespace StationResponse
} // namespace LOFAR
