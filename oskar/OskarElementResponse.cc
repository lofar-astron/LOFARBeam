#include "config.h"

#include "OskarElementResponse.h"
#include "OskarSphericalWaveCoeff.h"

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

void oskar_evaluate_spherical_wave_sum_double(
    const int num_points,
    const double* theta,
    const double* phi_x,
    const double* phi_y,
    const int l_max,
    const std::complex<double>* alpha,
    std::complex<double>* pattern);

void OskarElementResponseSphericalWave::element_response(
    double freq,
    double theta,
    double phi,
    std::complex<double> (&response)[2][2]) const
{
    int l_max = 1; // TODO
    std::complex<double>* response_ptr = (std::complex<double> *) response;
    std::complex<double>* alpha_ptr; // TODO

    double phi_x = phi;
    double phi_y = phi + M_PI/2;
    oskar_evaluate_spherical_wave_sum_double(1, &theta, &phi_x, &phi_y, l_max, alpha_ptr, response_ptr);
}

std::string OskarElementResponseSphericalWave::get_path(
    const char* filename) const
{
    std::stringstream ss;
    ss << LOFARBEAM_DATA_DIR << "/";
    ss << filename;
    return ss.str();
}

OskarElementResponseSphericalWave::OskarElementResponseSphericalWave()
{
    std::string path = get_path("oskar.h5");
    m_coeffs.reset(new OskarSphericalWaveCoefficients(path));
}