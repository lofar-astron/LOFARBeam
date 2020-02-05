#include "OSKARElementResponse.h"
#include "oskar.h"
#include "config.h"
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

OSKARElementResponseSphericalWave::OSKARElementResponseSphericalWave()
{
    std::string path = get_path("oskar.h5");
    m_coeffs.reset(new OskarSphericalWaveCoefficients(path));
}

void OSKARElementResponseSphericalWave::response(
    double freq,
    double theta,
    double phi,
    std::complex<double> (&response)[2][2]) const
{
    std::lock_guard<std::mutex> lock(m_mutex);

    // DEBUG
    // For now fix frequency
    freq = 140e6;

    if (!m_coeffs->read_frequency(freq)) {
        return;
    }

    int l_max = m_coeffs->get_l_max();
    int element = 0; // TODO
    std::complex<double>* response_ptr = (std::complex<double> *) response;
    std::complex<double>* alpha_ptr = m_coeffs->get_alpha_ptr(element);

    double phi_x = phi;
    double phi_y = phi + M_PI/2;
    oskar_evaluate_spherical_wave_sum_double(1, &theta, &phi_x, &phi_y, l_max, alpha_ptr, response_ptr);
}

std::string OSKARElementResponseSphericalWave::get_path(
    const char* filename) const
{
    std::stringstream ss;
    ss << LOFARBEAM_DATA_DIR << "/";
    ss << filename;
    return ss.str();
}


} // namespace StationResponse
} // namespace LOFAR
