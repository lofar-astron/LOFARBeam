#include "BeamFormer.h"

#include <cmath>

constexpr double speed_of_light = 299792458.0;

namespace LOFAR {
namespace StationResponse {

std::vector<std::complex<double>> BeamFormer::compute_geometric_response(double freq, const vector3r_t &direction) const
{
    std::vector<std::complex<double>> result;
    result.reserve(m_antennas.size());
    for (auto &antenna : m_antennas)
    {
        double dl = direction[0] * antenna->m_position[0] +
                    direction[1] * antenna->m_position[1] +
                    direction[2] * antenna->m_position[2];

        double phase = 2 * M_PI * dl * speed_of_light * freq;
        result.push_back({std::sin(phase), std::cos(phase)});
    }
    return result;
}

vector3r_t BeamFormer::compute_local_pointing(double time) const
{
    return {m_pointing[0], m_pointing[1], m_pointing[1]};
}


std::vector<std::pair<std::complex<double>,std::complex<double>>> BeamFormer::compute_weights(double time, double freq) const
{
    std::vector<std::pair<std::complex<double>,std::complex<double>>> result;

    vector3r_t pointing = compute_local_pointing(time);

    auto geometric_response = compute_geometric_response(freq, pointing);
    result.reserve(geometric_response.size());
    for (auto phasor : geometric_response)
    {
        result.push_back({std::conj(phasor), std::conj(phasor)});
    }
    return result;
}


matrix22c_t BeamFormer::response(
    real_t time,
    real_t freq,
    const vector3r_t &direction,
    const Options &options) const
{
    auto weights = compute_weights(time, freq);
    auto geometric_response = compute_geometric_response(freq, direction);

    matrix22c_t result = {0};
    for (unsigned int antenna_idx = 0; antenna_idx < m_antennas.size(); ++antenna_idx) {
        auto antenna = m_antennas[antenna_idx];
        auto antenna_weight = weights[antenna_idx];
        auto antenna_geometric_reponse = geometric_response[antenna_idx];

        matrix22c_t antenna_response = antenna->response(time, freq, direction, options);

        result[0][0] += antenna_weight.first * antenna_geometric_reponse * antenna_response[0][0];
        result[0][1] += antenna_weight.first * antenna_geometric_reponse * antenna_response[0][1];
        result[1][0] += antenna_weight.second * antenna_geometric_reponse * antenna_response[1][0];
        result[1][1] += antenna_weight.second * antenna_geometric_reponse * antenna_response[1][1];

    }
    return result;
}


} // namespace StationResponse
} // namespace LOFAR

