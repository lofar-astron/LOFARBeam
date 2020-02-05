#include "Element.h"
#include "MathUtil.h"

namespace LOFAR {
namespace StationResponse {

matrix22c_t Element::response(
    real_t time,
    real_t freq,
    const vector3r_t &direction,
    const Options &options) const
{
    vector2r_t thetaphi = cart2thetaphi(direction);
    thetaphi[1] -= m_orientation;
    matrix22c_t result;
    static_assert(sizeof(std::complex<double>[2][2]) == sizeof(matrix22c_t));
    m_element_response->response(m_id, freq, thetaphi[0], thetaphi[1], reinterpret_cast<std::complex<double> (&)[2][2]>(result));
    return result;
}

} // namespace StationResponse
} // namespace LOFAR
