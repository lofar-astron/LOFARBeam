#ifndef ANTENNA_H
#define ANTENNA_H

#include <complex>
#include <memory>

#include "Types.h"

namespace LOFAR {
namespace StationResponse {

class Antenna
{
public:

    typedef std::shared_ptr<Antenna> Ptr;

    struct Options
    {
        real_t freq0;
        const vector3r_t *station0;
        const vector3r_t *tile0;
    };

    Antenna() {}

    Antenna(double (&position)[3]) :
        m_position{position[0], position[1], position[2]}
    {
    }

    virtual matrix22c_t response(
        real_t time,
        real_t freq,
        const vector3r_t &direction,
        const Options &options = {}) const = 0;

    virtual diag22c_t arrayFactor(
        real_t time,
        real_t freq,
        const vector3r_t &direction,
        const Options &options = {}) const
    { return {1.0, 1.0}; }


    double m_position[3];
    bool m_enabled[2];

};

} // namespace StationResponse
} // namespace LOFAR
#endif
