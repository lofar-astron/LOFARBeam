//# TileAntenna.cc: Semi-analytical model of a LOFAR HBA tile.
//#
//# Copyright (C) 2013
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$

#include <lofar_config.h>
#include <StationResponse/TileAntenna.h>
#include <StationResponse/Constants.h>
#include <StationResponse/MathUtil.h>
#include <ElementResponse/ElementResponse.h>

namespace LOFAR
{
namespace StationResponse
{

TileAntenna::TileAntenna(const TileConfig &config)
    :   itsConfig(config)
{
}

void TileAntenna::setConfig(const TileConfig &config)
{
    itsConfig = config;
}

const TileAntenna::TileConfig &TileAntenna::config() const
{
    return itsConfig;
}

raw_array_factor_t TileAntenna::rawArrayFactor(real_t freq,
    const vector3r_t &direction, const vector3r_t &direction0) const
{
    // Angular wave number.
    real_t k = Constants::_2pi * freq / Constants::c;

    // We need to compute the phase difference between a signal arriving from
    // the target direction and a signal arriving from the reference direction,
    // both relative to the center of the tile.
    //
    // This phase difference can be computed using a single dot product per
    // dipole element by exploiting the fact that the dot product is
    // distributive over vector addition:
    //
    // a . b + a . c = a . (b + c)
    //
    vector3r_t difference = direction - direction0;

    complex_t af(0.0, 0.0);
    for(TileConfig::const_iterator element_it = itsConfig.begin(),
        element_end = itsConfig.end(); element_it != element_end; ++element_it)
    {
        // Compute the effective delay for a plane wave approaching from the
        // direction of interest with respect to the position of element i
        // when beam forming in the reference direction using time delays.
        real_t shift = k * dot(difference, *element_it);
        af += complex_t(cos(shift), sin(shift));
    }

    real_t size = itsConfig.size();
    raw_array_factor_t result = {{{af, af}}, {{size, size}}};
    return result;
}

matrix22c_t TileAntenna::elementResponse(real_t freq,
    const vector3r_t &direction) const
{
    // The positive X dipole direction is SW of the reference orientation,
    // which translates to a phi coordinate of 5/4*pi in the topocentric
    // spherical coordinate system. The phi coordinate is corrected for this
    // offset before evaluating the antenna model.
    vector2r_t thetaphi = cart2thetaphi(direction);
    thetaphi[1] -= 5.0 * Constants::pi_4;

    matrix22c_t response;
    element_response_hba(freq, thetaphi[0], thetaphi[1],
        reinterpret_cast<std::complex<double> (&)[2][2]>(response));
    return response;
}

} //# namespace StationResponse
} //# namespace LOFAR
