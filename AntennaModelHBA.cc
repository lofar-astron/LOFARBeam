//# AntennaModelHBA.cc: HBA antenna model interface definitions.
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
#include <StationResponse/AntennaModelHBA.h>
#include <StationResponse/MathUtil.h>

namespace LOFAR
{
namespace StationResponse
{

AntennaModelHBA::~AntennaModelHBA()
{
}

matrix22c_t AntennaModelHBA::response(real_t freq, const vector3r_t &direction,
    const vector3r_t &direction0) const
{
    return normalize(rawResponse(freq, direction, direction0));
}

diag22c_t AntennaModelHBA::arrayFactor(real_t freq, const vector3r_t &direction,
    const vector3r_t &direction0) const
{
    return normalize(rawArrayFactor(freq, direction, direction0));
}

raw_response_t AntennaModelHBA::rawResponse(real_t freq,
    const vector3r_t &direction, const vector3r_t &direction0) const
{
    raw_array_factor_t af = rawArrayFactor(freq, direction, direction0);

    raw_response_t result;
    result.response = elementResponse(freq, direction);
    result.response[0][0] *= af.factor[0];
    result.response[0][1] *= af.factor[0];
    result.response[1][0] *= af.factor[1];
    result.response[1][1] *= af.factor[1];
    result.weight = af.weight;
    return result;
}

} //# namespace StationResponse
} //# namespace LOFAR
