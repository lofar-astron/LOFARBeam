//# AntennaFieldHBA.cc: Representation of an HBA antenna field.
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
#include <StationResponse/AntennaFieldHBA.h>
#include <StationResponse/MathUtil.h>

namespace LOFAR
{
namespace StationResponse
{

AntennaFieldHBA::AntennaFieldHBA(const string &name,
    const CoordinateSystem &coordinates, const AntennaModelHBA::ConstPtr &model)
    :   AntennaField(name, coordinates),
        itsAntennaModel(model)
{
}

matrix22c_t AntennaFieldHBA::response(real_t time, real_t freq,
    const vector3r_t &direction, const vector3r_t &direction0) const
{
    return itsAntennaModel->response(freq, itrf2field(direction),
        itrf2field(direction0)) * rotation(time, direction);
}

diag22c_t AntennaFieldHBA::arrayFactor(real_t, real_t freq,
    const vector3r_t &direction, const vector3r_t &direction0) const
{
    return itsAntennaModel->arrayFactor(freq, itrf2field(direction),
        itrf2field(direction0));
}

raw_response_t AntennaFieldHBA::rawResponse(real_t time, real_t freq,
    const vector3r_t &direction, const vector3r_t &direction0) const
{
    raw_response_t result = itsAntennaModel->rawResponse(freq,
        itrf2field(direction), itrf2field(direction0));
    result.response = result.response * rotation(time, direction);
    return result;
}

raw_array_factor_t AntennaFieldHBA::rawArrayFactor(real_t, real_t freq,
    const vector3r_t &direction, const vector3r_t &direction0) const
{
    return itsAntennaModel->rawArrayFactor(freq, itrf2field(direction),
        itrf2field(direction0));
}

matrix22c_t AntennaFieldHBA::elementResponse(real_t time, real_t freq,
    const vector3r_t &direction) const
{
    return itsAntennaModel->elementResponse(freq, itrf2field(direction))
        * rotation(time, direction);
}

} //# namespace StationResponse
} //# namespace LOFAR
