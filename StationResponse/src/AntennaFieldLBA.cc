//# AntennaFieldLBA.cc: Representation of an LBA antenna field.
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
#include <StationResponse/AntennaFieldLBA.h>
#include <StationResponse/MathUtil.h>

namespace LOFAR
{
namespace StationResponse
{

AntennaFieldLBA::AntennaFieldLBA(const string &name,
    const CoordinateSystem &coordinates, const AntennaModelLBA::ConstPtr &model)
    :   AntennaField(name, coordinates),
        itsAntennaModel(model)
{
}

raw_array_factor_t AntennaFieldLBA::rawArrayFactor(real_t, real_t,
    const vector3r_t&, const vector3r_t&) const
{
    raw_array_factor_t af = {{{1.0, 1.0}}, {{1.0, 1.0}}};
    return af;
}

matrix22c_t AntennaFieldLBA::elementResponse(real_t time, real_t freq,
    const vector3r_t &direction) const
{
    return itsAntennaModel->response(freq, itrf2field(direction))
        * rotation(time, direction);
}

} //# namespace StationResponse
} //# namespace LOFAR
