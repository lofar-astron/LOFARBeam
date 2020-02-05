//# Station.cc: Representation of the station beam former.
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

#include "Station.h"
#include "MathUtil.h"

#include "hamaker/HamakerElementResponse.h"
#include "oskar/OSKARElementResponse.h"
#include "lobes/LOBESElementResponse.h"
// #include "DualDipoleAntenna.h"
// #include "TileAntenna.h"

namespace LOFAR {
namespace StationResponse
{

Station::Station(
    const std::string &name,
    const vector3r_t &position,
    const ElementResponseModel model)
    :   itsName(name),
        itsPosition(position),
        itsPhaseReference(position),
        itsElementResponse(nullptr)
{
    setModel(model);
}

void Station::setModel(const ElementResponseModel model)
{
    switch (model)
    {
        case Hamaker:
            itsElementResponse.set(HamakerElementResponse::getInstance(itsName));
            break;
        case OSKARDipole:
            itsElementResponse.set(OSKARElementResponseDipole::getInstance());
            break;
        case OSKARSphericalWave:
            itsElementResponse.set(OSKARElementResponseSphericalWave::getInstance());
            break;
        case LOBES:
            itsElementResponse.set(LOBESElementResponse::getInstance(itsName));
            break;
        default:
            std::stringstream message;
            message << "The requested element response model '"
                    << model << "' is not implemented.";
            throw std::runtime_error(message.str());
    }
}

const std::string &Station::name() const
{
    return itsName;
}

const vector3r_t &Station::position() const
{
    return itsPosition;
}

void Station::setPhaseReference(const vector3r_t &reference)
{
    itsPhaseReference = reference;
}

const vector3r_t &Station::phaseReference() const
{
    return itsPhaseReference;
}

// ========================================================

matrix22c_t Station::elementResponse(real_t time, real_t freq,
    const vector3r_t &direction, const bool rotate) const
{
    // TODO
//   if (rotate)
//     return itsAntenna->response(freq, itrf2field(direction))
//         * rotation(time, direction);
//   else
//     return itsAntenna->response(freq, itrf2field(direction));
//     return itsElementResponse->response(freq, direction);
}

matrix22c_t Station::response(real_t time, real_t freq,
    const vector3r_t &direction, real_t freq0, const vector3r_t &station0,
    const vector3r_t &tile0, const bool rotate) const
{
    Antenna::Options options = {
        .freq0 = freq0,
        .station0 = &station0,
        .tile0 = &tile0};
    return itsAntenna->response(time, freq, direction);
}

diag22c_t Station::arrayFactor(real_t time, real_t freq,
    const vector3r_t &direction, real_t freq0, const vector3r_t &station0,
    const vector3r_t &tile0) const
{
    Antenna::Options options = {
        .freq0 = freq0,
        .station0 = &station0,
        .tile0 = &tile0};
    return itsAntenna->arrayFactor(time, freq, direction);
}

} //# namespace StationResponse
} // namespace LOFAR
