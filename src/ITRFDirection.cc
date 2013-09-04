//# ITRFDirection.cc: Functor that maps time to an ITRF direction.
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
#include <StationResponse/ITRFDirection.h>

#include <measures/Measures/MPosition.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MEpoch.h>

namespace LOFAR
{
namespace StationResponse
{

ITRFDirection::ITRFDirection(const vector3r_t &position,
    const vector2r_t &direction)
{
    casa::MVPosition mvPosition(position[0], position[1], position[2]);
    casa::MPosition mPosition(mvPosition, casa::MPosition::ITRF);
    itsFrame = casa::MeasFrame(casa::MEpoch(), mPosition);

    // Order of angles seems to be longitude (along the equator), lattitude
    // (towards the pole).
    casa::MVDirection mvDirection(direction[0], direction[1]);
    casa::MDirection mDirection(mvDirection, casa::MDirection::J2000);
    itsConverter = casa::MDirection::Convert(mDirection,
        casa::MDirection::Ref(casa::MDirection::ITRF, itsFrame));
}

ITRFDirection::ITRFDirection(const vector3r_t &position,
    const vector3r_t &direction)
{
    casa::MVPosition mvPosition(position[0], position[1], position[2]);
    casa::MPosition mPosition(mvPosition, casa::MPosition::ITRF);
    itsFrame = casa::MeasFrame(casa::MEpoch(), mPosition);

    casa::MVDirection mvDirection(direction[0], direction[1], direction[2]);
    casa::MDirection mDirection(mvDirection, casa::MDirection::J2000);
    itsConverter = casa::MDirection::Convert(mDirection,
        casa::MDirection::Ref(casa::MDirection::ITRF, itsFrame));
}

vector3r_t ITRFDirection::at(real_t time) const
{
    // Cannot use MeasFrame::resetEpoch(Double), because that assumes the
    // argument is UTC in (fractional) days (MJD).
    itsFrame.resetEpoch(casa::Quantity(time, "s"));

    const casa::MDirection &mITRF = itsConverter();
    const casa::MVDirection &mvITRF = mITRF.getValue();

    vector3r_t itrf = {{mvITRF(0), mvITRF(1), mvITRF(2)}};
    return itrf;
}

} //# namespace StationResponse
} //# namespace LOFAR
