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

#include "ITRFDirection.h"

#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>

namespace LOFAR {
namespace StationResponse {


    //ITRF position of CS002LBA, just to use a fixed reference
    const vector3r_t ITRFDirection::itsLOFARPosition = {{ 826577.022720000,461022.995082000,5064892.814 }};

  //TODO: initialize converter with a time (and fixed position) and convert specific directions. Needed for wslean as well as for the makestationresponse executable.


ITRFDirection::ITRFDirection(const vector3r_t &position,
    const vector2r_t &direction)
{
    casacore::MVPosition mvPosition(position[0], position[1], position[2]);
    casacore::MPosition mPosition(mvPosition, casacore::MPosition::ITRF);
    itsFrame = casacore::MeasFrame(casacore::MEpoch(), mPosition);

    // Order of angles seems to be longitude (along the equator), lattitude
    // (towards the pole).
    casacore::MVDirection mvDirection(direction[0], direction[1]);
    casacore::MDirection mDirection(mvDirection, casacore::MDirection::J2000);
    itsConverter = casacore::MDirection::Convert(mDirection,
        casacore::MDirection::Ref(casacore::MDirection::ITRF, itsFrame));
}

ITRFDirection::ITRFDirection(const vector2r_t &direction):
  ITRFDirection(itsLOFARPosition, direction)
{
    //create ITRF Direction from fixed stationposition
}

ITRFDirection::ITRFDirection(const vector3r_t &position,
    const vector3r_t &direction)
{
    casacore::MVPosition mvPosition(position[0], position[1], position[2]);
    casacore::MPosition mPosition(mvPosition, casacore::MPosition::ITRF);
    itsFrame = casacore::MeasFrame(casacore::MEpoch(), mPosition);

    casacore::MVDirection mvDirection(direction[0], direction[1], direction[2]);
    casacore::MDirection mDirection(mvDirection, casacore::MDirection::J2000);
    itsConverter = casacore::MDirection::Convert(mDirection,
        casacore::MDirection::Ref(casacore::MDirection::ITRF, itsFrame));
}

ITRFDirection::ITRFDirection(const vector3r_t &direction):
  ITRFDirection(itsLOFARPosition, direction)

{
    //create ITRF Direction from fixed stationposition
}


vector3r_t ITRFDirection::at(real_t time) const
{
    // Cannot use MeasFrame::resetEpoch(Double), because that assumes the
    // argument is UTC in (fractional) days (MJD).
    itsFrame.resetEpoch(casacore::Quantity(time, "s"));

    const casacore::MDirection &mITRF = itsConverter();
    const casacore::MVDirection &mvITRF = mITRF.getValue();

    vector3r_t itrf = {{mvITRF(0), mvITRF(1), mvITRF(2)}};
    return itrf;
}


} //# namespace StationResponse
} // namespace LOFAR
