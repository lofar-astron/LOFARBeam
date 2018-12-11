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
#include "Dir2ITRF.h"

#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>

namespace LOFAR
{
namespace StationResponse
{
  //TODO: initialize converter with a time (and fixed position) and convert specific directions. Needed for wslean as well as for the makestationresponse executable.


Dir2ITRF::Dir2ITRF(real_t time)
{
    //create ITRF Direction from fixed stationposition
    casacore::MVPosition mvPosition(lofarposition[0], lofarposition[1], lofarposition[2]);
    casacore::MPosition mPosition(mvPosition, casacore::MPosition::ITRF);
    casacore::MEpoch timeEpoch(casacore::Quantity(time, "s"));
    itsFrame = casacore::MeasFrame(timeEpoch, mPosition);

    // Order of angles seems to be longitude (along the equator), lattitude
    // (towards the pole).
    itsConverter = casacore::MDirection::Convert(
        casacore::MDirection::J2000,
        casacore::MDirection::Ref(casacore::MDirection::ITRF, itsFrame));
}

void Dir2ITRF::setTime(real_t time)
{
    // Cannot use MeasFrame::resetEpoch(Double), because that assumes the
    // argument is UTC in (fractional) days (MJD).
    itsFrame.resetEpoch(casacore::Quantity(time, "s"));
}

vector3r_t Dir2ITRF::convert_v(const vector2r_t &direction) const
{

    casacore::MVDirection mvDirection(direction[0], direction[1]);
    casacore::MDirection mDirection(mvDirection, casacore::MDirection::J2000);
    const casacore::MDirection &mITRF = itsConverter(mDirection);
    const casacore::MVDirection &mvITRF = mITRF.getValue();

    vector3r_t itrf = {{mvITRF(0), mvITRF(1), mvITRF(2)}};
    return itrf;
}

vector3r_t Dir2ITRF::convert_v(const vector3r_t &direction) const
{
    casacore::MVDirection mvDirection(direction[0], direction[1], direction[2]);
    casacore::MDirection mDirection(mvDirection, casacore::MDirection::J2000);

    const casacore::MDirection &mITRF = itsConverter(mDirection);
    const casacore::MVDirection &mvITRF = mITRF.getValue();

    vector3r_t itrf = {{mvITRF(0), mvITRF(1), mvITRF(2)}};
    return itrf;
}

 
vector3r_t Dir2ITRF::convert_v(const casa::MDirection &direction) const{
    const casacore::MDirection &mITRF = itsConverter(direction);
    const casacore::MVDirection &mvITRF = mITRF.getValue();

    vector3r_t itrf = {{mvITRF(0), mvITRF(1), mvITRF(2)}};
    return itrf;
   
    
  }

casacore::MDirection Dir2ITRF::convert(const vector2r_t &direction) const
{

    casacore::MVDirection mvDirection(direction[0], direction[1]);
    casacore::MDirection mDirection(mvDirection, casacore::MDirection::J2000);
    const casacore::MDirection mITRF = itsConverter(mDirection);
    return mITRF;
}

casacore::MDirection Dir2ITRF::convert(const vector3r_t &direction) const
{
    casacore::MVDirection mvDirection(direction[0], direction[1], direction[2]);
    casacore::MDirection mDirection(mvDirection, casacore::MDirection::J2000);

    const casacore::MDirection mITRF = itsConverter(mDirection);
    return mITRF;
}

 
casacore::MDirection Dir2ITRF::convert(const casa::MDirection &direction) const{
    const casacore::MDirection &mITRF = itsConverter(direction);
    return mITRF;
   
    
  }


} //# namespace StationResponse
} //# namespace LOFAR
