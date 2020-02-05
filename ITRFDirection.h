//# ITRFDirection.h: Functor that maps time to an ITRF direction.
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

#ifndef LOFAR_STATIONRESPONSE_ITRFDIRECTION_H
#define LOFAR_STATIONRESPONSE_ITRFDIRECTION_H

// \file
// Functor that maps time to an ITRF direction.

#include "Types.h"

#include <casacore/measures/Measures/MeasFrame.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>

#include <memory>

namespace LOFAR {
namespace StationResponse {

// \addtogroup StationResponse
// @{

class ITRFDirection
{
public:
    typedef std::shared_ptr<ITRFDirection>       Ptr;
    typedef std::shared_ptr<const ITRFDirection> ConstPtr;

    ITRFDirection(const vector3r_t &position, const vector2r_t &direction);
    ITRFDirection(const vector3r_t &position, const vector3r_t &direction);
    ITRFDirection(const vector2r_t &direction);
    ITRFDirection(const vector3r_t &direction);

    vector3r_t at(real_t time) const;

    const static vector3r_t& LOFARPosition() { return itsLOFARPosition; }

private:
    //ITRF position of CS002LBA, just to use a fixed reference
    const static vector3r_t itsLOFARPosition;

    mutable casacore::MeasFrame             itsFrame;
    mutable casacore::MDirection::Convert   itsConverter;
};

// @}

} //# namespace StationResponse
} // namespace LOFAR

#endif
