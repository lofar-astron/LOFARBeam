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

#include <StationResponse/Types.h>
#include <Common/lofar_smartptr.h>

#include <measures/Measures/MeasFrame.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MCDirection.h>

namespace LOFAR
{
namespace StationResponse
{

// \addtogroup StationResponse
// @{

class ITRFDirection
{
public:
    typedef shared_ptr<ITRFDirection>       Ptr;
    typedef shared_ptr<const ITRFDirection> ConstPtr;

    ITRFDirection(const vector3r_t &position, const vector2r_t &direction);
    ITRFDirection(const vector3r_t &position, const vector3r_t &direction);

    vector3r_t at(real_t time) const;

private:
    mutable casa::MeasFrame             itsFrame;
    mutable casa::MDirection::Convert   itsConverter;
};

// @}

} //# namespace BBS
} //# namespace LOFAR

#endif
