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

#ifndef LOFAR_STATIONRESPONSE_DIR2ITRF_H
#define LOFAR_STATIONRESPONSE_DIR2ITRF_H

// \file
// Functor that maps J2000 to an ITRF direction.

#include "Types.h"

#include <casacore/measures/Measures/MeasFrame.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>

#include <memory>

namespace LOFAR
{
namespace StationResponse
{

// \addtogroup StationResponse
// @{

class Dir2ITRF
{
public:
    typedef std::shared_ptr<ITRFDirection>       Ptr;
    typedef std::shared_ptr<const ITRFDirection> ConstPtr;

    Dir2ITRF(real_t time);
    
    void setTime(real_t time);
    vector3r_t convert_v(const vector2r_t &direction) const;
    vector3r_t convert_v(const vector3r_t &direction) const;
    vector3r_t convert_v(const casa::MDirection &direction) const;
    casa::MDirection convert(const vector2r_t &direction) const;
    casa::MDirection convert(const vector3r_t &direction) const;
    casa::MDirection convert(const casa::MDirection &direction) const;


private:
    mutable casacore::MeasFrame             itsFrame;
    mutable casacore::MDirection::Convert   itsConverter;
};

// @}

} //# namespace BBS
} //# namespace LOFAR

#endif
