//# Constants.h: %Constants used in this library.
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

#ifndef LOFAR_STATIONRESPONSE_CONSTANTS_H
#define LOFAR_STATIONRESPONSE_CONSTANTS_H

// \file
// %Constants used in this library.

#include <StationResponse/Types.h>

namespace LOFAR
{
namespace StationResponse
{

// \addtogroup StationResponse
// @{

/** %Constants used in this library. */
namespace Constants
{
/** 2.0 * pi */
const real_t _2pi = 6.283185307179586476925286;

/** pi / 2.0 */
const real_t pi_2 = 1.570796326794896619231322;

/** pi / 4.0 */
const real_t pi_4 = 0.7853981633974483096156608;

/** Speed of light (m/s) */
const real_t c = 2.99792458e+08;
} //# namespace Constants

// @}

} //# namespace StationResponse
} //# namespace LOFAR

#endif
