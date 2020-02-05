//# Types.h: Types used in this library.
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

#ifndef LOFAR_STATIONRESPONSE_TYPES_H
#define LOFAR_STATIONRESPONSE_TYPES_H

// \file
// Types used in this library.

#include <array>
#include <cstring>
#include <ostream>
#include <complex>

namespace LOFAR {
namespace StationResponse
{

// \addtogroup StationResponse
// @{

/** Print the contents of a static array. */
template <typename T, size_t N>
std::ostream &operator<<(std::ostream &out, const std::array<T, N> &obj);

/** Type used for real scalars. */
typedef double                                      real_t;

/** Type used for complex scalars. */
typedef std::complex<double>                        complex_t;

/** Type used for 2-dimensional real vectors. */
typedef std::array<real_t, 2>                     vector2r_t;

/** Type used for 3-dimensional real vectors. */
typedef std::array<real_t, 3>                     vector3r_t;

/** Type used for 2x2 real diagonal matrices. */
typedef std::array<real_t, 2>                     diag22r_t;

/** Type used for 2x2 complex diagonal matrices. */
typedef std::array<complex_t, 2>                  diag22c_t;

/** Type used for 2x2 real matrices. */
typedef std::array<std::array<real_t, 2>, 2>    matrix22r_t;

/** Type used for 2x2 complex matrices. */
typedef std::array<std::array<complex_t, 2>, 2> matrix22c_t;

/** Response of an array of antenna elements. */
struct raw_response_t
{
    /** Combined response of all (enabled) antenna elements in the array. */
    matrix22c_t response;

    /** Number of antenna elements contributing to the combined response, per
     *  polarization.
     */
    diag22r_t   weight;
};

/** Array factor of an array of antenna elements. A wave front of an incident
 *  plane wave will arrive at each antenna element at a potentially different
 *  time. The time of arrival depends on the location of the antenna element and
 *  the direction of arrival of the plane wave. With respect to a pre-defined
 *  phase reference location, there is a (possibly negative) delay between the
 *  arrival of a wave front at a given antenna element and the arrival of the
 *  same wave front at the phase reference location. The array factor is the sum
 *  of the phase shifts due to these delays. It describes the "sensitivity" of
 *  the array as a function of direction.
 */
struct raw_array_factor_t
{
    /** Array factor due to all (enabled) antenna elements in the array. */
    diag22c_t   factor;

    /** Number of antenna elements contributing to the array factor, per
     *  polarization.
     */
    diag22r_t   weight;
};

// @}

template <typename T, size_t N>
std::ostream &operator<<(std::ostream &out, const std::array<T, N> &obj)
{
  print(out, obj.begin(), obj.end());
  return out;
}

} //# namespace StationResponse
} // namespace LOFAR

#endif
