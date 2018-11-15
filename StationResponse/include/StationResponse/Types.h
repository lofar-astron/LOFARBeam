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

#include <Common/lofar_complex.h>
#include <Common/StreamUtil.h>

namespace LOFAR
{
namespace StationResponse
{

// \addtogroup StationResponse
// @{

/**
 *  \brief Array of static size.
 *
 *  This struct wraps a (unidimensional) C array. Like C arrays, it is possible
 *  to use initializers, for example:
 *  \code
 *  int foo[3] = {1, 2, 3};
 *  static_array<int, 3> bar = {{1, 2, 3}};
 *  \endcode
 *  Note the \e double curly braces.
 *
 *  Unlike C arrays, it is possible to return instances of this struct.
 *  \code
 *  static_array<int, 3> foo()
 *  {
 *      static_array<int, 3> bar = {{1, 2, 3}};
 *      return bar;
 *  }
 *  \endcode
 *  \see boost::array
 */
template <typename T, size_t N>
struct static_array
{
    typedef T       *iterator;
    typedef const T *const_iterator;

    T  __data[N];

    /** Returns the size of the array. */
    static size_t size();

    /**
     *  \brief Read-only access to a specific array element.
     *  \param n Index (0-based) of the element to access.
     *  \return A constant reference to the element at index \p n.
     *
     *  This operator provides array-style element access. No range check is
     *  performed on the index \p n. The result of this operator is undefined
     *  for out of range indices.
     */
    const T &operator[](size_t n) const;

    /**
     *  \brief Access to a specific array element.
     *  \param n Index (0-based) of the element to access.
     *  \return A non-constant reference to the element at index \p n.
     *
     *  This operator provides array-style element access. No range check is
     *  performed on the index \p n. The result of this operator is undefined
     *  for out of range indices.
     */
    T &operator[](size_t n);

    /** Returns an iterator that points to the first element in the array. */
    iterator begin();

    /** Returns an iterator that points one past the last  element in the
     *  array.
     */
    iterator end();

    /** Returns a read-only iterator that points to the first element in the
     *  array.
     */
    const_iterator begin() const;

    /** Returns a read-only iterator that points one past the last element in
     *  the array.
     */
    const_iterator end() const;
};

/** Print the contents of a static array. */
template <typename T, size_t N>
ostream &operator<<(ostream &out, const static_array<T,N> &obj);

/** Type used for real scalars. */
typedef double                                      real_t;

/** Type used for complex scalars. */
typedef dcomplex                                    complex_t;

/** Type used for 2-dimensional real vectors. */
typedef static_array<real_t, 2>                     vector2r_t;

/** Type used for 3-dimensional real vectors. */
typedef static_array<real_t, 3>                     vector3r_t;

/** Type used for 2x2 real diagonal matrices. */
typedef static_array<real_t, 2>                     diag22r_t;

/** Type used for 2x2 complex diagonal matrices. */
typedef static_array<complex_t, 2>                  diag22c_t;

/** Type used for 2x2 real matrices. */
typedef static_array<static_array<real_t, 2>, 2>    matrix22r_t;

/** Type used for 2x2 complex matrices. */
typedef static_array<static_array<complex_t, 2>, 2> matrix22c_t;

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

//# ------------------------------------------------------------------------- //
//# - Implementation: static_array                                          - //
//# ------------------------------------------------------------------------- //

template <typename T, size_t N>
inline const T &static_array<T, N>::operator[](size_t n) const
{
    return __data[n];
}

template <typename T, size_t N>
inline T &static_array<T, N>::operator[](size_t n)
{
    return __data[n];
}

template <typename T, size_t N>
inline typename static_array<T, N>::iterator static_array<T, N>::begin()
{
    return __data;
}

template <typename T, size_t N>
inline typename static_array<T, N>::iterator static_array<T, N>::end()
{
    return __data + N;
}

template <typename T, size_t N>
inline typename static_array<T, N>::const_iterator static_array<T, N>::begin()
    const
{
    return __data;
}

template <typename T, size_t N>
inline typename static_array<T, N>::const_iterator static_array<T, N>::end()
    const
{
    return __data + N;
}

template <typename T, size_t N>
inline size_t static_array<T, N>::size()
{
    return N;
}

template <typename T, size_t N>
ostream &operator<<(ostream &out, const static_array<T,N> &obj)
{
  print(out, obj.begin(), obj.end());
  return out;
}

} //# namespace StationResponse
} //# namespace LOFAR

#endif
