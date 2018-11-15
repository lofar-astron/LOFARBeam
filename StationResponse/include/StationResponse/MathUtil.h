//# MathUtil.h: Various mathematical operations on vectors and matrices.
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

#ifndef LOFAR_STATIONRESPONSE_MATHUTIL_H
#define LOFAR_STATIONRESPONSE_MATHUTIL_H

// \file
// Various mathematical operations on vectors and matrices.

#include <StationResponse/Constants.h>
#include <StationResponse/Types.h>
#include <Common/lofar_math.h>

namespace LOFAR
{
namespace StationResponse
{

// \addtogroup StationResponse
// @{

inline real_t dot(const vector3r_t &arg0, const vector3r_t &arg1)
{
    return arg0[0] * arg1[0] + arg0[1] * arg1[1] +  arg0[2] * arg1[2];
}

inline double norm(const vector3r_t &arg0)
{
    return sqrt(dot(arg0, arg0));
}

inline vector3r_t operator*(real_t arg0, const vector3r_t arg1)
{
    vector3r_t result = {{arg0 * arg1[0], arg0 * arg1[1], arg0 * arg1[2]}};
    return result;
}

inline vector3r_t normalize(const vector3r_t &arg0)
{
    return 1.0 / norm(arg0) * arg0;
}

inline vector2r_t cart2thetaphi(const vector3r_t &cart)
{
    real_t r = sqrt(cart[0] * cart[0] + cart[1] * cart[1]);
    vector2r_t thetaphi = {{Constants::pi_2 - atan2(cart[2], r), atan2(cart[1],
        cart[0])}};
    return thetaphi;
}

inline vector3r_t thetaphi2cart(const vector2r_t &thetaphi)
{
    real_t r = sin(thetaphi[0]);
    vector3r_t cart = {{r * cos(thetaphi[1]), r * sin(thetaphi[1]),
        cos(thetaphi[0])}};
    return cart;
}

// returns az, el, r.
inline vector3r_t cart2sph(const vector3r_t &cart)
{
    real_t r = sqrt(cart[0] * cart[0] + cart[1] * cart[1]);

    vector3r_t sph;
    sph[0] = atan2(cart[1], cart[0]);
    sph[1] = atan2(cart[2], r);
    sph[2] = norm(cart);
    return sph;
}

// expects az, el, r.
inline vector3r_t sph2cart(const vector3r_t &sph)
{
    vector3r_t cart = {{sph[2] * cos(sph[1]) * cos(sph[0]), sph[2] * cos(sph[1])
        * sin(sph[0]), sph[2] * sin(sph[1])}};
    return cart;
}

inline matrix22c_t operator*(const matrix22c_t &arg0, const matrix22r_t &arg1)
{
    matrix22c_t result;
    result[0][0] = arg0[0][0] * arg1[0][0] + arg0[0][1] * arg1[1][0];
    result[0][1] = arg0[0][0] * arg1[0][1] + arg0[0][1] * arg1[1][1];
    result[1][0] = arg0[1][0] * arg1[0][0] + arg0[1][1] * arg1[1][0];
    result[1][1] = arg0[1][0] * arg1[0][1] + arg0[1][1] * arg1[1][1];
    return result;
}

inline vector3r_t cross(const vector3r_t &arg0, const vector3r_t &arg1)
{
    vector3r_t result;
    result[0] = arg0[1] * arg1[2] - arg0[2] * arg1[1];
    result[1] = arg0[2] * arg1[0] - arg0[0] * arg1[2];
    result[2] = arg0[0] * arg1[1] - arg0[1] * arg1[0];
    return result;
}

inline vector3r_t operator+(const vector3r_t &arg0, const vector3r_t &arg1)
{
    vector3r_t result = {{arg0[0] + arg1[0], arg0[1] + arg1[1],
        arg0[2] + arg1[2]}};
    return result;
}

inline vector3r_t operator-(const vector3r_t &arg0, const vector3r_t &arg1)
{
    vector3r_t result = {{arg0[0] - arg1[0], arg0[1] - arg1[1],
        arg0[2] - arg1[2]}};
    return result;
}

inline matrix22c_t normalize(const raw_response_t &raw)
{
    matrix22c_t response = {{ {{}}, {{}} }};

    if(raw.weight[0] != 0.0)
    {
        response[0][0] = raw.response[0][0] / raw.weight[0];
        response[0][1] = raw.response[0][1] / raw.weight[0];
    }

    if(raw.weight[1] != 0.0)
    {
        response[1][0] = raw.response[1][0] / raw.weight[1];
        response[1][1] = raw.response[1][1] / raw.weight[1];
    }

    return response;
}

inline diag22c_t normalize(const raw_array_factor_t &raw)
{
    diag22c_t af = {{}};

    if(raw.weight[0] != 0.0)
    {
        af[0] = raw.factor[0] / raw.weight[0];
    }

    if(raw.weight[1] != 0.0)
    {
        af[1] = raw.factor[1] / raw.weight[1];
    }

    return af;
}

// @}

} //# namespace StationResponse
} //# namespace LOFAR

#endif
