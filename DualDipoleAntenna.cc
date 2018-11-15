//# DualDipoleAntenna.cc: Semi-analytical model of a LOFAR LBA dual dipole
//# antenna.
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
#include <StationResponse/DualDipoleAntenna.h>
#include <StationResponse/Constants.h>
#include <StationResponse/MathUtil.h>
#include <ElementResponse/ElementResponse.h>

namespace LOFAR
{
namespace StationResponse
{

matrix22c_t DualDipoleAntenna::response(real_t freq,
    const vector3r_t &direction) const
{
    // The positive X dipole direction is SW of the reference orientation,
    // which translates to a phi coordinate of 5/4*pi in the topocentric
    // spherical coordinate system. The phi coordinate is corrected for this
    // offset before evaluating the antenna model.
    vector2r_t thetaphi = cart2thetaphi(direction);
    thetaphi[1] -= 5.0 * Constants::pi_4;
    matrix22c_t response;
    element_response_lba(freq, thetaphi[0], thetaphi[1],
        reinterpret_cast<std::complex<double> (&)[2][2]>(response));
    return response;
}

} //# namespace StationResponse
} //# namespace LOFAR
