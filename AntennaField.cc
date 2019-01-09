//# AntennaField.cc: Representation of a LOFAR antenna field, with methods to
//# compute its response to incoming radiation.
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

#include "AntennaField.h"
#include "Constants.h"
#include "MathUtil.h"
#include "ElementResponse.h"

#include <casacore/measures/Measures/MeasFrame.h>

namespace LOFAR
{
namespace StationResponse
{

AntennaField::AntennaField(const string &name,
    const CoordinateSystem &coordinates)
    :   itsName(name),
        itsCoordinateSystem(coordinates),
        itsNCPCacheTime(-1)
{
    vector3r_t ncp = {{0.0, 0.0, 1.0}};
    itsNCP.reset(new ITRFDirection(ncp));
    vector3r_t ncppol0 = {{1.0, 0.0, 0.0}};
    itsNCPPol0.reset(new ITRFDirection(ncppol0));
}

AntennaField::~AntennaField()
{
}

const string &AntennaField::name() const
{
    return itsName;
}

const vector3r_t &AntennaField::position() const
{
    return itsCoordinateSystem.origin;
}

const AntennaField::CoordinateSystem &AntennaField::coordinates() const
{
    return itsCoordinateSystem;
}

void AntennaField::addAntenna(const Antenna &antenna)
{
    itsAntennae.push_back(antenna);
}

size_t AntennaField::nAntennae() const
{
    return itsAntennae.size();
}

const AntennaField::Antenna &AntennaField::antenna(size_t n) const
{
    return itsAntennae[n];
}

AntennaField::Antenna &AntennaField::antenna(size_t n)
{
    return itsAntennae[n];
}

AntennaField::AntennaList::const_iterator AntennaField::beginAntennae() const
{
    return itsAntennae.begin();
}

AntennaField::AntennaList::const_iterator AntennaField::endAntennae() const
{
    return itsAntennae.end();
}

vector3r_t AntennaField::ncp(real_t time) const
{
    if(time != itsNCPCacheTime)
    {
        itsNCPCacheDirection = itsNCP->at(time);
        itsNCPCacheTime = time;
    }

    return itsNCPCacheDirection;
}

vector3r_t AntennaField::ncppol0(real_t time) const
{
    if(time != itsNCPPol0CacheTime)
    {
        itsNCPPol0CacheDirection = itsNCPPol0->at(time);
        itsNCPPol0CacheTime = time;
    }

    return itsNCPPol0CacheDirection;
}

vector3r_t AntennaField::itrf2field(const vector3r_t &itrf) const
{
    const CoordinateSystem::Axes &axes = itsCoordinateSystem.axes;
    vector3r_t station = {{dot(axes.p, itrf), dot(axes.q, itrf),
        dot(axes.r, itrf)}};
    return station;
}

matrix22r_t AntennaField::rotation(real_t time, const vector3r_t &direction)
    const
{
    //rotation needs to be optional, normally you only want to rotate your coordinatesytem for the center of your (mosaiced) image    
    // Compute the cross product of the NCP and the target direction. This
    // yields a vector tangent to the celestial sphere at the target
    // direction, pointing towards the East (the direction of +Y in the IAU
    // definition, or positive right ascension).
    // Test if the direction is equal to the NCP. If it is, take a random
    // vector orthogonal to v1 (the east is not defined here).
    vector3r_t v1;
    if (std::abs(ncp(time)[0]-direction[0])<1e-9 &&
        std::abs(ncp(time)[1]-direction[1])<1e-9 &&
        std::abs(ncp(time)[2]-direction[2])<1e-9) {
        // Make sure v1 is orthogonal to ncp(time). In the direction of the meridian
        v1 = normalize(ncppol0(time));
    } else {
        v1 = normalize(cross(ncp(time), direction));
    }

    // Compute the cross product of the antenna field normal (R) and the
    // target direction. This yields a vector tangent to the topocentric
    // spherical coordinate system at the target direction, pointing towards
    // the direction of positive phi (which runs East over North around the
    // pseudo zenith).
    // Test if the normal is equal to the target direction. If it is, take
    // a random vector orthogonal to the normal.
    vector3r_t v2;
    if (std::abs(itsCoordinateSystem.axes.r[0]-direction[0])<1e-9 &&
        std::abs(itsCoordinateSystem.axes.r[1]-direction[1])<1e-9 &&
        std::abs(itsCoordinateSystem.axes.r[2]-direction[2])<1e-9)
    {
        // Nothing to be rotated if the direction is equal to zenith
        v2 = v1;
    } else {
        v2 = normalize(cross(itsCoordinateSystem.axes.r, direction));
    }

    // Compute the cosine and sine of the parallactic angle, i.e. the angle
    // between v1 and v2, both tangent to a latitude circle of their
    // respective spherical coordinate systems.
    real_t coschi = dot(v1, v2);
    real_t sinchi;
    if (coschi==1.0)
      sinchi = 0.0;
    else
      sinchi = dot(cross(v1, v2), direction);

    // The input coordinate system is a right handed system with its third
    // axis along the direction of propagation (IAU +Z). The output
    // coordinate system is right handed as well, but its third axis points
    // in the direction of arrival (i.e. exactly opposite).
    //
    // Because the electromagnetic field is always perpendicular to the
    // direction of propagation, we only need to relate the (X, Y) axes of
    // the input system to the corresponding (theta, phi) axes of the output
    // system.
    //
    // To this end, we first rotate the input system around its third axis
    // to align the Y axis with the phi axis. The X and theta axis are
    // parallel after this rotation, but point in opposite directions. To
    // align the X axis with the theta axis, we flip it.
    //
    // The Jones matrix to align the Y axis with the phi axis when these are
    // separated by an angle phi (measured counter-clockwise around the
    // direction of propagation, looking towards the origin), is given by:
    //
    // [ cos(phi)  sin(phi)]
    // [-sin(phi)  cos(phi)]
    //
    // Here, cos(phi) and sin(phi) can be computed directly, without having
    // to compute phi first (see the computation of coschi and sinchi
    // above).
    //
    // Now, sinchi as computed above is opposite to sin(phi), because the
    // direction used in the computation is the direction of arrival instead
    // of the direction of propagation. Therefore, the sign of sinchi needs
    // to be reversed. Furthermore, as explained above, the X axis has to be
    // flipped to align with the theta axis. The Jones matrix returned from
    // this function is therefore given by:
    //
    // [-coschi  sinchi]
    // [ sinchi  coschi]
    matrix22r_t rotation = {{{{-coschi, sinchi}}, {{sinchi, coschi}}}};
    return rotation;
}

matrix22c_t AntennaField::response(real_t time, real_t freq,
    const vector3r_t &direction, const vector3r_t &direction0, const bool rotate) const
{
    return normalize(rawResponse(time, freq, direction, direction0, rotate));
}

diag22c_t AntennaField::arrayFactor(real_t time, real_t freq,
    const vector3r_t &direction, const vector3r_t &direction0) const
{
    return normalize(rawArrayFactor(time, freq, direction, direction0));
}

raw_response_t AntennaField::rawResponse(real_t time, real_t freq,
    const vector3r_t &direction, const vector3r_t &direction0, const bool rotate) const
{
    raw_array_factor_t af = rawArrayFactor(time, freq, direction, direction0);

    raw_response_t result;
    result.response = elementResponse(time, freq, direction, rotate);
    result.response[0][0] *= af.factor[0];
    result.response[0][1] *= af.factor[0];
    result.response[1][0] *= af.factor[1];
    result.response[1][1] *= af.factor[1];
    result.weight = af.weight;
    return result;
}

} //# namespace StationResponse
} //# namespace LOFAR
