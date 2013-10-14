//# AntennaField.h: Representation of a LOFAR antenna field, with methods to
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

#ifndef LOFAR_STATIONRESPONSE_ARRAYFIELD_H
#define LOFAR_STATIONRESPONSE_ARRAYFIELD_H

// \file
// Representation of a LOFAR antenna field, with methods to compute its response
// to incoming radiation.

#include <Common/lofar_smartptr.h>
#include <Common/lofar_string.h>
#include <Common/lofar_vector.h>
#include <StationResponse/Constants.h>
#include <StationResponse/Types.h>
#include <StationResponse/ITRFDirection.h>

namespace LOFAR
{
namespace StationResponse
{

// \addtogroup StationResponse
// @{

/**
 *  \brief Base class that represents a generic LOFAR antenna field.
 */
class AntennaField
{
public:
    typedef shared_ptr<AntennaField>        Ptr;
    typedef shared_ptr<const AntennaField>  ConstPtr;

    /**
     *  \brief Antenna field coordinate system.
     *
     *  A right handed, cartesian, local coordinate system with coordinate axes
     *  \p p, \p q, and \p r is associated with each antenna field.
     *
     *  The r-axis is orthogonal to the antenna field, and points towards the
     *  local pseudo zenith.
     *
     *  The q-axis is the northern bisector of the \p X and \p Y dipoles, i.e.
     *  it is the reference direction from which the orientation of the dual
     *  dipole antennae is determined. The q-axis points towards the North at
     *  the core. At remote sites it is defined as the intersection of the
     *  antenna field plane and a plane parallel to the meridian plane at the
     *  core. This ensures the reference directions at all sites are similar.
     *
     *  The p-axis is orthogonal to both other axes, and points towards the East
     *  at the core.
     *
     *  The axes and origin of the anntena field coordinate system are expressed
     *  as vectors in the geocentric, cartesian, ITRF coordinate system, in
     *  meters.
     *
     *  \sa "LOFAR Reference Plane and Reference Direction", M.A. Brentjens,
     *  LOFAR-ASTRON-MEM-248.
     */
    struct CoordinateSystem
    {
        struct Axes
        {
            vector3r_t  p;
            vector3r_t  q;
            vector3r_t  r;
        };

        vector3r_t  origin;
        Axes        axes;
    };

    /** A single antenna. */
    struct Antenna
    {
        /**
         *  \brief Position of the antenna relative to the antenna field center
         *  (origin). This is a vector expressed in the geocentric, cartesian,
         *  ITRF coordinate system, in meters.
         */
        vector3r_t  position;

        /**
         *  \brief Status of the \p X and \p Y signal paths of the antenna,
         *  respectively.
         */
        bool        enabled[2];
    };

    typedef vector<Antenna> AntennaList;

    AntennaField(const string &name, const CoordinateSystem &coordinates)
        :   itsName(name),
            itsCoordinateSystem(coordinates),
            itsNCPCacheTime(-1)
    {
        vector3r_t ncp = {{0.0, 0.0, 1.0}};
        itsNCP.reset(new ITRFDirection(position(), ncp));
    }

    virtual ~AntennaField()
    {
    }

    const string &name() const
    {
        return itsName;
    }

    const vector3r_t &position() const
    {
        return itsCoordinateSystem.origin;
    }

    const CoordinateSystem &coordinates() const
    {
        return itsCoordinateSystem;
    }

    void addAntenna(const Antenna &antenna)
    {
        itsAntennae.push_back(antenna);
    }

    size_t nAntennae() const
    {
        return itsAntennae.size();
    }

    const Antenna &antenna(size_t n) const
    {
        return itsAntennae[n];
    }

    Antenna &antenna(size_t n)
    {
        return itsAntennae[n];
    }

    AntennaList::const_iterator beginAntennae() const
    {
        return itsAntennae.begin();
    }

    AntennaList::const_iterator endAntennae() const
    {
        return itsAntennae.end();
    }

    /**
     *  \brief Compute the response of the antenna field for a plane wave of a
     *  given frequency arriving from a given direction.
     *  \param time Time (UTC, s).
     *  \param freq Frequency of the plane wave (Hz).
     *  \param direction Direction of arrival (ITRF, m).
     *  \param freq0 Beamformer reference frequency (Hz).
     *  \param station0 Beamformer reference direction (ITRF, m).
     *  \return A Jones matrix that represents the %antenna field response.
     *
     *  For any given sub band, the %LOFAR station beam former computes weights
     *  for a single reference frequency. Usually, this reference frequency is
     *  the center frequency of the sub band. For any frequency except the
     *  reference frequency, these weights are an approximation. This aspect of
     *  the system is taken into account in the computation of the response.
     *  Therefore, both the frequency of interest \p freq and the reference
     *  frequency \p freq0 need to be specified.
     *
     *  Both \p target and \p reference are vectors that represent a direction
     *  of \e arrival, i.e. these are vectors of unit length that point \e from
     *  the antenna field \e towards the direction from which the plane wave
     *  arrives.
     */
    virtual matrix22c_t response(real_t time, real_t freq,
        const vector3r_t &direction, const vector3r_t &direction0) const;

    virtual diag22c_t arrayFactor(real_t time, real_t freq,
        const vector3r_t &direction, const vector3r_t &direction0) const;

    virtual raw_response_t rawResponse(real_t time, real_t freq,
        const vector3r_t &direction, const vector3r_t &direction0) const;

    virtual raw_array_factor_t rawArrayFactor(real_t time, real_t freq,
        const vector3r_t &direction, const vector3r_t &direction0) const = 0;

    virtual matrix22c_t singleElementResponse(real_t time, real_t freq,
        const vector3r_t &direction) const = 0;

protected:
    matrix22r_t rotation(real_t time, const vector3r_t &direction) const;
    vector3r_t itrf2station(const vector3r_t &itrf) const;

private:
    vector3r_t ncp(real_t time) const;

    string              itsName;
    CoordinateSystem    itsCoordinateSystem;
    AntennaList         itsAntennae;
    ITRFDirection::Ptr  itsNCP;
    mutable real_t      itsNCPCacheTime;
    mutable vector3r_t  itsNCPCacheDirection;
};

// @}

} //# namespace StationResponse
} //# namespace LOFAR

#endif
