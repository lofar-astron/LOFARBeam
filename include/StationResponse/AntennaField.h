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


    AntennaField(const string &name, const CoordinateSystem &coordinates);

    virtual ~AntennaField();

    /** Return the name of the antenna field. */
    const string &name() const;

    /** Return the phase reference position of the antenna field. */
    const vector3r_t &position() const;

    /** Return the antenna field coordinate system. */
    const CoordinateSystem &coordinates() const;

    /** Add an antenna to the antenna field. */
    void addAntenna(const Antenna &antenna);

    /** Return the number of antennae in the antenna field. */
    size_t nAntennae() const;

    /*!
     *  \name Antennae accessors
     *  These member functions provide access to the antennae that are part of
     *  the antenna field.
     */
    // @{

    /** Return a read-only reference to the antenna with the requested index. */
    const Antenna &antenna(size_t n) const;

    /** Return a writeable reference to the antenna with the requested index. */
    Antenna &antenna(size_t n);

    /** Return a read-only iterator that points to the first antenna of the
     *  antenna field.
     */
    AntennaList::const_iterator beginAntennae() const;

    /** Return a read-only iterator that points one position past the last
     *  antenna of the antenna field.
     */
    AntennaList::const_iterator endAntennae() const;

    // @}

    /*!
     *  \brief Compute the response of the antenna field for a plane wave of
     *  frequency \p freq, arriving from direction \p direction, with the analog
     *  %tile beam former steered towards \p direction0. For LBA antenna fields,
     *  \p direction0 has no effect.
     *
     *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
     *  \param freq Frequency of the plane wave (Hz).
     *  \param direction Direction of arrival (ITRF, m).
     *  \param direction0 Tile beam former reference direction (ITRF, m).
     *  \return Jones matrix that represents the response of the antenna field.
     *
     *  The directions \p direction, and \p direction0 are vectors that
     *  represent a direction of \e arrival. These vectors have unit length and
     *  point \e from the ground \e towards the direction from which the plane
     *  wave arrives.
     */
    virtual matrix22c_t response(real_t time, real_t freq,
        const vector3r_t &direction, const vector3r_t &direction0) const;

    /*!
     *  \brief Compute the array factor of the antenna field for a plane wave of
     *  frequency \p freq, arriving from direction \p direction, analog %tile
     *  beam former steered towards \p direction0. For LBA antenna fields,
     *  \p direction0 has no effect.
     *
     *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
     *  \param freq Frequency of the plane wave (Hz).
     *  \param direction Direction of arrival (ITRF, m).
     *  \param direction0 Tile beam former reference direction (ITRF, m).
     *  \return A diagonal matrix with the array factor of the X and Y antennae.
     *
     *  The directions \p direction, and \p direction0 are vectors that
     *  represent a direction of \e arrival. These vectors have unit length and
     *  point \e from the ground \e towards the direction from which the plane
     *  wave arrives.
     */
    virtual diag22c_t arrayFactor(real_t time, real_t freq,
        const vector3r_t &direction, const vector3r_t &direction0) const;

    /*!
     *  \brief Compute the response of the antenna field for a plane wave of
     *  frequency \p freq, arriving from direction \p direction, with the analog
     *  %tile beam former steered towards \p direction0. For LBA antenna fields,
     *  \p direction0 has no effect.
     *
     *  This method returns the non-normalized (raw) response. This allows the
     *  response of several antenna fields to be summed together, followed by
     *  normalization of the sum.
     *
     *  \see response(real_t time, real_t freq, const vector3r_t &direction,
     *  const vector3r_t &direction0) const
     */
    virtual raw_response_t rawResponse(real_t time, real_t freq,
        const vector3r_t &direction, const vector3r_t &direction0) const;

    /*!
     *  \brief Compute the array factor of the antenna field for a plane wave of
     *  frequency \p freq, arriving from direction \p direction, analog %tile
     *  beam former steered towards \p direction0. For LBA antenna fields,
     *  \p direction0 has no effect.
     *
     *  This method returns the non-normalized (raw) array factor. This allows
     *  the array factor of several antenna fields to be summed together,
     *  followed by normalization of the sum.
     *
     *  \see diag22c_t arrayFactor(real_t time, real_t freq,
     *  const vector3r_t &direction, const vector3r_t &direction0) const
     */
    virtual raw_array_factor_t rawArrayFactor(real_t time, real_t freq,
        const vector3r_t &direction, const vector3r_t &direction0) const = 0;

    /*!
     *  \brief Compute the response of a single antenna for a plane wave of
     *  frequency \p freq, arriving from direction \p direction.
     *
     */
    virtual matrix22c_t elementResponse(real_t time, real_t freq,
        const vector3r_t &direction) const = 0;

protected:
    /** Compute the parallactic rotation. */
    matrix22r_t rotation(real_t time, const vector3r_t &direction) const;

    /** Transform a vector from ITRF to antenna field coordinates. */
    vector3r_t itrf2field(const vector3r_t &itrf) const;

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
