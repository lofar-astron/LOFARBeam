//# Station.h: Representation of the station beam former.
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

#ifndef LOFAR_STATIONRESPONSE_STATION_H
#define LOFAR_STATIONRESPONSE_STATION_H

// \file
// Representation of the station beam former.

#include "ElementResponse.h"
#include "Antenna.h"
#include "BeamFormer.h"
#include "Types.h"

#include <memory>
#include <vector>

namespace LOFAR {
namespace StationResponse
{

// \addtogroup StationResponse
// @{

class Station
{
public:
    typedef std::shared_ptr<Station>             Ptr;
    typedef std::shared_ptr<const Station>       ConstPtr;
//     typedef std::vector<AntennaField::ConstPtr>  FieldList;

    /*!
     *  \brief Construct a new Station instance.
     *
     *  \param name Name of the station.
     *  \param position Position of the station (ITRF, m).
     */
    Station(
        const std::string &name,
        const vector3r_t &position,
        const ElementResponseModel model);

    void setModel(const ElementResponseModel model);

    /*!
     *  \brief Return the name of the station.
     */
    const std::string &name() const;

    /*!
     *  \brief Return the position of the station (ITRF, m).
     */
    const vector3r_t &position() const;

    /*!
     *  \brief Set the phase reference position. This is the position where the
     *  delay of the incoming plane wave is assumed to be zero.
     *
     *  \param reference Phase reference position (ITRF, m).
     *
     *  By default, it is assumed the position of the station is also the phase
     *  reference position. Use this method to set the phase reference position
     *  explicitly when this assumption is false.
     */
    void setPhaseReference(const vector3r_t &reference);

    /*!
     *  \brief Return the phase reference position (ITRF, m).
     *
     *  \see Station::setPhaseReference()
     */
    const vector3r_t &phaseReference() const;

    /*!
     *  \brief Add an antenna field to the station.
     *
     *  Physical %LOFAR stations consist of an LBA field, and either one (remote
     *  and international stations) or two (core stations) HBA fields. Virtual
     *  %LOFAR stations can consist of a combination of the antenna fields of
     *  several physical stations.
     *
     *  Use this method to add the appropriate antenna fields to the station.
     */
//     void addField(const AntennaField::ConstPtr &field);


    /*!
     *  \brief Return the number of available antenna fields.
     */
    size_t nFields() const;

    /*!
     *  \brief Return the requested antenna field.
     *
     *  \param i Antenna field number (0-based).
     *  \return An AntennaField::ConstPtr to the requested AntennaField
     *  instance, or an empty AntennaField::ConstPtr if \p i is out of bounds.
     */
//     AntennaField::ConstPtr field(size_t i) const;

    /*!
     *  \brief Return an iterator that points to the beginning of the list of
     *  antenna fields.
     */
//     FieldList::const_iterator beginFields() const;

    /*!
     *  \brief Return an iterator that points to the end of the list of antenna
     *  fields.
     */
//     FieldList::const_iterator endFields() const;

    /*!
     *  \brief Compute the response of the station for a plane wave of frequency
     *  \p freq, arriving from direction \p direction, with the %station beam
     *  former steered towards \p station0, and, for HBA stations, the analog
     *  %tile beam former steered towards \p tile0. For LBA stations, \p tile0
     *  has no effect.
     *
     *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
     *  \param freq Frequency of the plane wave (Hz).
     *  \param direction Direction of arrival (ITRF, m).
     *  \param freq0 %Station beam former reference frequency (Hz).
     *  \param station0 %Station beam former reference direction (ITRF, m).
     *  \param tile0 Tile beam former reference direction (ITRF, m).
     *  \param rotate Boolean deciding if paralactic rotation should be applied.
     *  \return Jones matrix that represents the %station response.
     *
     *  For any given sub-band, the %LOFAR station beam former computes weights
     *  for a single reference frequency. Usually, this reference frequency is
     *  the center frequency of the sub-band. For any frequency except the
     *  reference frequency, these weights are an approximation. This aspect of
     *  the system is taken into account in the computation of the response.
     *  Therefore, both the frequency of interest \p freq and the reference
     *  frequency \p freq0 need to be specified.
     *
     *  The directions \p direction, \p station0, and \p tile0 are vectors that
     *  represent a direction of \e arrival. These vectors have unit length and
     *  point \e from the ground \e towards the direction from which the plane
     *  wave arrives.
     */
    matrix22c_t response(real_t time, real_t freq, const vector3r_t &direction,
        real_t freq0, const vector3r_t &station0, const vector3r_t &tile0, const bool rotate = true)
        const;

    /*!
     *  \brief Compute the array factor of the station for a plane wave of
     *  frequency \p freq, arriving from direction \p direction, with the
     *  %station beam former steered towards \p station0, and, for HBA stations
     *  the analog %tile beam former steered towards \p tile0. For LBA stations,
     *  \p tile0 has no effect.
     *
     *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
     *  \param freq Frequency of the plane wave (Hz).
     *  \param direction Direction of arrival (ITRF, m).
     *  \param freq0 %Station beam former reference frequency (Hz).
     *  \param station0 %Station beam former reference direction (ITRF, m).
     *  \param tile0 Tile beam former reference direction (ITRF, m).
     *  \param rotate Boolean deciding if paralactic rotation should be applied.
     *  \return A diagonal matrix with the array factor of the X and Y antennae.
     *
     *  For any given sub-band, the %LOFAR station beam former computes weights
     *  for a single reference frequency. Usually, this reference frequency is
     *  the center frequency of the sub-band. For any frequency except the
     *  reference frequency, these weights are an approximation. This aspect of
     *  the system is taken into account in the computation of the response.
     *  Therefore, both the frequency of interest \p freq and the reference
     *  frequency \p freq0 need to be specified.
     *
     *  The directions \p direction, \p station0, and \p tile0 are vectors that
     *  represent a direction of \e arrival. These vectors have unit length and
     *  point \e from the ground \e towards the direction from which the plane
     *  wave arrives.
     */
    diag22c_t arrayFactor(real_t time, real_t freq, const vector3r_t &direction,
        real_t freq0, const vector3r_t &station0, const vector3r_t &tile0)
        const;

    /*!
     *  \name Convenience member functions
     *  These member functions perform the same function as the corresponding
     *  non-template member functions, for a list of frequencies or (frequency,
     *  reference frequency) pairs.
     */
    // @{

    /*!
     *  \brief Convenience method to compute the response of the station for a
     *  list of frequencies, and a fixed reference frequency.
     *
     *  \param count Number of frequencies.
     *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
     *  \param freq Input iterator for a list of frequencies (Hz) of length
     *  \p count.
     *  \param direction Direction of arrival (ITRF, m).
     *  \param freq0 %Station beam former reference frequency (Hz).
     *  \param station0 %Station beam former reference direction (ITRF, m).
     *  \param tile0 Tile beam former reference direction (ITRF, m).
     *  \param rotate Boolean deciding if paralactic rotation should be applied.
     *  \param buffer Output iterator with room for \p count instances of type
     *  ::matrix22c_t.
     *
     *  \see response(real_t time, real_t freq, const vector3r_t &direction,
     *  real_t freq0, const vector3r_t &station0, const vector3r_t &tile0) const
     */
    template <typename T, typename U>
    void response(unsigned int count, real_t time, T freq,
        const vector3r_t &direction, real_t freq0, const vector3r_t &station0,
        const vector3r_t &tile0, U buffer, const bool rotate=true) const;

    /*!
     *  \brief Convenience method to compute the array factor of the station for
     *  a list of frequencies, and a fixed reference frequency.
     *
     *  \param count Number of frequencies.
     *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
     *  \param freq Input iterator for a list of frequencies (Hz) of length
     *  \p count.
     *  \param direction Direction of arrival (ITRF, m).
     *  \param freq0 %Station beam former reference frequency (Hz).
     *  \param station0 %Station beam former reference direction (ITRF, m).
     *  \param tile0 Tile beam former reference direction (ITRF, m).
     *  \param rotate Boolean deciding if paralactic rotation should be applied.
     *  \param buffer Output iterator with room for \p count instances of type
     *  ::diag22c_t.
     *
     *  \see arrayFactor(real_t time, real_t freq, const vector3r_t &direction,
     *  real_t freq0, const vector3r_t &station0, const vector3r_t &tile0) const
     */
    template <typename T, typename U>
    void arrayFactor(unsigned int count, real_t time, T freq,
        const vector3r_t &direction, real_t freq0, const vector3r_t &station0,
        const vector3r_t &tile0, U buffer) const;

    /*!
     *  \brief Convenience method to compute the response of the station for a
     *  list of (frequency, reference frequency) pairs.
     *
     *  \param count Number of frequencies.
     *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
     *  \param freq Input iterator for a list of frequencies (Hz) of length
     *  \p count.
     *  \param direction Direction of arrival (ITRF, m).
     *  \param freq0 Input iterator for a list of %Station beam former reference
     *  frequencies (Hz) of length \p count.
     *  \param station0 %Station beam former reference direction (ITRF, m).
     *  \param tile0 Tile beam former reference direction (ITRF, m).
     *  \param rotate Boolean deciding if paralactic rotation should be applied.
     *  \param buffer Output iterator with room for \p count instances of type
     *  ::matrix22c_t.
     *
     *  \see response(real_t time, real_t freq, const vector3r_t &direction,
     *  real_t freq0, const vector3r_t &station0, const vector3r_t &tile0) const
     */
    template <typename T, typename U>
    void response(unsigned int count, real_t time, T freq,
        const vector3r_t &direction, T freq0, const vector3r_t &station0,
        const vector3r_t &tile0, U buffer, const bool rotate=true) const;

    /*!
     *  \brief Convenience method to compute the array factor of the station for
     *  list of (frequency, reference frequency) pairs.
     *
     *  \param count Number of frequencies.
     *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
     *  \param freq Input iterator for a list of frequencies (Hz) of length
     *  \p count.
     *  \param direction Direction of arrival (ITRF, m).
     *  \param freq0 %Station beam former reference frequency (Hz).
     *  \param station0 %Station beam former reference direction (ITRF, m).
     *  \param tile0 Tile beam former reference direction (ITRF, m).
     *  \param rotate Boolean deciding if paralactic rotation should be applied.
     *  \param buffer Output iterator with room for \p count instances of type
     *  ::diag22c_t.
     *
     *  \see arrayFactor(real_t time, real_t freq, const vector3r_t &direction,
     *  real_t freq0, const vector3r_t &station0, const vector3r_t &tile0) const
     */
    template <typename T, typename U>
    void arrayFactor(unsigned int count, real_t time, T freq,
        const vector3r_t &direction, T freq0, const vector3r_t &station0,
        const vector3r_t &tile0, U buffer) const;

    // @}

    // ===================================================================
    // New methods introduced in refactor
    // ==================================================================
    /**
     *  \brief Station coordinate system.
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


    const ElementResponse::Ptr get_element_response() {return itsElementResponse;}

    matrix22c_t elementResponse(real_t time, real_t freq,
        const vector3r_t &direction, const bool rotate = true) const;

    matrix22c_t response(
        real_t time,
        real_t freq,
        const vector3r_t &direction) const
    {
        return itsAntenna->response(time, freq, direction);
    }

    void set_antenna(Antenna::Ptr antenna) {itsAntenna = antenna;}

private:
    std::string itsName;
    vector3r_t  itsPosition;
    vector3r_t  itsPhaseReference;
    ElementResponseModel itsElementResponseModel = ElementResponseModel::Unknown;
    ElementResponse::Ptr itsElementResponse;
    Element::Ptr itsElement;

    Antenna::Ptr itsAntenna;
};

// @}

//# ------------------------------------------------------------------------- //
//# - Implementation: Station                                               - //
//# ------------------------------------------------------------------------- //

template <typename T, typename U>
void Station::response(unsigned int count, real_t time, T freq,
    const vector3r_t &direction, real_t freq0, const vector3r_t &station0,
    const vector3r_t &tile0, U buffer, const bool rotate) const
{
    for(unsigned int i = 0; i < count; ++i)
    {
        *buffer++ = response(time, *freq++, direction, freq0, station0, tile0,
            rotate);
    }
}

template <typename T, typename U>
void Station::arrayFactor(unsigned int count, real_t time, T freq,
    const vector3r_t &direction, real_t freq0, const vector3r_t &station0,
    const vector3r_t &tile0, U buffer) const
{
    for(unsigned int i = 0; i < count; ++i)
    {
        *buffer++ = arrayFactor(time, *freq++, direction, freq0, station0,
            tile0);
    }
}

template <typename T, typename U>
void Station::response(unsigned int count, real_t time, T freq,
    const vector3r_t &direction, T freq0, const vector3r_t &station0,
    const vector3r_t &tile0, U buffer, const bool rotate) const
{
    for(unsigned int i = 0; i < count; ++i)
    {
        *buffer++ = response(time, *freq++, direction, *freq0++, station0,
            tile0, rotate);
    }
}

template <typename T, typename U>
void Station::arrayFactor(unsigned int count, real_t time, T freq,
    const vector3r_t &direction, T freq0, const vector3r_t &station0,
    const vector3r_t &tile0, U buffer) const
{
    for(unsigned int i = 0; i < count; ++i)
    {
        *buffer++ = arrayFactor(time, *freq++, direction, *freq0++, station0,
            tile0);
    }
}

} //# namespace StationResponse
} // namespace LOFAR

#endif
