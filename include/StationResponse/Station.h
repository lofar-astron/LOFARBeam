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

#include <Common/lofar_smartptr.h>
#include <Common/lofar_string.h>
#include <Common/lofar_vector.h>
#include <StationResponse/AntennaField.h>
#include <StationResponse/Types.h>

namespace LOFAR
{
namespace StationResponse
{

// \addtogroup StationResponse
// @{

class Station
{
public:
    typedef shared_ptr<Station>         Ptr;
    typedef shared_ptr<const Station>   ConstPtr;
    typedef vector<AntennaField::Ptr>   FieldList;

    Station(const string &name, const vector3r_t &position);
    const string &name() const;
    const vector3r_t &position() const;

    void setPhaseReference(const vector3r_t &reference);
    const vector3r_t &phaseReference() const;

    void addAntennaField(const AntennaField::Ptr &field);
    FieldList::const_iterator beginFields() const;
    FieldList::const_iterator endFields() const;

    /**
     *  \brief Compute the response of the station for a plane wave of frequency
     *  \p freq, arriving from direction \p direction, with the %station beam
     *  former steered towards \p direction0, and, for HBA stations, the analog
     *  %tile beam former steered towards \p tile0. For LBA stations, \p tile0
     *  has no effect.
     *
     *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
     *  \param freq Frequency of the plane wave (Hz).
     *  \param direction Direction of arrival (ITRF, m).
     *  \param freq0 %Station beam former reference frequency (Hz).
     *  \param direction0 %Station beam former reference direction (ITRF, m).
     *  \param tile0 Tile beam former reference direction (ITRF, m).
     *  \return Jones matrix that represents the %station response.
     *
     *  For any given sub band, the %LOFAR station beam former computes weights
     *  for a single reference frequency. Usually, this reference frequency is
     *  the center frequency of the sub band. For any frequency except the
     *  reference frequency, these weights are an approximation. This aspect of
     *  the system is taken into account in the computation of the response.
     *  Therefore, both the frequency of interest \p freq and the reference
     *  frequency \p freq0 need to be specified.
     *
     *  The directions \p direction, \p direction0, and \p tile0 are vectors
     *  that represent a direction of \e arrival. These vectors have unit length
     *  and point \e from the ground \e towards the direction from which the
     *  plane wave arrives.
     */
    matrix22c_t response(real_t time, real_t freq, const vector3r_t &direction,
        real_t freq0, const vector3r_t direction0, const vector3r_t &tile0)
        const;

    /**
     *  \brief Compute the array factor of the station for a plane wave of
     *  frequency \p freq, arriving from direction \p direction, with the
     *  %station beam former steered towards \p direction0, and, for HBA
     *  stations, the analog %tile beam former steered towards \p tile0. For LBA
     *  stations, \p tile0 has no effect.
     *
     *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
     *  \param freq Frequency of the plane wave (Hz).
     *  \param direction Direction of arrival (ITRF, m).
     *  \param freq0 %Station beam former reference frequency (Hz).
     *  \param direction0 %Station beam former reference direction (ITRF, m).
     *  \param tile0 Tile beam former reference direction (ITRF, m).
     *  \return A diagonal matrix with the array factor the the X and Y
     *          antennae.
     *
     *  For any given sub band, the %LOFAR station beam former computes weights
     *  for a single reference frequency. Usually, this reference frequency is
     *  the center frequency of the sub band. For any frequency except the
     *  reference frequency, these weights are an approximation. This aspect of
     *  the system is taken into account in the computation of the response.
     *  Therefore, both the frequency of interest \p freq and the reference
     *  frequency \p freq0 need to be specified.
     *
     *  The directions \p direction, \p direction0, and \p tile0 are vectors
     *  that represent a direction of \e arrival. These vectors have unit length
     *  and point \e from the ground \e towards the direction from which the
     *  plane wave arrives.
     */
    diag22c_t arrayFactor(real_t time, real_t freq, const vector3r_t &direction,
        real_t freq0, const vector3r_t direction0, const vector3r_t &tile0)
        const;

    // @{
    template <typename T, typename U>
    void response(unsigned int count, real_t time, T freq,
        const vector3r_t &direction, real_t freq0, const vector3r_t direction0,
        const vector3r_t &tile0, U out) const;

    template <typename T, typename U>
    void arrayFactor(unsigned int count, real_t time, T freq,
        const vector3r_t &direction, real_t freq0, const vector3r_t direction0,
        const vector3r_t &tile0, U out) const;

    template <typename T, typename U>
    void response(unsigned int count, real_t time, T freq,
        const vector3r_t &direction, T freq0, const vector3r_t direction0,
        const vector3r_t &tile0, U out) const;

    template <typename T, typename U>
    void arrayFactor(unsigned int count, real_t time, T freq,
        const vector3r_t &direction, T freq0, const vector3r_t direction0,
        const vector3r_t &tile0, U out) const;
    // @}

private:
    raw_array_factor_t fieldArrayFactor(const AntennaField::ConstPtr &field,
        real_t time, real_t freq, const vector3r_t &direction, real_t freq0,
        const vector3r_t &position0, const vector3r_t &direction0) const;

private:
    string      itsName;
    vector3r_t  itsPosition;
    vector3r_t  itsPhaseReference;
    FieldList   itsFields;
};

// @}

//# ------------------------------------------------------------------------- //
//# - Implementation: Station                                               - //
//# ------------------------------------------------------------------------- //

template <typename T, typename U>
void Station::response(unsigned int count, real_t time, T freq,
    const vector3r_t &direction, real_t freq0, const vector3r_t direction0,
    const vector3r_t &tile0, U out) const
{
    for(unsigned int i = 0; i < count; ++i)
    {
        *out++ = response(time, *freq++, direction, freq0, direction0, tile0);
    }
}

template <typename T, typename U>
void Station::arrayFactor(unsigned int count, real_t time, T freq,
    const vector3r_t &direction, real_t freq0, const vector3r_t direction0,
    const vector3r_t &tile0, U out) const
{
    for(unsigned int i = 0; i < count; ++i)
    {
        *out++ = arrayFactor(time, *freq++, direction, freq0, direction0,
            tile0);
    }
}

template <typename T, typename U>
void Station::response(unsigned int count, real_t time, T freq,
    const vector3r_t &direction, T freq0, const vector3r_t direction0,
    const vector3r_t &tile0, U out) const
{
    for(unsigned int i = 0; i < count; ++i)
    {
        *out++ = response(time, *freq++, direction, *freq0++, direction0,
            tile0);
    }
}

template <typename T, typename U>
void Station::arrayFactor(unsigned int count, real_t time, T freq,
    const vector3r_t &direction, T freq0, const vector3r_t direction0,
    const vector3r_t &tile0, U out) const
{
    for(unsigned int i = 0; i < count; ++i)
    {
        *out++ = arrayFactor(time, *freq++, direction, *freq0++, direction0,
            tile0);
    }
}

} //# namespace StationResponse
} //# namespace LOFAR

#endif
