//# LofarMetaDataUtil.h: Utility functions to read the meta data relevant for
//# simulating the beam from LOFAR observations stored in MS format.
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

#ifndef LOFAR_STATIONRESPONSE_LOFARMETADATAUTIL_H
#define LOFAR_STATIONRESPONSE_LOFARMETADATAUTIL_H

// \file
// Utility functions to read the meta data relevant for simulating the beam from
// LOFAR observations stored in MS format.

#include <StationResponse/Station.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MeasurementSets/MSAntennaColumns.h>

namespace LOFAR
{
namespace StationResponse
{

// \addtogroup StationResponse
// @{

Station::Ptr readStation(const casa::MeasurementSet &ms, unsigned int id);

template <typename T>
void readStations(const casa::MeasurementSet &ms, T out_it)
{
    casa::ROMSAntennaColumns antenna(ms.antenna());
    for(unsigned int i = 0; i < antenna.nrow(); ++i)
    {
        *out_it++ = readStation(ms, i);
    }
}

// @}

} //# namespace StationResponse
} //# namespace LOFAR

#endif
