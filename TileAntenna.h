//# TileAntenna.h: Semi-analytical model of a LOFAR HBA tile.
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

#ifndef LOFAR_STATIONRESPONSE_TILEANTENNA_H
#define LOFAR_STATIONRESPONSE_TILEANTENNA_H

// \file
// Semi-analytical model of a LOFAR HBA tile.

#include <StationResponse/AntennaModelHBA.h>

namespace LOFAR
{
namespace StationResponse
{

// \addtogroup StationResponse
// @{

class TileAntenna: public AntennaModelHBA
{
public:
    typedef shared_ptr<TileAntenna>         Ptr;
    typedef shared_ptr<const TileAntenna>   ConstPtr;

    typedef static_array<vector3r_t, 16>    TileConfig;

    explicit TileAntenna(const TileConfig &config);

    void setConfig(const TileConfig &config);

    const TileConfig &config() const;

    virtual raw_array_factor_t rawArrayFactor(real_t freq,
        const vector3r_t &direction, const vector3r_t &direction0) const;

    virtual matrix22c_t elementResponse(real_t freq,
        const vector3r_t &direction) const;

private:
    TileConfig  itsConfig;
};

// @}

} //# namespace StationResponse
} //# namespace LOFAR

#endif
