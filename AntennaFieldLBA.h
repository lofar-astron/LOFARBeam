//# AntennaFieldLBA.h: Representation of an LBA antenna field.
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

#ifndef LOFAR_STATIONRESPONSE_ANTENNAFIELDLBA_H
#define LOFAR_STATIONRESPONSE_ANTENNAFIELDLBA_H

// \file
// Representation of an LBA antenna field.

#include <StationResponse/AntennaField.h>
#include <StationResponse/AntennaModelLBA.h>

namespace LOFAR
{
namespace StationResponse
{

// \addtogroup StationResponse
// @{

class AntennaFieldLBA: public AntennaField
{
public:
    typedef shared_ptr<AntennaFieldLBA>         Ptr;
    typedef shared_ptr<const AntennaFieldLBA>   ConstPtr;

    AntennaFieldLBA(const string &name, const CoordinateSystem &coordinates,
        const AntennaModelLBA::ConstPtr &model);

    virtual raw_array_factor_t rawArrayFactor(real_t time, real_t freq,
        const vector3r_t &direction, const vector3r_t &direction0) const;

    virtual matrix22c_t elementResponse(real_t time, real_t freq,
        const vector3r_t &direction) const;

private:
    AntennaModelLBA::ConstPtr   itsAntennaModel;
};

// @}

} //# namespace StationResponse
} //# namespace LOFAR

#endif
