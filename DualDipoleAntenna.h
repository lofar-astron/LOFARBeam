//# DualDipoleAntenna.h: Semi-analytical model of a LOFAR LBA dual dipole
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

#ifndef LOFAR_STATIONRESPONSE_DUALDIPOLEANTENNA_H
#define LOFAR_STATIONRESPONSE_DUALDIPOLEANTENNA_H

// \file
// Semi-analytical model of a LOFAR LBA dual dipole antenna.

#include <StationResponse/AntennaModelLBA.h>

namespace LOFAR
{
namespace StationResponse
{

// \addtogroup StationResponse
// @{

class DualDipoleAntenna: public AntennaModelLBA
{
public:
    typedef shared_ptr<DualDipoleAntenna>       Ptr;
    typedef shared_ptr<const DualDipoleAntenna> ConstPtr;

    virtual matrix22c_t response(real_t freq, const vector3r_t &direction)
        const;
};

// @}

} //# namespace StationResponse
} //# namespace LOFAR

#endif
