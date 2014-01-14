//# Station.cc: Representation of the station beam former.
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
#include <StationResponse/Station.h>
#include <StationResponse/MathUtil.h>

namespace LOFAR
{
namespace StationResponse
{

Station::Station(const string &name, const vector3r_t &position)
    :   itsName(name),
        itsPosition(position),
        itsPhaseReference(position)
{
}

const string &Station::name() const
{
    return itsName;
}

const vector3r_t &Station::position() const
{
    return itsPosition;
}

void Station::setPhaseReference(const vector3r_t &reference)
{
    itsPhaseReference = reference;
}

const vector3r_t &Station::phaseReference() const
{
    return itsPhaseReference;
}

void Station::addField(const AntennaField::ConstPtr &field)
{
    itsFields.push_back(field);
}

size_t Station::nFields() const
{
    return itsFields.size();
}

AntennaField::ConstPtr Station::field(size_t i) const
{
    return (i < itsFields.size() ? itsFields[i] : AntennaField::ConstPtr());
}

Station::FieldList::const_iterator Station::beginFields() const
{
    return itsFields.begin();
}

Station::FieldList::const_iterator Station::endFields() const
{
    return itsFields.end();
}

matrix22c_t Station::response(real_t time, real_t freq,
    const vector3r_t &direction, real_t freq0, const vector3r_t &station0,
    const vector3r_t &tile0) const
{
    raw_response_t result = {{{{{}}, {{}}}}, {{}}};
    for(FieldList::const_iterator field_it = beginFields(),
        field_end = endFields(); field_it != field_end; ++field_it)
    {
        raw_array_factor_t field = fieldArrayFactor(*field_it, time, freq,
            direction, freq0, phaseReference(), station0);

        raw_response_t antenna = (*field_it)->rawResponse(time, freq,
            direction, tile0);

        result.response[0][0] += field.factor[0] * antenna.response[0][0];
        result.response[0][1] += field.factor[0] * antenna.response[0][1];
        result.response[1][0] += field.factor[1] * antenna.response[1][0];
        result.response[1][1] += field.factor[1] * antenna.response[1][1];

        result.weight[0] += field.weight[0] * antenna.weight[0];
        result.weight[1] += field.weight[1] * antenna.weight[1];
    }

    return normalize(result);
}

diag22c_t Station::arrayFactor(real_t time, real_t freq,
    const vector3r_t &direction, real_t freq0, const vector3r_t &station0,
    const vector3r_t &tile0) const
{
    raw_array_factor_t af = {{{}}, {{}}};
    for(FieldList::const_iterator field_it = beginFields(),
        field_end = endFields(); field_it != field_end; ++field_it)
    {
        raw_array_factor_t field = fieldArrayFactor(*field_it, time, freq,
            direction, freq0, phaseReference(), station0);

        raw_array_factor_t antenna = (*field_it)->rawArrayFactor(time, freq,
            direction, tile0);

        af.factor[0] += field.factor[0] * antenna.factor[0];
        af.factor[1] += field.factor[1] * antenna.factor[1];
        af.weight[0] += field.weight[0] * antenna.weight[0];
        af.weight[1] += field.weight[1] * antenna.weight[1];
    }

    return normalize(af);
}

raw_array_factor_t
Station::fieldArrayFactor(const AntennaField::ConstPtr &field,
    real_t, real_t freq, const vector3r_t &direction, real_t freq0,
    const vector3r_t &position0, const vector3r_t &direction0) const
{
    real_t k = Constants::_2pi * freq / Constants::c;
    real_t k0 = Constants::_2pi * freq0 / Constants::c;

    vector3r_t offset = field->position() - position0;

    raw_array_factor_t af = {{{}}, {{}}};
    typedef AntennaField::AntennaList AntennaList;
    for(AntennaList::const_iterator antenna_it = field->beginAntennae(),
        antenna_end = field->endAntennae(); antenna_it != antenna_end;
        ++antenna_it)
    {
        if(!antenna_it->enabled[0] && !antenna_it->enabled[1])
        {
            continue;
        }

        vector3r_t position = offset + antenna_it->position;
        real_t phase = k * dot(position, direction) - k0 * dot(position,
            direction0);
        complex_t shift = complex_t(cos(phase), sin(phase));

        if(antenna_it->enabled[0])
        {
            af.factor[0] += shift;
            ++af.weight[0];
        }

        if(antenna_it->enabled[1])
        {
            af.factor[1] += shift;
            ++af.weight[1];
        }
    }

    return af;
}

} //# namespace StationResponse
} //# namespace LOFAR
