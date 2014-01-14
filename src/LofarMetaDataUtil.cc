//# LofarMetaDataUtil.cc: Utility functions to read the meta data relevant for
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

#include <lofar_config.h>
#include <StationResponse/LofarMetaDataUtil.h>
#include <StationResponse/AntennaFieldLBA.h>
#include <StationResponse/AntennaFieldHBA.h>
#include <StationResponse/MathUtil.h>
#include <StationResponse/TileAntenna.h>
#include <StationResponse/DualDipoleAntenna.h>

#include <measures/Measures/MDirection.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MCPosition.h>
#include <measures/Measures/MeasTable.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

#include <ms/MeasurementSets/MSAntenna.h>
#include <ms/MeasurementSets/MSAntennaParse.h>
#include <ms/MeasurementSets/MSAntennaColumns.h>
#include <ms/MeasurementSets/MSDataDescription.h>
#include <ms/MeasurementSets/MSDataDescColumns.h>
#include <ms/MeasurementSets/MSField.h>
#include <ms/MeasurementSets/MSFieldColumns.h>
#include <ms/MeasurementSets/MSObservation.h>
#include <ms/MeasurementSets/MSObsColumns.h>
#include <ms/MeasurementSets/MSPolarization.h>
#include <ms/MeasurementSets/MSPolColumns.h>
#include <ms/MeasurementSets/MSSpectralWindow.h>
#include <ms/MeasurementSets/MSSpWindowColumns.h>
#include <ms/MeasurementSets/MSSelection.h>

namespace LOFAR
{
namespace StationResponse
{

using namespace casa;

bool hasColumn(const Table &table, const string &column)
{
    return table.tableDesc().isColumn(column);
}

bool hasSubTable(const Table &table, const string &name)
{
    return table.keywordSet().isDefined(name);
}

Table getSubTable(const Table &table, const string &name)
{
    return table.keywordSet().asTable(name);
}

TileAntenna::TileConfig readTileConfig(const Table &table, unsigned int row)
{
    ROArrayQuantColumn<Double> c_tile_offset(table, "TILE_ELEMENT_OFFSET", "m");

    // Read tile configuration for HBA antenna fields.
    Matrix<Quantity> aips_offset = c_tile_offset(row);
    assert(aips_offset.ncolumn() == TileAntenna::TileConfig::size());

    TileAntenna::TileConfig config;
    for(unsigned int i = 0; i < config.size(); ++i)
    {
        config[i][0] = aips_offset(0, i).getValue();
        config[i][1] = aips_offset(1, i).getValue();
        config[i][2] = aips_offset(2, i).getValue();
    }
    return config;
}

void transformToFieldCoordinates(TileAntenna::TileConfig &config,
    const AntennaField::CoordinateSystem::Axes &axes)
{
    for(unsigned int i = 0; i < config.size(); ++i)
    {
        const vector3r_t position = config[i];
        config[i][0] = dot(position, axes.p);
        config[i][1] = dot(position, axes.q);
        config[i][2] = dot(position, axes.r);
    }
}

AntennaField::CoordinateSystem readCoordinateSystem(const Table &table,
    unsigned int id)
{
    ROArrayQuantColumn<Double> c_position(table, "POSITION", "m");
    ROArrayQuantColumn<Double> c_axes(table, "COORDINATE_AXES", "m");

    // Read antenna field center (ITRF).
    Vector<Quantity> aips_position = c_position(id);
    assert(aips_position.size() == 3);

    vector3r_t position = {{aips_position(0).getValue(),
        aips_position(1).getValue(), aips_position(2).getValue()}};

    // Read antenna field coordinate axes (ITRF).
    Matrix<Quantity> aips_axes = c_axes(id);
    assert(aips_axes.shape().isEqual(IPosition(2, 3, 3)));

    vector3r_t p = {{aips_axes(0, 0).getValue(), aips_axes(1, 0).getValue(),
        aips_axes(2, 0).getValue()}};
    vector3r_t q = {{aips_axes(0, 1).getValue(), aips_axes(1, 1).getValue(),
        aips_axes(2, 1).getValue()}};
    vector3r_t r = {{aips_axes(0, 2).getValue(), aips_axes(1, 2).getValue(),
        aips_axes(2, 2).getValue()}};

    AntennaField::CoordinateSystem system = {position, {p, q, r}};
    return system;
}

void readAntennae(const Table &table, unsigned int id,
    const AntennaField::Ptr &field)
{
    ROArrayQuantColumn<Double> c_offset(table, "ELEMENT_OFFSET", "m");
    ROArrayColumn<Bool> c_flag(table, "ELEMENT_FLAG");

    // Read element offsets and flags.
    Matrix<Quantity> aips_offset = c_offset(id);
    assert(aips_offset.shape().isEqual(IPosition(2, 3, aips_offset.ncolumn())));

    Matrix<Bool> aips_flag = c_flag(id);
    assert(aips_flag.shape().isEqual(IPosition(2, 2, aips_offset.ncolumn())));

    for(size_t i = 0; i < aips_offset.ncolumn(); ++i)
    {
        AntennaField::Antenna antenna;
        antenna.position[0] = aips_offset(0, i).getValue();
        antenna.position[1] = aips_offset(1, i).getValue();
        antenna.position[2] = aips_offset(2, i).getValue();
        antenna.enabled[0] = !aips_flag(0, i);
        antenna.enabled[1] = !aips_flag(1, i);
        field->addAntenna(antenna);
    }
}

AntennaField::Ptr readAntennaField(const Table &table, unsigned int id)
{
    AntennaField::Ptr field;
    AntennaField::CoordinateSystem system = readCoordinateSystem(table, id);

    ROScalarColumn<String> c_name(table, "NAME");
    const string &name = c_name(id);

    if(name == "LBA")
    {
        DualDipoleAntenna::Ptr model(new DualDipoleAntenna());
        field = AntennaField::Ptr(new AntennaFieldLBA(name, system, model));
    }
    else
    {
        TileAntenna::TileConfig config = readTileConfig(table, id);
        transformToFieldCoordinates(config, system.axes);

        TileAntenna::Ptr model(new TileAntenna(config));
        field = AntennaField::Ptr(new AntennaFieldHBA(name, system, model));
    }

    readAntennae(table, id, field);
    return field;
}

void readStationPhaseReference(const Table &table, unsigned int id,
    const Station::Ptr &station)
{
    const string columnName("LOFAR_PHASE_REFERENCE");
    if(hasColumn(table, columnName))
    {
        ROScalarMeasColumn<MPosition> c_reference(table, columnName);
        MPosition mReference = MPosition::Convert(c_reference(id),
            MPosition::ITRF)();
        MVPosition mvReference = mReference.getValue();
        vector3r_t reference = {{mvReference(0), mvReference(1),
            mvReference(2)}};

        station->setPhaseReference(reference);
    }
}

Station::Ptr readStation(const MeasurementSet &ms, unsigned int id)
{
    ROMSAntennaColumns antenna(ms.antenna());
    assert(antenna.nrow() > id && !antenna.flagRow()(id));

    // Get station name.
    const string name(antenna.name()(id));

    // Get station position (ITRF).
    MPosition mPosition = MPosition::Convert(antenna.positionMeas()(id),
        MPosition::ITRF)();
    MVPosition mvPosition = mPosition.getValue();
    const vector3r_t position = {{mvPosition(0), mvPosition(1), mvPosition(2)}};

    // Create station.
    Station::Ptr station(new Station(name, position));

    // Read phase reference position (if available).
    readStationPhaseReference(ms.antenna(), id, station);

    // Read antenna field information.
    Table tab_field = getSubTable(ms, "LOFAR_ANTENNA_FIELD");
    tab_field = tab_field(tab_field.col("ANTENNA_ID") == static_cast<Int>(id));

    for(size_t i = 0; i < tab_field.nrow(); ++i)
    {
        station->addField(readAntennaField(tab_field, i));
    }

    return station;
}

} //# namespace StationResponse
} //# namespace LOFAR
