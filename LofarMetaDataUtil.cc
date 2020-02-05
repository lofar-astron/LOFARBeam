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

#include "LofarMetaDataUtil.h"
// #include "AntennaFieldLBA.h"
// #include "AntennaFieldHBA.h"
#include "MathUtil.h"
// #include "TileAntenna.h"
// #include "DualDipoleAntenna.h"

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/Measures/MeasTable.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

#include <cassert>
#include <stdexcept>

#include <casacore/ms/MeasurementSets/MSAntenna.h>
#include <casacore/ms/MSSel/MSSelection.h>
#include <casacore/ms/MSSel/MSAntennaParse.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>
#include <casacore/ms/MeasurementSets/MSDataDescription.h>
#include <casacore/ms/MeasurementSets/MSDataDescColumns.h>
#include <casacore/ms/MeasurementSets/MSField.h>
#include <casacore/ms/MeasurementSets/MSFieldColumns.h>
#include <casacore/ms/MeasurementSets/MSObservation.h>
#include <casacore/ms/MeasurementSets/MSObsColumns.h>
#include <casacore/ms/MeasurementSets/MSPolarization.h>
#include <casacore/ms/MeasurementSets/MSPolColumns.h>
#include <casacore/ms/MeasurementSets/MSSpectralWindow.h>
#include <casacore/ms/MeasurementSets/MSSpWindowColumns.h>

namespace LOFAR {
namespace StationResponse {

using namespace casacore;

typedef std::array<vector3r_t, 16>    TileConfig;

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

TileConfig readTileConfig(const Table &table, unsigned int row)
{
    ROArrayQuantColumn<Double> c_tile_offset(table, "TILE_ELEMENT_OFFSET", "m");

    // Read tile configuration for HBA antenna fields.
    Matrix<Quantity> aips_offset = c_tile_offset(row);
    assert(aips_offset.ncolumn() == TileAntenna::TileConfig::size());

    TileConfig config;
    for(unsigned int i = 0; i < config.size(); ++i)
    {
        config[i][0] = aips_offset(0, i).getValue();
        config[i][1] = aips_offset(1, i).getValue();
        config[i][2] = aips_offset(2, i).getValue();
    }
    return config;
}

void transformToFieldCoordinates(TileConfig &config,
    const Station::CoordinateSystem::Axes &axes)
{
    for(unsigned int i = 0; i < config.size(); ++i)
    {
        const vector3r_t position = config[i];
        config[i][0] = dot(position, axes.p);
        config[i][1] = dot(position, axes.q);
        config[i][2] = dot(position, axes.r);
    }
}

// AntennaField::CoordinateSystem readCoordinateSystemAartfaac(
//     const Table &table, unsigned int id)
// {
//     ROArrayQuantColumn<Double> c_position(table, "POSITION", "m");
//
//     // Read antenna field center (ITRF).
//     Vector<Quantity> aips_position = c_position(id);
//     assert(aips_position.size() == 3);
//
//     vector3r_t position = {{aips_position(0).getValue(),
//         aips_position(1).getValue(), aips_position(2).getValue()}};
//
//     TableRecord keywordset = table.keywordSet();
//     Matrix<double> aips_axes;
//     keywordset.get("AARTFAAC_COORDINATE_AXES", aips_axes);
//     assert(aips_axes.shape().isEqual(IPosition(2, 3, 3)));
//
//     vector3r_t p = {{aips_axes(0, 0), aips_axes(1, 0), aips_axes(2, 0)}};
//     vector3r_t q = {{aips_axes(0, 1), aips_axes(1, 1), aips_axes(2, 1)}};
//     vector3r_t r = {{aips_axes(0, 2), aips_axes(1, 2), aips_axes(2, 2)}};
//
//     AntennaField::CoordinateSystem system = {position, {p, q, r}};
//
//     return system;
// }

Station::CoordinateSystem readCoordinateSystem(const Table &table,
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

    Station::CoordinateSystem coordinate_system = {position, {p, q, r}};
    return coordinate_system;
}

BeamFormer::Ptr make_tile(std::string name, Station::CoordinateSystem coordinate_system,
                          TileConfig tile_config, ElementResponse::Ptr)
{
}

BeamFormer::Ptr readAntennaField(const Table &table, unsigned int id, ElementResponse::Ptr element_response)
{
    BeamFormer::Ptr beam_former(new BeamFormer());
    Station::CoordinateSystem coordinate_system = readCoordinateSystem(table, id);

    ROScalarColumn<String> c_name(table, "NAME");
    ROArrayQuantColumn<Double> c_offset(table, "ELEMENT_OFFSET", "m");
    ROArrayColumn<Bool> c_flag(table, "ELEMENT_FLAG");

    const string &name = c_name(id);

    // Read element offsets and flags.
    Matrix<Quantity> aips_offset = c_offset(id);
    assert(aips_offset.shape().isEqual(IPosition(2, 3, aips_offset.ncolumn())));

    Matrix<Bool> aips_flag = c_flag(id);
    assert(aips_flag.shape().isEqual(IPosition(2, 2, aips_offset.ncolumn())));

    TileConfig tile_config;
    if(name != "LBA") readTileConfig(table, id);

    transformToFieldCoordinates(tile_config, coordinate_system.axes);

    for(size_t i = 0; i < aips_offset.ncolumn(); ++i)
    {
        Antenna::Ptr antenna;
        if(name == "LBA") {
            antenna = Element::Ptr(new Element(id, element_response));
        } else {
            // HBA, HBA0, HBA1
            antenna = make_tile(name, coordinate_system, tile_config, element_response);
        }
        antenna->m_position[0] = aips_offset(0, i).getValue();
        antenna->m_position[1] = aips_offset(1, i).getValue();
        antenna->m_position[2] = aips_offset(2, i).getValue();
        antenna->m_enabled[0] = !aips_flag(0, i);
        antenna->m_enabled[1] = !aips_flag(1, i);
        beam_former->add_antenna(antenna);
    }
    return beam_former;
}



BeamFormer::Ptr readAntennaFieldAartfaac(const Table &table, const string &ant_type,
                                           unsigned int id)
{
    BeamFormer::Ptr field;
//     AntennaField::CoordinateSystem system = readCoordinateSystemAartfaac(table, id);
//
//     if (ant_type == "LBA")
//     {
//         DualDipoleAntenna::Ptr model(new DualDipoleAntenna());
//         field = AntennaField::Ptr(new AntennaFieldLBA(ant_type, system, model));
//     }
//     else // HBA
//     {
//          // TODO: implement this
//          throw std::runtime_error("HBA for Aartfaac is not implemented yet.");
//     }
//
//     // Add only one antenna to the field (no offset, always enabled)
//     AntennaField::Antenna antenna;
//     antenna.position[0] = 0.;
//     antenna.position[1] = 0.;
//     antenna.position[2] = 0.;
//     antenna.enabled[0] = true;
//     antenna.enabled[1] = true;
//
//     field->addAntenna(antenna);

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

Station::Ptr readStation(
    const MeasurementSet &ms,
    unsigned int id,
    const ElementResponseModel model)
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
    Station::Ptr station(new Station(name, position, model));

    // Read phase reference position (if available).
    readStationPhaseReference(ms.antenna(), id, station);

    // Read antenna field information.
    ROScalarColumn<String> telescope_name_col(getSubTable(ms, "OBSERVATION"),
                                              "TELESCOPE_NAME");
    string telescope_name = telescope_name_col(0);

    if (telescope_name == "LOFAR")
    {
        Table tab_field = getSubTable(ms, "LOFAR_ANTENNA_FIELD");
        tab_field = tab_field(tab_field.col("ANTENNA_ID") == static_cast<Int>(id));

        if (tab_field.nrow() == 1)
        {
            // There is only one field
            // The Station will consist of the BeamFormer returned by readAntennaField
            station->set_antenna(readAntennaField(tab_field, 0, station->get_element_response()));
        }
        else
        {
            // There are multiple fields
            // The Station will consist of a BeamFormer that combines the fields
            auto beam_former = BeamFormer::Ptr(new BeamFormer());
            for(size_t i = 0; i < tab_field.nrow(); ++i)
            {
                beam_former->add_antenna(readAntennaField(tab_field, i, station->get_element_response()));
            }
            station->set_antenna(beam_former);
        }
    }
    else if (telescope_name == "AARTFAAC")
    {
        ROScalarColumn<String> ant_type_col(getSubTable(ms, "OBSERVATION"),
                                            "AARTFAAC_ANTENNA_TYPE");
        string ant_type = ant_type_col(0);

        Table tab_field = getSubTable(ms, "ANTENNA");
        station->set_antenna(readAntennaFieldAartfaac(tab_field, ant_type, id));
    }

    return station;
}

MDirection readTileBeamDirection(const casacore::MeasurementSet &ms) {
    MDirection tileBeamDir;

    Table fieldTable = getSubTable(ms, "FIELD");

    if (fieldTable.nrow() != 1) {
        throw std::runtime_error("MS has multiple fields, this does not work with the LOFAR beam library.");
    }

    if (hasColumn(fieldTable, "LOFAR_TILE_BEAM_DIR"))
    {
        ROArrayMeasColumn<MDirection> tileBeamCol(fieldTable,
                                                  "LOFAR_TILE_BEAM_DIR");
        tileBeamDir = *(tileBeamCol(0).data());
    }
    else
    {
      ROArrayMeasColumn<MDirection> tileBeamCol(fieldTable,
                                                "DELAY_DIR");
      tileBeamDir = *(tileBeamCol(0).data());
    }

    return tileBeamDir;
}

} //# namespace StationResponse
} // namespace LOFAR
