//# pystationresponse.cc: python module for StationResponse object.
//# Copyright (C) 2007
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
//# $Id: pystationresponse.cc 33141 2015-12-16 15:10:19Z dijkema $


#include "ITRFDirection.h"
#include "LofarMetaDataUtil.h"
#include "Station.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Containers/ValueHolder.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSColumns.h>
#include <casacore/measures/Measures/MeasTable.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/Measures/MCDirection.h>

#include <casacore/python/Converters/PycExcp.h>
#include <casacore/python/Converters/PycBasicData.h>
#include <casacore/python/Converters/PycValueHolder.h>

#include <boost/python.hpp>
#include <boost/python/args.hpp>

using namespace casacore;
using namespace boost::python;
using namespace LOFAR::StationResponse;

namespace LOFAR
{
namespace BBS
{
  namespace
  {
    /*!
     *  \brief Convert an ITRF position given as a StationResponse::vector3r_t
     *  instance to a casacore::MPosition.
     */
    MPosition toMPositionITRF(const vector3r_t &position);

    /*!
     *  \brief Convert a casacore::MPosition instance to a
     #  StationResponse::vector3r_t instance.
     */
    vector3r_t fromMPosition(const MPosition &position);

    /*!
     *  \brief Convert a casacore::MDirection instance to a
     *  StationResponse::vector3r_t instance.
     */
    vector3r_t fromMDirection(const MDirection &direction);

    /*!
    *  \brief Check if the specified column exists as a column of the specified
    *  table.
    *
    *  \param table The Table instance to check.
    *  \param column The name of the column.
    */
    bool hasColumn(const Table &table, const string &column);

    /*!
    *  \brief Check if the specified sub-table exists as a sub-table of the
    *  specified table.
    *
    *  \param table The Table instance to check.
    *  \param name The name of the sub-table.
    */
    bool hasSubTable(const Table &table, const string &name);

    /*!
    *  \brief Provide access to a sub-table by name.
    *
    *  \param table The Table instance to which the sub-table is associated.
    *  \param name The name of the sub-table.
    */
    Table getSubTable(const Table &table, const string &name);

    /*!
    *  \brief Attempt to read the position of the observatory. If the
    *  observatory position is unknown, the specified default position is
    *  returned.
    *
    *  \param ms MeasurementSet to read the observatory position from.
    *  \param idObservation Identifier that determines of which observation the
    *  observatory position should be read.
    *  \param defaultPosition The position that will be returned if the
    *  observatory position is unknown.
    */
    MPosition readObservatoryPosition(const MeasurementSet &ms,
      unsigned int idObservation, const MPosition &defaultPosition);

    /*!
     *  \brief Read the phase reference direction.
     *
     *  \param ms MeasurementSet to read the phase reference direction from.
     *  \param idField Identifier of the field of which the phase reference
     *  direction should be read.
     */
    MDirection readPhaseReference(const MeasurementSet &ms,
      unsigned int idField);

    /*!
    *  \brief Read the station beam former reference direction.
    *
    *  \param ms MeasurementSet to read the station beam former reference
    *  direction from.
    *  \param idField Identifier of the field of which the station beam former
    *  reference direction should be read.
    */
    MDirection readDelayReference(const MeasurementSet &ms,
      unsigned int idField);

    /*!
    *  \brief Read the station beam former reference direction.
    *
    *  \param ms MeasurementSet to read the tile beam former reference direction
    *  from.
    *  \param idField Identifier of the field of which the tile beam former
    *  reference direction should be read.
    */
    MDirection readTileReference(const MeasurementSet &ms,
      unsigned int idField);
  } //# unnamed namespace

  class PyStationResponse
  {
  public:
    PyStationResponse(const string& msName, bool inverse = false,
        bool useElementResponse = true, bool useArrayFactor = true,
        bool useChanFreq = false);

    // Get the software version.
    string version(const string& type) const;

    // Set the delay reference direction in radians, J2000. The delay reference
    // direction is the direction used by the station beamformer.
    void setRefDelay(double ra, double dec);

    // Get the delay reference direction in meters, ITRF. The delay reference
    // direction is the direction used by the station beamformer.
    ValueHolder getRefDelay(real_t time);

    // Set the tile reference direction in radians, J2000. The tile reference
    // direction is the direction used by the analog tile beamformer and is
    // relevant only for HBA observations.
    void setRefTile(double ra, double dec);

    // Get the tile reference direction in meters, ITRF. The delay reference
    // direction is the direction used by the analog tile beamformer and is
    // relevant only for HBA observations.
    ValueHolder getRefTile(real_t time);

    // Set the direction of interest in radians, J2000. Can and often will be
    // different than the delay and/or tile reference direction.
    void setDirection(double ra, double dec);

    // Get the direction of intereset in meters, ITRF.
    ValueHolder getDirection(real_t time);

    // Compute the LOFAR beam Jones matrices for the given time, station, and/or
    // channel.
    ValueHolder evaluate0(double time);
    ValueHolder evaluate1(double time, int station);
    ValueHolder evaluate2(double time, int station, int channel);
    ValueHolder evaluate3(double time, int station, double freq);
    ValueHolder evaluate4(double time, int station, double freq, const ValueHolder& direction, const ValueHolder& station0, const ValueHolder& tile0); 

  private:
    Matrix<DComplex> evaluate_itrf(
      const Station::ConstPtr &station, double time, double freq, double freq0,
      const vector3r_t &direction, const vector3r_t &station0, 
      const vector3r_t &tile0) const;

    Matrix<DComplex> evaluate(const Station::ConstPtr &station, double time,
      double freq,  double freq0) const;

    Cube<DComplex> evaluate(const Station::ConstPtr &station, double time,
      const Vector<Double> &freq, const Vector<Double> &freq0) const;

    void invert(matrix22c_t &in) const;
    void invert(diag22c_t &in) const;

    //# Data members.
    bool                  itsInverse;
    bool                  itsUseElementResponse;
    bool                  itsUseArrayFactor;
    bool                  itsUseChanFreq;
    Vector<Station::Ptr>  itsStations;
    Vector<Double>        itsChanFreq;
    Vector<Double>        itsRefFreq;

    vector3r_t            itsRefPosition;
    ITRFDirection::Ptr    itsRefDelay;
    ITRFDirection::Ptr    itsRefTile;

    ITRFDirection::Ptr    itsDirection;
  };

  PyStationResponse::PyStationResponse(const string &name, bool inverse,
    bool useElementResponse, bool useArrayFactor, bool useChanFreq)
    : itsInverse(inverse),
      itsUseElementResponse(useElementResponse),
      itsUseArrayFactor(useArrayFactor),
      itsUseChanFreq(useChanFreq)
  {
    MeasurementSet ms(name);

    // Read spectral window id.
    const unsigned int idDataDescription = 0;
    ROMSDataDescColumns desc(ms.dataDescription());
    assert(desc.nrow() > idDataDescription);
    assert(!desc.flagRow()(idDataDescription));

    // Read the spectral information.
    const unsigned int idWindow = desc.spectralWindowId()(idDataDescription);
    ROMSSpWindowColumns window(ms.spectralWindow());
    assert(window.nrow() > idWindow);
    assert(!window.flagRow()(idWindow));

    itsChanFreq = window.chanFreq()(idWindow);
    itsRefFreq = Vector<Double>(itsChanFreq.size(),
      window.refFrequency()(idWindow));

    // Read the station information.
    ROMSAntennaColumns antenna(ms.antenna());
    itsStations.resize(antenna.nrow());
    for(unsigned int i = 0; i < antenna.nrow(); ++i)
    {
      itsStations(i) = readStation(ms, i);
    }

    // Read observatory position. If unknown, default to the position of the
    // first station.
    unsigned int idObservation = 0;
    MPosition refPosition = readObservatoryPosition(ms, idObservation,
      toMPositionITRF(itsStations(0)->position()));
    itsRefPosition = fromMPosition(MPosition::Convert(refPosition,
      MPosition::ITRF)());

    // Read the reference directions.
    unsigned int idField = 0;
    itsRefDelay.reset(new ITRFDirection(itsRefPosition,
      fromMDirection(MDirection::Convert(readDelayReference(ms, idField),
      MDirection::J2000)())));

    itsRefTile.reset(new ITRFDirection(itsRefPosition,
      fromMDirection(MDirection::Convert(readTileReference(ms, idField),
      MDirection::J2000)())));

    itsDirection.reset(new ITRFDirection(itsRefPosition,
      fromMDirection(MDirection::Convert(readPhaseReference(ms, idField),
      MDirection::J2000)())));
  }

  string PyStationResponse::version(const string& type) const
  {
//     return Version::getInfo<pystationresponseVersion>("stationresponse", type);
      return "0.1";
  }

  void PyStationResponse::setRefDelay(double ra, double dec)
  {
    vector2r_t direction = {{ra, dec}};
    itsRefDelay.reset(new ITRFDirection(itsRefPosition, direction));
  }

  ValueHolder PyStationResponse::getRefDelay(real_t time)
  {
    vector3r_t refDelay=itsRefDelay->at(time);
    Vector<Double> result(3);
    result(0)=refDelay[0]; result(1)=refDelay[1]; result(2)=refDelay[2];

    return ValueHolder(result);
  }

  void PyStationResponse::setRefTile(double ra, double dec)
  {
    vector2r_t direction = {{ra, dec}};
    itsRefTile.reset(new ITRFDirection(itsRefPosition, direction));
  }

  ValueHolder PyStationResponse::getRefTile(real_t time)
  {
    vector3r_t refTile=itsRefTile->at(time);
    Vector<Double> result(3);
    result(0)=refTile[0]; result(1)=refTile[1]; result(2)=refTile[2];

    return ValueHolder(result);
  }

  void PyStationResponse::setDirection(double ra, double dec)
  {
    vector2r_t direction = {{ra, dec}};
    itsDirection.reset(new ITRFDirection(itsRefPosition, direction));
  }

  ValueHolder PyStationResponse::getDirection(real_t time)
  {
    vector3r_t direction=itsDirection->at(time);
    Vector<Double> result(3);
    result(0)=direction[0]; result(1)=direction[1]; result(2)=direction[2];

    return ValueHolder(result);
  }

  ValueHolder PyStationResponse::evaluate0(double time)
  {
    Array<DComplex> result(IPosition(4, 2, 2, itsChanFreq.size(),
      itsStations.size()));

    for(unsigned int i = 0; i < itsStations.size(); ++i)
    {
      IPosition start(4, 0, 0, 0, i);
      IPosition end(4, 1, 1, itsChanFreq.size() - 1, i);
      Cube<DComplex> slice = result(start, end).nonDegenerate();
      if(itsUseChanFreq)
      {
        slice = evaluate(itsStations(i), time, itsChanFreq, itsChanFreq);
      }
      else
      {
        slice = evaluate(itsStations(i), time, itsChanFreq, itsRefFreq);
      }
    }

    return ValueHolder(result);
  }

  ValueHolder PyStationResponse::evaluate1(double time, int station)
  {
//     assertSTR(station >= 0 && static_cast<size_t>(station)
//       < itsStations.size(), "invalid station number: " << station);

    if(itsUseChanFreq)
    {
      return ValueHolder(evaluate(itsStations(station), time, itsChanFreq,
        itsChanFreq));
    }

    return ValueHolder(evaluate(itsStations(station), time, itsChanFreq,
      itsRefFreq));
  }

  ValueHolder PyStationResponse::evaluate2(double time, int station,
    int channel)
  {
//     assertSTR(station >= 0 && static_cast<size_t>(station)
//       < itsStations.size(), "invalid station number: " << station);
//     assertSTR(channel >= 0 && static_cast<size_t>(channel)
//       < itsChanFreq.size(), "invalid channel number: " << channel);


    double freq = itsChanFreq(channel);
    if(itsUseChanFreq)
    {
      return ValueHolder(evaluate(itsStations(station), time, freq, freq));
    }

    double freq0 = itsRefFreq(channel);
    return ValueHolder(evaluate(itsStations(station), time, freq, freq0));
  }

  ValueHolder PyStationResponse::evaluate3(double time, int station,
    double freq)
  {
//     assertSTR(station >= 0 && static_cast<size_t>(station)
//       < itsStations.size(), "invalid station number: " << station);

    if(itsUseChanFreq)
    {
      return ValueHolder(evaluate(itsStations(station), time, freq, freq));
    }

    double freq0 = itsRefFreq(0);
    return ValueHolder(evaluate(itsStations(station), time, freq, freq0));
  }

  ValueHolder PyStationResponse::evaluate4(double time, int station, double freq, const ValueHolder& vh_direction, const ValueHolder& vh_station0, const ValueHolder& vh_tile0)
  {
    assert (vh_direction.dataType() == TpArrayDouble);
    assert (vh_station0.dataType() == TpArrayDouble);
    assert (vh_tile0.dataType() == TpArrayDouble);
    Array<Double> arr_dir(vh_direction.asArrayDouble());
    Array<Double> st0_dir(vh_station0.asArrayDouble());
    Array<Double> tile_dir(vh_tile0.asArrayDouble());
    vector3r_t direction={{arr_dir.data()[0],arr_dir.data()[1],arr_dir.data()[2]}};
    vector3r_t station0 ={{st0_dir.data()[0],st0_dir.data()[1],st0_dir.data()[2]}};
    vector3r_t tile0    ={{tile_dir.data()[0],tile_dir.data()[1],tile_dir.data()[2]}};

    if(itsUseChanFreq)
    {
      return ValueHolder(evaluate_itrf(itsStations(station), time, freq, freq, 
                                       direction, station0, tile0));
    }

    double freq0 = itsRefFreq(0);
    return ValueHolder(evaluate_itrf(itsStations(station), time, freq, freq0,
                                     direction, station0, tile0));
  }

  Cube<DComplex> PyStationResponse::evaluate(const Station::ConstPtr &station,
    double time, const Vector<Double> &freq, const Vector<Double> &freq0) const
  {
    Cube<DComplex> result(2, 2, freq.size(), 0.0);
    if(itsUseArrayFactor)
    {
      vector3r_t direction = itsDirection->at(time);
      vector3r_t station0 = itsRefDelay->at(time);
      vector3r_t tile0 = itsRefTile->at(time);

      if(itsUseElementResponse)
      {
        for(unsigned int i = 0; i < freq.size(); ++i)
        {
          matrix22c_t response = station->response(time, freq(i), direction,
            freq0(i), station0, tile0);

          if(itsInverse)
          {
            invert(response);
          }

          result(0, 0, i) = response[0][0];
          result(1, 0, i) = response[0][1];
          result(0, 1, i) = response[1][0];
          result(1, 1, i) = response[1][1];
        }
      }
      else
      {
        for(unsigned int i = 0; i < freq.size(); ++i)
        {
          diag22c_t af = station->arrayFactor(time, freq(i), direction,
            freq0(i), station0, tile0);

          if(itsInverse)
          {
            invert(af);
          }

          result(0, 0, i) = af[0];
          result(1, 1, i) = af[1];
        }
      }
    }
    else if(itsUseElementResponse)
    {
      // For a station with multiple antenna fields, need to select for which
      // field the element response will be evaluated. Here the first field of the
      // station is always selected.
      AntennaField::ConstPtr field = *station->beginFields();

      vector3r_t direction = itsDirection->at(time);
      for(unsigned int i = 0; i < freq.size(); ++i)
      {
        matrix22c_t response = field->elementResponse(time, freq(i),
          direction);

        if(itsInverse)
        {
          invert(response);
        }

        result(0, 0, i) = response[0][0];
        result(1, 0, i) = response[0][1];
        result(0, 1, i) = response[1][0];
        result(1, 1, i) = response[1][1];
      }
    }
    else
    {
      for(unsigned int i = 0; i < freq.size(); ++i)
      {
        result(0, 0, i) = 1.0;
        result(1, 1, i) = 1.0;
      }
    }

    return result;
  }

  Matrix<DComplex> PyStationResponse::evaluate_itrf(
    const Station::ConstPtr &station, double time, double freq, double freq0,
    const vector3r_t &direction, const vector3r_t &station0, 
    const vector3r_t &tile0) const
  {
    Matrix<DComplex> result(2, 2, 0.0);
    if(itsUseArrayFactor)
    {
      if(itsUseElementResponse)
      {
        matrix22c_t response = station->response(time, freq, direction, freq0,
          station0, tile0);

        if(itsInverse)
        {
          invert(response);
        }

        result(0, 0) = response[0][0];
        result(1, 0) = response[0][1];
        result(0, 1) = response[1][0];
        result(1, 1) = response[1][1];
      }
      else
      {
        diag22c_t af = station->arrayFactor(time, freq, direction, freq0,
          station0, tile0);

        if(itsInverse)
        {
          invert(af);
        }

        result(0, 0) = af[0];
        result(1, 1) = af[1];
      }
    }
    else if(itsUseElementResponse)
    {
      // For a station with multiple antenna fields, need to select for which
      // field the element response will be evaluated. Here the first field of
      // the station is always selected.
      AntennaField::ConstPtr field = *station->beginFields();

      matrix22c_t response = field->elementResponse(time, freq,
        direction);

      if(itsInverse)
      {
        invert(response);
      }

      result(0, 0) = response[0][0];
      result(1, 0) = response[0][1];
      result(0, 1) = response[1][0];
      result(1, 1) = response[1][1];
    }
    else
    {
      result(0, 0) = 1.0;
      result(1, 1) = 1.0;
    }

    return result;
  }

  Matrix<DComplex> PyStationResponse::evaluate(const Station::ConstPtr &station,
    double time, double freq,  double freq0) const
  {
    vector3r_t direction;
    vector3r_t station0;
    vector3r_t tile0;
    if (itsUseArrayFactor) {
      direction = itsDirection->at(time);
      station0 = itsRefDelay->at(time);
      tile0 = itsRefTile->at(time);
    } else if (itsUseElementResponse) {
      direction = itsDirection->at(time);
    }
    return evaluate_itrf(station, time, freq, freq0, direction, station0, tile0);
  }

  void PyStationResponse::invert(matrix22c_t &in) const
  {
    complex_t invDet = 1.0 / (in[0][0] * in[1][1] - in[0][1] * in[1][0]);

    complex_t tmp = in[1][1];
    in[1][1] = in[0][0];
    in[0][0] = tmp;

    in[0][0] *= invDet;
    in[0][1] *= -invDet;
    in[1][0] *= -invDet;
    in[1][1] *= invDet;
  }

  void PyStationResponse::invert(diag22c_t &in) const
  {
    DComplex invDet = 1.0 / (in[0] * in[1]);
    DComplex tmp = in[1];
    in[1] = in[0];
    in[0] = tmp;

    in[0] *= invDet;
    in[1] *= invDet;
  }

  // Now define the interface in Boost-Python.
  void pystationresponse()
  {
    class_<PyStationResponse> ("StationResponse",
        init<std::string, bool, bool, bool, bool>())
      .def ("version", &PyStationResponse::version,
        (boost::python::arg("type")="other"))
      .def ("setRefDelay", &PyStationResponse::setRefDelay,
        (boost::python::arg("ra"), boost::python::arg("dec")))
      .def ("getRefDelay", &PyStationResponse::getRefDelay,
        (boost::python::arg("time")))
      .def ("setRefTile", &PyStationResponse::setRefTile,
        (boost::python::arg("ra"), boost::python::arg("dec")))
      .def ("getRefTile", &PyStationResponse::getRefTile,
        (boost::python::arg("time")))
      .def ("setDirection", &PyStationResponse::setDirection,
        (boost::python::arg("ra"), boost::python::arg("dec")))
      .def ("getDirection", &PyStationResponse::getDirection,
        (boost::python::arg("time")))
      .def ("evaluate0", &PyStationResponse::evaluate0,
        (boost::python::arg("time")))
      .def ("evaluate1", &PyStationResponse::evaluate1,
      (boost::python::arg("time"), boost::python::arg("station")))
      .def ("evaluate2", &PyStationResponse::evaluate2,
      (boost::python::arg("time"), boost::python::arg("station"),
         boost::python::arg("channel")))
      .def ("evaluate3", &PyStationResponse::evaluate3,
      (boost::python::arg("time"), boost::python::arg("station"),
         boost::python::arg("freq")))
      .def ("evaluate4", &PyStationResponse::evaluate4,
      (boost::python::arg("time"), boost::python::arg("station"),
         boost::python::arg("freq"), boost::python::arg("direction"),
         boost::python::arg("station0"), boost::python::arg("tile0")))
      ;
  }

  namespace
  {
    MPosition toMPositionITRF(const vector3r_t &position)
    {
      MVPosition mvITRF(position[0], position[1], position[2]);
      return MPosition(mvITRF, MPosition::ITRF);
    }

    vector3r_t fromMPosition(const MPosition &position)
    {
      MVPosition mvPosition = position.getValue();
      vector3r_t result = {{mvPosition(0), mvPosition(1), mvPosition(2)}};
      return result;
    }

    vector3r_t fromMDirection(const MDirection &direction)
    {
      MVDirection mvDirection = direction.getValue();
      vector3r_t result = {{mvDirection(0), mvDirection(1), mvDirection(2)}};
      return result;
    }

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

    MPosition readObservatoryPosition(const MeasurementSet &ms,
      unsigned int idObservation, const MPosition &defaultPosition)
    {
      // Get the instrument position in ITRF coordinates, or use the centroid
      // of the station positions if the instrument position is unknown.
      ROMSObservationColumns observation(ms.observation());
      assert(observation.nrow() > idObservation);
      assert(!observation.flagRow()(idObservation));

      // Read observatory name and try to look-up its position.
      const string observatory = observation.telescopeName()(idObservation);

      // Look-up observatory position, default to specified default position.
      MPosition position(defaultPosition);
      MeasTable::Observatory(position, observatory);
      return position;
    }

    MDirection readPhaseReference(const MeasurementSet &ms,
      unsigned int idField)
    {
      ROMSFieldColumns field(ms.field());
      assert(field.nrow() > idField);
      assert(!field.flagRow()(idField));

      return field.phaseDirMeas(idField);
    }

    MDirection readDelayReference(const MeasurementSet &ms,
      unsigned int idField)
    {
      ROMSFieldColumns field(ms.field());
      assert(field.nrow() > idField);
      assert(!field.flagRow()(idField));

      return field.delayDirMeas(idField);
    }

    MDirection readTileReference(const MeasurementSet &ms,
      unsigned int idField)
    {
      // The MeasurementSet class does not support LOFAR specific columns, so we
      // use ROArrayMeasColumn to read the tile beam reference direction.
      Table tab_field = getSubTable(ms, "FIELD");

      static const String columnName = "LOFAR_TILE_BEAM_DIR";
      if(hasColumn(tab_field, columnName))
      {
        ROArrayMeasColumn<MDirection> c_direction(tab_field, columnName);
        if(c_direction.isDefined(idField))
        {
          return c_direction(idField)(IPosition(1, 0));
        }
      }

      // By default, the tile beam reference direction is assumed to be equal
      // to the station beam reference direction (for backward compatibility,
      // and for non-HBA measurements).
      return readDelayReference(ms, idField);
    }
  } //# unnamed namespace

} //# namespace BBS
} //# namespace LOFAR

// Define the python module itself.
BOOST_PYTHON_MODULE(_stationresponse)
{
  casacore::python::register_convert_excp();
  casacore::python::register_convert_basicdata();
  casacore::python::register_convert_casa_valueholder();

  LOFAR::BBS::pystationresponse();
}
