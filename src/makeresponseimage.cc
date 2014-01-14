//# makeresponseimage.cc: Generate images of the station response for a given
//# MeasurementSet.
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

#include <StationResponse/Package__Version.h>
#include <StationResponse/LofarMetaDataUtil.h>
#include <Common/InputParSet.h>
#include <Common/lofar_sstream.h>
#include <Common/LofarLogger.h>
#include <Common/SystemUtil.h>
#include <Common/Version.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>
#include <coordinates/Coordinates/StokesCoordinate.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>
#include <images/Images/PagedImage.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MCPosition.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MeasTable.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MPosition.h>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MeasurementSets/MSDataDescription.h>
#include <ms/MeasurementSets/MSDataDescColumns.h>
#include <ms/MeasurementSets/MSField.h>
#include <ms/MeasurementSets/MSFieldColumns.h>
#include <ms/MeasurementSets/MSObservation.h>
#include <ms/MeasurementSets/MSObsColumns.h>
#include <ms/MeasurementSets/MSSpectralWindow.h>
#include <ms/MeasurementSets/MSSpWindowColumns.h>
#include <tables/Tables/ExprNode.h>

// There is no wrapped include file lofar_iterator.h.
#include <iterator>

using namespace casa;
using namespace LOFAR;
using namespace LOFAR::StationResponse;
using LOFAR::operator<<;

namespace
{
    /*!
     *  \brief Map of ITRF directions required to compute an image of the
     *  station beam.
     *
     *  The station beam library uses the ITRF coordinate system to express
     *  station positions and source directions. Since the Earth moves with
     *  respect to the sky, the ITRF coordinates of a source vary with time.
     *  This structure stores the ITRF coordinates for the station and tile beam
     *  former reference directions, as well as for a grid of points on the sky,
     *  along with the time for which these ITRF coordinates are valid.
     */
    struct ITRFDirectionMap
    {
        /*!
         *  \brief The time for which this ITRF direction map is valid (MJD(UTC)
         *  in seconds).
         */
        double_t            time0;

        /*!
         *  \brief Station beam former reference direction expressed in ITRF
         *  coordinates.
         */
        vector3r_t          station0;

        /*!
         *  \brief Tile beam former reference direction expressed in ITRF
         *  coordinates.
         */
        vector3r_t          tile0;

        /*!
         *  \brief ITRF coordinates for a grid of points on the sky.
         */
        Matrix<vector3r_t>  directions;
    };

    /*!
     *  \brief Create an ITRFDirectionMap.
     *
     *  \param coordinates Sky coordinate system definition.
     *  \param shape Number of points along the RA and DEC axis.
     *  \param epoch Time for which to compute the ITRF coordinates.
     *  \param position0 Station beam former reference position (phase
     *  reference).
     *  \param station0 Station beam former reference direction (pointing).
     *  \param tile0 Tile beam former reference direction (pointing).
     */
    ITRFDirectionMap makeDirectionMap(const DirectionCoordinate &coordinates,
        const IPosition &shape, const MEpoch &epoch, const MPosition &position0,
        const MDirection &station0, const MDirection &tile0);

    /*!
     *  \brief Create a DirectionCoordinate instance that defines an image
     *  coordinate system on the sky (J2000).
     *
     *  \param reference Direction that corresponds to the origin of the
     *  coordinate system.
     *  \param size Number of points along each axis (RA, DEC). The index of the
     *  origin of the coordinate system is (size / 2, size / 2).
     *  \param delta Angular step size in radians (assumed to be the same for
     *  both axes).
     */
    DirectionCoordinate makeCoordinates(const MDirection &reference,
        unsigned int size, double delta);

    /*!
     *  \brief Convert an ITRF position given as a StationResponse::vector3r_t
     *  instance to a casa::MPosition.
     */
    MPosition toMPositionITRF(const vector3r_t &position);

    /*!
     *  \brief Convert a casa::MPosition instance to a
     *  StationResponse::vector3r_t instance.
     */
    vector3r_t fromMPosition(const MPosition &position);

    /*!
     *  \brief Convert a casa::MDirection instance to a
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
     *  \brief Read the list of unique timestamps.
     *
     *  \param ms MeasurementSet to read the list of unique timestamps from.
     */
    Vector<Double> readUniqueTimes(const MeasurementSet &ms);

    /*!
     *  \brief Read the reference frequency of the subband associated to the
     *  specified data description identifier.
     *
     *  \param ms MeasurementSet to read the reference frequency from.
     *  \param idDataDescription Identifier that determines of which subband the
     *  reference frequency should be read.
     */
    double readFreqReference(const MeasurementSet &ms,
        unsigned int idDataDescription);

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
     *  \param ms MeasurementSet to read the tile beam former reference
     *  direction from.
     *  \param idField Identifier of the field of which the tile beam former
     *  reference direction should be read.
     */
    MDirection readTileReference(const MeasurementSet &ms,
        unsigned int idField);

    /*!
     *  \brief Store the specified cube of pixel data as a CASA image.
     *
     *  \param data Pixel data. The third dimension is assumed to be of length
     *  4, referring to the correlation products XX, XY, YX, YY (in this order).
     *  \param coordinates Sky coordinate system definition.
     *  \param frequency Frequency for which the pixel data is valid (Hz).
     *  \param name File name of the output image.
     */
    template <class T>
    void store(const Cube<T> &data, const DirectionCoordinate &coordinates,
        double frequency, const string &name);

    /*!
     *  \brief Convert a string to a CASA Quantity (value with unit).
     */
    Quantity readQuantity(const String &in);

    /*!
     *  \brief Remove all elements from the range [first, last) that fall
     *  outside the interval [min, max].
     *
     *  This function returns an iterator new_last such that the range [first,
     *  new_last) contains no elements that fall outside the interval [min,
     *  max]. The iterators in the range [new_last, last) are all still
     *  dereferenceable, but the elements that they point to are unspecified.
     *  The order of the elements that are not removed is unchanged.
     */
    template <typename T>
    T filter(T first, T last, int min, int max);
} //# unnamed namespace


int main(int argc, char *argv[])
{
    TEST_SHOW_VERSION(argc, argv, StationResponse);
    INIT_LOGGER(basename(string(argv[0])));
    Version::show<StationResponseVersion>(cout);

    // Parse inputs.
    LOFAR::InputParSet inputs;
    inputs.create("ms", "", "Name of input MeasurementSet", "string");
    inputs.create("stations", "0", "IDs of stations to process", "int vector");
    inputs.create("cellsize", "60arcsec", "Angular pixel size",
        "quantity string");
    inputs.create("size", "256", "Number of pixels along each axis", "int");
    inputs.create("offset", "0s", "Time offset from the start of the MS",
        "quantity string");
    inputs.create("frames", "0", "Number of images that will be generated for"
        " each station (equally spaced over the duration of the MS)", "int");
    inputs.create("abs", "false", "If set to true, store the absolute value of"
        " the beam response instead of the complex value (intended for use with"
        " older versions of casaviewer)", "bool");
    inputs.readArguments(argc, argv);

    vector<int> stationIDs(inputs.getIntVector("stations"));
    Quantity cellsize = readQuantity(inputs.getString("cellsize"));
    unsigned int size = max(inputs.getInt("size"), 1);
    Quantity offset = readQuantity(inputs.getString("offset"));
    unsigned int nFrames = max(inputs.getInt("frames"), 1);
    Bool abs = inputs.getBool("abs");

    // Open MS.
    MeasurementSet ms(inputs.getString("ms"));
    unsigned int idObservation = 0, idField = 0, idDataDescription = 0;

    // Read station meta-data.
    vector<Station::Ptr> stations;
    readStations(ms, std::back_inserter(stations));

    // Remove illegal station indices.
    stationIDs.erase(filter(stationIDs.begin(), stationIDs.end(), 0,
        static_cast<int>(stations.size()) - 1), stationIDs.end());

    // Read unique timestamps
    Table selection =
        ms(ms.col("OBSERVATION_ID") == static_cast<Int>(idObservation)
            && ms.col("FIELD_ID") == static_cast<Int>(idField)
            && ms.col("DATA_DESC_ID") == static_cast<Int>(idDataDescription));
    Vector<Double> time = readUniqueTimes(selection);

    // Read reference frequency.
    double refFrequency = readFreqReference(ms, idDataDescription);

    // Use the position of the first selected station as the array reference
    // position if the observatory position cannot be found.
    MPosition refPosition = readObservatoryPosition(ms, idField,
        toMPositionITRF(stations.front()->position()));

    // Read phase reference direction.
    MDirection refPhase = readPhaseReference(ms, idField);

    // Read delay reference direction.
    MDirection refDelay = readDelayReference(ms, idField);

    // Read tile reference direction.
    MDirection refTile = readTileReference(ms, idField);

    // Create image coordinate system.
    IPosition shape(2, size, size);
    DirectionCoordinate coordinates = makeCoordinates(refPhase, size,
        cellsize.getValue("rad"));

    // Compute station response images.
    Cube<Complex> response(size, size, 4);

    MEpoch refEpoch;
    Quantity refTime(time(0), "s");
    refTime = refTime + offset;
    Quantity deltaTime((time(time.size() - 1) - time(0) - offset.getValue("s"))
        / (nFrames - 1), "s");

    for(size_t j = 0; j < nFrames; ++j)
    {
        cout << "[ Frame: " << (j + 1) << "/" << nFrames << " Offset: +"
            << refTime.getValue() - time(0) << " s ]" << endl;

        // Update reference epoch.
        refEpoch.set(refTime);

        cout << "Creating ITRF direction map... " << flush;
        ITRFDirectionMap directionMap = makeDirectionMap(coordinates, shape,
            refEpoch, refPosition, refDelay, refTile);
        cout << "done." << endl;

        cout << "Computing response images... " << flush;
        for(vector<int>::const_iterator it = stationIDs.begin(),
            end = stationIDs.end(); it != end; ++it)
        {
            Station::ConstPtr station = stations[*it];
            cout << *it << ":" << station->name() << " " << flush;

            for(size_t y = 0; y < size; ++y)
            {
                for(size_t x = 0; x < size; ++x)
                {
                    matrix22c_t E = station->response(directionMap.time0,
                        refFrequency, directionMap.directions(x, y),
                        refFrequency, directionMap.station0,
                        directionMap.tile0);

                    response(x, y, 0) = E[0][0];
                    response(x, y, 1) = E[0][1];
                    response(x, y, 2) = E[1][0];
                    response(x, y, 3) = E[1][1];
                }
            }

            std::ostringstream oss;
            oss << "response-" << station->name() << "-frame-" << (j + 1)
                << ".img";

            if(abs)
            {
                store(Cube<Float>(amplitude(response)), coordinates,
                    refFrequency, oss.str());
            }
            else
            {
                store(response, coordinates, refFrequency, oss.str());
            }
        }
        cout << endl;

        refTime = refTime + deltaTime;
    }

    return 0;
}

namespace
{
    ITRFDirectionMap makeDirectionMap(const DirectionCoordinate &coordinates,
        const IPosition &shape, const MEpoch &epoch, const MPosition &position0,
        const MDirection &station0, const MDirection &tile0)
    {
        ITRFDirectionMap map;

        // Convert from MEpoch to a time in MJD(UTC) in seconds.
        MEpoch mEpochUTC = MEpoch::Convert(epoch, MEpoch::Ref(MEpoch::UTC))();
        MVEpoch mvEpochUTC = mEpochUTC.getValue();
        Quantity qEpochUTC = mvEpochUTC.getTime();
        map.time0 = qEpochUTC.getValue("s");

        // Create conversion engine J2000 => ITRF at the specified epoch.
        MDirection::Convert convertor = MDirection::Convert(MDirection::J2000,
          MDirection::Ref(MDirection::ITRF, MeasFrame(epoch, position0)));

        // Compute station and tile beam former reference directions in ITRF at
        // the specified epoch.
        map.station0 = fromMDirection(convertor(station0));
        map.tile0 = fromMDirection(convertor(tile0));

        // Pre-allocate space for the grid of ITRF directions.
        map.directions.resize(shape);

        // Compute ITRF directions.
        MDirection world;
        Vector<Double> pixel = coordinates.referencePixel();
        for(pixel(1) = 0.0; pixel(1) < shape(1); ++pixel(1))
        {
            for(pixel(0) = 0.0; pixel(0) < shape(0); ++pixel(0))
            {
                // CoordinateSystem::toWorld(): RA range [-pi,pi], DEC range
                // [-pi/2,pi/2].
                if(coordinates.toWorld(world, pixel))
                {
                    map.directions(pixel(0), pixel(1)) =
                        fromMDirection(convertor(world));
                }
            }
        }

        return map;
    }

    DirectionCoordinate makeCoordinates(const MDirection &reference,
        unsigned int size, double delta)
    {
        MDirection referenceJ2000 = MDirection::Convert(reference,
            MDirection::J2000)();
        Quantum<Vector<Double> > referenceAngles = referenceJ2000.getAngle();
        double ra = referenceAngles.getBaseValue()(0);
        double dec = referenceAngles.getBaseValue()(1);

        Matrix<Double> xform(2,2);
        xform = 0.0; xform.diagonal() = 1.0;
        return DirectionCoordinate(MDirection::J2000,
            Projection(Projection::SIN), ra, dec, -delta, delta, xform,
            size / 2, size / 2);
    }

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
        ASSERT(observation.nrow() > idObservation);
        ASSERT(!observation.flagRow()(idObservation));

        // Read observatory name and try to look-up its position.
        const string observatory = observation.telescopeName()(idObservation);

        // Look-up observatory position, default to specified default position.
        MPosition position(defaultPosition);
        MeasTable::Observatory(position, observatory);
        return position;
    }

    Vector<Double> readUniqueTimes(const MeasurementSet &ms)
    {
        Table tab_sorted = ms.sort("TIME", Sort::Ascending, Sort::HeapSort
            | Sort::NoDuplicates);

        ROScalarColumn<Double> c_time(tab_sorted, "TIME");
        return c_time.getColumn();
    }

    double readFreqReference(const MeasurementSet &ms,
        unsigned int idDataDescription)
    {
        ROMSDataDescColumns desc(ms.dataDescription());
        ASSERT(desc.nrow() > idDataDescription);
        ASSERT(!desc.flagRow()(idDataDescription));
        uInt idWindow = desc.spectralWindowId()(idDataDescription);

        ROMSSpWindowColumns window(ms.spectralWindow());
        ASSERT(window.nrow() > idWindow);
        ASSERT(!window.flagRow()(idWindow));

        return window.refFrequency()(idWindow);
    }

    MDirection readPhaseReference(const MeasurementSet &ms,
        unsigned int idField)
    {
        ROMSFieldColumns field(ms.field());
        ASSERT(field.nrow() > idField);
        ASSERT(!field.flagRow()(idField));

        return field.phaseDirMeas(idField);
    }

    MDirection readDelayReference(const MeasurementSet &ms,
        unsigned int idField)
    {
        ROMSFieldColumns field(ms.field());
        ASSERT(field.nrow() > idField);
        ASSERT(!field.flagRow()(idField));

        return field.delayDirMeas(idField);
    }

    MDirection readTileReference(const MeasurementSet &ms, unsigned int idField)
    {
        // The MeasurementSet class does not support LOFAR specific columns, so
        // we use ROArrayMeasColumn to read the tile beam reference direction.
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

    template <class T>
    void store(const Cube<T> &data, const DirectionCoordinate &coordinates,
        double frequency, const string &name)
    {
        ASSERT(data.shape()(2) == 4);

        Vector<Int> stokes(4);
        stokes(0) = Stokes::XX;
        stokes(1) = Stokes::XY;
        stokes(2) = Stokes::YX;
        stokes(3) = Stokes::YY;

        CoordinateSystem csys;
        csys.addCoordinate(coordinates);
        csys.addCoordinate(StokesCoordinate(stokes));
        csys.addCoordinate(SpectralCoordinate(MFrequency::TOPO, frequency, 0.0,
            0.0));

        PagedImage<T> im(TiledShape(IPosition(4, data.shape()(0),
            data.shape()(1), 4, 1)), csys, name);
        im.putSlice(data, IPosition(4, 0, 0, 0, 0));
    }

    Quantity readQuantity(const String &in)
    {
        Quantity result;
        bool status = Quantity::read(result, in);
        ASSERT(status);
        return result;
    }

    template <typename T>
    T filter(T first, T last, int min, int max)
    {
        T new_last = first;
        for(; first != last; ++first)
        {
            if(*first >= min && *first <= max)
            {
                *new_last++ = *first;
            }
        }

        return new_last;
    }
} // unnamed namespace.
