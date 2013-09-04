//# makebeamimage.cc: Generate images of the beam response of multiple stations
//# for a given MS.
//#
//# Copyright (C) 2011
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//# more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$

#include <lofar_config.h>

#include <StationResponse/Package__Version.h>
//#include <StationResponse/StationResponse.h>
#include <StationResponse/LofarMetaDataUtil.h>

//#include <LofarFT/LofarConvolutionFunction.h>

#include <Common/InputParSet.h>
#include <Common/LofarLogger.h>
#include <Common/SystemUtil.h>
#include <Common/Version.h>

#include <images/Images/PagedImage.h>
#include <images/Images/HDF5Image.h>
#include <images/Images/ImageFITSConverter.h>

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

#include <measures/Measures/MDirection.h>
#include <measures/Measures/MEpoch.h>

#include <casa/Arrays/ArrayUtil.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayIter.h>
#include <casa/Utilities/Regex.h>
#include <casa/Utilities/Assert.h>
#include <casa/OS/Directory.h>
#include <casa/OS/File.h>
#include <casa/Exceptions/Error.h>
#include <casa/OS/Timer.h>
#include <casa/OS/PrecTimer.h>
#include <casa/iostream.h>
#include <casa/sstream.h>

#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>
#include <coordinates/Coordinates/StokesCoordinate.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>

#include <iterator>

using namespace casa;
using namespace LOFAR;
using LOFAR::operator<<;

template <class T>
void store(const Matrix<T> &data, const string &name)
{
    Matrix<Double> xform(2, 2);
    xform = 0.0;
    xform.diagonal() = 1.0;
    Quantum<Double> incLon((8.0 / data.shape()(0)) * C::pi / 180.0, "rad");
    Quantum<Double> incLat((8.0 / data.shape()(1)) * C::pi / 180.0, "rad");
    Quantum<Double> refLatLon(45.0 * C::pi / 180.0, "rad");
    DirectionCoordinate dir(MDirection::J2000, Projection(Projection::SIN),
        refLatLon, refLatLon, incLon, incLat, xform, data.shape()(0) / 2,
        data.shape()(1) / 2);
    store(dir, data, name);
}

template <class T>
void store (const DirectionCoordinate &dir, const Matrix<T> &data,
const string &name)
{
    //cout<<"Saving... "<<name<<endl;
    Vector<Int> stokes(1);
    stokes(0) = Stokes::I;
    CoordinateSystem csys;
    csys.addCoordinate(dir);
    csys.addCoordinate(StokesCoordinate(stokes));
    csys.addCoordinate(SpectralCoordinate(casa::MFrequency::TOPO, 60e6, 0.0,
        0.0, 60e6));
    PagedImage<T> im(TiledShape(IPosition(4, data.shape()(0), data.shape()(1),
        1, 1)), csys, name);
    im.putSlice(data, IPosition(4, 0, 0, 0, 0));
}

template <class T>
void store(const Cube<T> &data, const string &name)
{
    Matrix<Double> xform(2, 2);
    xform = 0.0;
    xform.diagonal() = 1.0;
    Quantum<Double> incLon((8.0 / data.shape()(0)) * C::pi / 180.0, "rad");
    Quantum<Double> incLat((8.0 / data.shape()(1)) * C::pi / 180.0, "rad");
    Quantum<Double> refLatLon(45.0 * C::pi / 180.0, "rad");
    DirectionCoordinate dir(MDirection::J2000, Projection(Projection::SIN),
        refLatLon, refLatLon, incLon, incLat, xform, data.shape()(0) / 2,
        data.shape()(1) / 2);
    store(dir, data, name);
}

template <class T>
void store(const DirectionCoordinate &dir, const Cube<T> &data,
const string &name)
{
    ASSERT(data.shape()(2) == 4);
    //cout<<"Saving... "<<name<<endl;
    Vector<Int> stokes(4);
    stokes(0) = Stokes::XX;
    stokes(1) = Stokes::XY;
    stokes(2) = Stokes::YX;
    stokes(3) = Stokes::YY;
    CoordinateSystem csys;
    csys.addCoordinate(dir);
    csys.addCoordinate(StokesCoordinate(stokes));
    csys.addCoordinate(SpectralCoordinate(casa::MFrequency::TOPO, 60e6, 0.0,
        0.0, 60e6));
    PagedImage<T> im(TiledShape(IPosition(4, data.shape()(0), data.shape()(1),
        4, 1)), csys, name);
    im.putSlice(data, IPosition(4, 0, 0, 0, 0));
}

Quantity readQuantity (const String& in)
{
  Quantity res;
  bool status = Quantity::read(res, in);
  ASSERT(status);
  return res;
}

struct ITRFDirectionMap
{
    casa::MEpoch                epoch;
    StationResponse::vector3r_t refDelay;
    StationResponse::vector3r_t refTile;
    casa::Cube<casa::Double>    directions;
};

ITRFDirectionMap makeDirectionMap(const MEpoch &epoch,
    const MPosition &refPosition, const MDirection &refDelay,
    const MDirection &refTile,
    const DirectionCoordinate &coordinates,
    const IPosition &shape)
{
    ASSERT(shape[0] > 0 && shape[1] > 0);

    ITRFDirectionMap map;
    map.epoch = epoch;

    // Create conversion engine J2000 -> ITRF at epoch.
    MDirection::Convert convertor = MDirection::Convert(MDirection::J2000,
    MDirection::Ref(MDirection::ITRF,
    MeasFrame(epoch, refPosition)));

    MVDirection mvRefDelay = convertor(refDelay).getValue();
    map.refDelay[0] = mvRefDelay(0);
    map.refDelay[1] = mvRefDelay(1);
    map.refDelay[2] = mvRefDelay(2);

    MVDirection mvRefTile = convertor(refTile).getValue();
    map.refTile[0] = mvRefTile(0);
    map.refTile[1] = mvRefTile(1);
    map.refTile[2] = mvRefTile(2);

    MDirection world;
    casa::Vector<Double> pixel = coordinates.referencePixel();

    Cube<Double> mapITRF(3, shape[0], shape[1], 0.0);
    for(pixel[1] = 0.0; pixel(1) < shape[1]; ++pixel[1])
    {
        for(pixel[0] = 0.0; pixel[0] < shape[0]; ++pixel[0])
        {
            // CoodinateSystem::toWorld()
            // DEC range [-pi/2,pi/2]
            // RA range [-pi,pi]
            if(coordinates.toWorld(world, pixel))
            {
                MVDirection mvITRF = convertor(world).getValue();
                mapITRF(0, pixel[0], pixel[1]) = mvITRF(0);
                mapITRF(1, pixel[0], pixel[1]) = mvITRF(1);
                mapITRF(2, pixel[0], pixel[1]) = mvITRF(2);
            }
        }
    }

    map.directions.reference(mapITRF);
    return map;
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

MDirection readPhaseReference(const MeasurementSet &ms, unsigned int idField)
{
    ROMSFieldColumns field(ms.field());
    ASSERT(field.nrow() > idField);
    ASSERT(!field.flagRow()(idField));

    return field.phaseDirMeas(idField);
}

MDirection readDelayReference(const MeasurementSet &ms, unsigned int idField)
{
    ROMSFieldColumns field(ms.field());
    ASSERT(field.nrow() > idField);
    ASSERT(!field.flagRow()(idField));

    return field.delayDirMeas(idField);
}

MDirection readTileReference(const MeasurementSet &ms, unsigned int idField)
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

int main(int argc, char *argv[])
{
    TEST_SHOW_VERSION(argc, argv, StationResponse);
    INIT_LOGGER(basename(string(argv[0])));
    Version::show<StationResponseVersion>(cout);

    LOFAR::InputParSet inputs;
    inputs.create("ms", "", "Name of input MeasurementSet", "string");
    inputs.create("stations", "[0]", "IDs of stations to process",
        "int vector");
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

    String msName = inputs.getString("ms");
    ASSERT(!msName.empty());

    vector<int> stationID(inputs.getIntVector("stations"));
    Quantity cellsize = readQuantity(inputs.getString("cellsize"));
    size_t size = max(inputs.getInt("size"), 1);
    Quantity offset = readQuantity(inputs.getString("offset"));
    size_t nFrames = max(inputs.getInt("frames"), 1);
    Bool abs = inputs.getBool("abs");

    // ---------------------------------------------------------------------- //

    MeasurementSet ms(msName);

    uInt idObservation = 0;
    uInt idField = 0;
    uInt idDataDescription = 0;

    // Read number of stations.
    ROMSAntennaColumns antenna(ms.antenna());
    uInt nStation = antenna.nrow();

    // Filter invalid station IDs.
    vector<unsigned int> filteredID;
    for(vector<int>::const_iterator it = stationID.begin(),
        end = stationID.end(); it != end; ++it)
    {
        if(*it >= 0 && static_cast<size_t>(*it) < nStation)
        {
            filteredID.push_back(*it);
        }
    }

    // Read phase reference direction.
    MDirection refDir = readPhaseReference(ms, idField);

    // Read delay reference direction.
    MDirection refDelay = readDelayReference(ms, idField);

    // Read tile reference direction.
    MDirection refTile = readTileReference(ms, idField);

    // Read reference frequency.
    ROMSDataDescColumns desc(ms.dataDescription());
    ASSERT(desc.nrow() > idDataDescription);
    ASSERT(!desc.flagRow()(idDataDescription));
    uInt idWindow = desc.spectralWindowId()(idDataDescription);

    ROMSSpWindowColumns window(ms.spectralWindow());
    ASSERT(window.nrow() > idWindow);
    ASSERT(!window.flagRow()(idWindow));

    double refFreq = window.refFrequency()(idWindow);

    // Read reference time.
    Table msView =
        ms(ms.col("OBSERVATION_ID") == static_cast<Int>(idObservation)
        && ms.col("FIELD_ID") == static_cast<Int>(idField)
        && ms.col("DATA_DESC_ID") == static_cast<Int>(idDataDescription));

    Table tab_sorted = msView.sort("TIME", Sort::Ascending,
        Sort::HeapSort | Sort::NoDuplicates);

    ROScalarColumn<Double> c_time(tab_sorted, "TIME");
    Vector<Double> time = c_time.getColumn();

    // ---------------------------------------------------------------------- //

    MDirection refDirJ2000(MDirection::Convert(refDir, MDirection::J2000)());
    Quantum<Vector<Double> > angles = refDirJ2000.getAngle();

    double ra = angles.getBaseValue()(0);
    double dec = angles.getBaseValue()(1);
    double delta = cellsize.getValue("rad");

    // Construct DirectionCoordinate instance.
    Matrix<Double> xform(2,2);
    xform = 0.0; xform.diagonal() = 1.0;
    DirectionCoordinate coordinates(MDirection::J2000,
                            Projection(Projection::SIN),
                            ra, dec,
                            -delta, delta,
                            xform,
                            size / 2, size / 2);

    // ---------------------------------------------------------------------- //

//    StationResponse::Instrument::Ptr instrument = StationResponse::readInstrument(ms,
//        idObservation);


    vector<StationResponse::Station::Ptr> antennae;
    StationResponse::readStations(ms, std::back_inserter(antennae));

//    ROMSAntennaColumns aips_antenna(ms.antenna());
//    for(unsigned int i = 0; i < aips_antenna.nrow(); ++i)
//    {
//        antennae.push_back(StationResponse::readStationLBA(ms, i));
//    }

    MVPosition mvRefPos(0.0, 0.0, 0.0);
    for(size_t i = 0; i < antennae.size(); ++i)
    {
        mvRefPos(0) += antennae[i]->position()[0];
        mvRefPos(1) += antennae[i]->position()[1];
        mvRefPos(2) += antennae[i]->position()[2];
    }
    mvRefPos(0) /= static_cast<double>(antennae.size());
    mvRefPos(1) /= static_cast<double>(antennae.size());
    mvRefPos(2) /= static_cast<double>(antennae.size());

    MPosition refPos(mvRefPos, MPosition::ITRF);

//    MPosition refPos(MVPosition(instrument->position()[0], instrument->position()[1], instrument->position()[2]), MPosition::ITRF);
//    LofarATerm aTerm(ms, Record());

    IPosition shape(2, size, size);
//    aTerm.setDirection(coordinates, shape);

    Vector<Double> freq(1, refFreq);

    Quantity refTime(time(0), "s");
    refTime = refTime + offset;

    Quantity deltaTime((time(time.size() - 1) - time(0) - offset.getValue("s"))
        / (nFrames - 1), "s");

//    StationResponse::ITRFDirectionFuncCasa::Ptr itrfdircasa(new StationResponse::ITRFDirectionFuncCasa(refDelay, refPos));

//    cout << "computing..." << endl;
    for(size_t j = 0; j < nFrames; ++j)
    {
        cout << "frame: " << j << " " << flush;

        MEpoch refEpoch;
        refEpoch.set(refTime);
//        aTerm.setEpoch(refEpoch);

        cout << "map..." << flush;
        ITRFDirectionMap _map = makeDirectionMap(refEpoch, refPos, refDelay,
            refTile, coordinates, shape);

        cout << "done. " << flush;

//        StationResponse::vector3r_t itrfdircasaat = itrfdircasa->at(refTime.getValue());
//        cout << "ITRF: " << itrfdircasaat[0] << " " << itrfdircasaat[1] << " " << itrfdircasaat[2] << " ";
//        cout << "ITRF: " << _map.refDelay[0] << " " << _map.refDelay[1] << " " << _map.refDelay[2] << " ";
//        cout << flush;

        for(vector<unsigned int>::const_iterator it = filteredID.begin(),
            end = filteredID.end(); it != end; ++it)

        {
            cout << "station: " << *it << " response..." << flush;
//            StationResponse::Station::ConstPtr station = instrument->station(*it);
//            StationResponse::AntennaField::ConstPtr field = station->field(0);
//            vector<Cube<Complex> > response = aTerm.evaluate(*it, freq, freq);

            StationResponse::Station::ConstPtr antenna = antennae[*it];

            StationResponse::matrix22c_t resp;
            Cube<Complex> response(size, size, 4);
//            Matrix<Complex> response(size, size);
            for(size_t y = 0; y < size; ++y)
            for(size_t x = 0; x < size; ++x)
            {
                StationResponse::vector3r_t target = {{_map.directions(0, x, y), _map.directions(1, x, y), _map.directions(2, x, y)}};

//                StationResponse::diag22c_t af = array_factor_test(field,
//                    refFreq, target, refFreq, _map.refDelay);
//                response(x, y) = af[0];

//                StationResponse::matrix22c_t resp = station_response(station,
//                    refFreq, target, refFreq, _map.refDelay);

//                spectral_station_response(station, 1, freq.data(), target,
//                    freq.data(), _map.refDelay, &resp);
//                antenna->response(1,target, _map.refDelay, _map.refTile,
//                    freq.data(), freq.data(), &resp);

                resp = antenna->response(refTime.getValue(), refFreq, target,
                    refFreq, _map.refDelay, _map.refTile);

                response(x, y, 0) = resp[0][0];
                response(x, y, 1) = resp[0][1];
                response(x, y, 2) = resp[1][0];
                response(x, y, 3) = resp[1][1];
            }

            cout << response(0, 0, 0) << " " << flush;

            cout << "store..." << flush;

            std::ostringstream oss;
            oss << "response-id-" << *it;
            if(nFrames > 1)
            {
                oss << "-frame-" << j;
            }
            oss << ".img";

            if(abs)
            {
                Cube<Float> ampl(amplitude(response));
                store(coordinates, ampl, oss.str());
            }
            else
            {
                store(coordinates, response, oss.str());
            }

            cout << "done. " << flush;
        }
        cout << endl;

        refTime = refTime + deltaTime;
    }
    cout << " done." << endl;

    return 0;
}
