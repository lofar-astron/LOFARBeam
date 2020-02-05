//# Dir2ITRF.cc: Convertor that maps time to an ITRF direction.

#include "ITRFDirection.h"
#include "ITRFConverter.h"

#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>

namespace LOFAR {
namespace StationResponse {

  //TODO: Initialize converter with a time (and fixed position) and convert specific directions.
  //      Needed for wslean as well as for the makestationresponse executable.


ITRFConverter::ITRFConverter(real_t time)
{
    //create ITRF Direction from fixed stationposition
    casacore::MVPosition mvPosition(ITRFDirection::LOFARPosition()[0], ITRFDirection::LOFARPosition()[1], ITRFDirection::LOFARPosition()[2]);
    casacore::MPosition mPosition(mvPosition, casacore::MPosition::ITRF);
    casacore::MEpoch timeEpoch(casacore::Quantity(time, "s"));
    itsFrame = casacore::MeasFrame(timeEpoch, mPosition);

    // Order of angles seems to be longitude (along the equator), lattitude
    // (towards the pole).
    itsConverter = casacore::MDirection::Convert(
        casacore::MDirection::J2000,
        casacore::MDirection::Ref(casacore::MDirection::ITRF, itsFrame));
}

void ITRFConverter::setTime(real_t time)
{
    // Cannot use MeasFrame::resetEpoch(Double), because that assumes the
    // argument is UTC in (fractional) days (MJD).
    itsFrame.resetEpoch(casacore::Quantity(time, "s"));
}

vector3r_t ITRFConverter::j2000ToITRF(const vector2r_t &j2000Direction) const
{
    casacore::MVDirection mvDirection(j2000Direction[0], j2000Direction[1]);
    casacore::MDirection mDirection(mvDirection, casacore::MDirection::J2000);
    const casacore::MVDirection mvITRF = itsConverter(mDirection).getValue();

    return vector3r_t{{mvITRF(0), mvITRF(1), mvITRF(2)}};
}

vector3r_t ITRFConverter::j2000ToITRF(const vector3r_t &j2000Direction) const
{
    casacore::MVDirection mvDirection(j2000Direction[0], j2000Direction[1], j2000Direction[2]);
    casacore::MDirection mDirection(mvDirection, casacore::MDirection::J2000);

    const casacore::MVDirection mvITRF = itsConverter(mDirection).getValue();

    return vector3r_t{{mvITRF(0), mvITRF(1), mvITRF(2)}};
}

vector3r_t ITRFConverter::toITRF(const casacore::MDirection &direction) const{
    const casacore::MVDirection mvITRF = itsConverter(direction).getValue();
    return vector3r_t{{mvITRF(0), mvITRF(1), mvITRF(2)}};
}

casacore::MDirection ITRFConverter::toDirection(const vector2r_t &j2000Direction) const
{
    casacore::MVDirection mvDirection(j2000Direction[0], j2000Direction[1]);
    casacore::MDirection mDirection(mvDirection, casacore::MDirection::J2000);
    return itsConverter(mDirection);
}

casacore::MDirection ITRFConverter::toDirection(const vector3r_t &j2000Direction) const
{
    casacore::MVDirection mvDirection(j2000Direction[0], j2000Direction[1], j2000Direction[2]);
    casacore::MDirection mDirection(mvDirection, casacore::MDirection::J2000);

    return itsConverter(mDirection);
}

casacore::MDirection ITRFConverter::toDirection(const casacore::MDirection &direction) const
{
    return itsConverter(direction);
}


} //# namespace StationResponse
} // namespace LOFAR
