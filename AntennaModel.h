#ifndef ANTENNAMODEL_H
#define ANTENNAMODEL_H

#include <complex>
#include <ostream>

#include "MutablePtr.h"

namespace StationResponse
{

enum AntennaModelID {
    Unknown,
    Hamaker,
    LOBES,
    OSKAR
};

std::ostream& operator<<(std::ostream& os, AntennaModelID model);

class AntennaModel
{
public:

    typedef MutablePtr<AntennaModel> Ptr;

    virtual void response(
        double freq,
        double theta,
        double phi,
        std::complex<double> (&result)[2][2]) const = 0;

    virtual void response(
        int    element_id,
        double freq,
        double theta,
        double phi,
        std::complex<double> (&result)[2][2]) const
    {
        response(freq, theta, phi, result);
    }

};

} // namespace StationResponse

#endif
