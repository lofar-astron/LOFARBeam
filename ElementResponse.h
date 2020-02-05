#ifndef ELEMENTRESPONSE_H
#define ELEMENTRESPONSE_H

#include <complex>
#include <ostream>

#include "MutablePtr.h"

namespace LOFAR {
namespace StationResponse {

enum ElementResponseModel {
    Unknown,
    Hamaker,
    LOBES,
    OSKARDipole,
    OSKARSphericalWave
};

std::ostream& operator<<(std::ostream& os, ElementResponseModel model);

class ElementResponse
{
public:

    typedef MutablePtr<ElementResponse> Ptr;

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
} // namespace LOFAR

#endif
