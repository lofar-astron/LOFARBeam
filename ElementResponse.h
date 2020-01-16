#ifndef ELEMENT_RESPONSE_H
#define ELEMENT_RESPONSE_H

#include <complex>

class ElementResponse
{
public:
    virtual void element_response(
        double freq,
        double theta,
        double phi,
        std::complex<double> (&response)[2][2]) const = 0;
};

#endif
