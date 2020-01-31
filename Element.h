#ifndef ELEMENT_H
#define ELEMENT_H

#include <complex>

class Element
{
public:

    virtual void response(
        double time,
        double freq,
        double theta,
        double phi,
        std::complex<double> (&result)[2][2]) const = 0;

    double m_position[3];

};

#endif
