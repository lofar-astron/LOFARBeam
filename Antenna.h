#ifndef ANTENNA_H
#define ANTENNA_H

#include <complex>

class Antenna : public Element
{
public:

    Antenna(int id, AntennaModel::Ptr model) :
        m_id(id),
        m_model(model)
    {}

    virtual void response(
        double time,
        double freq,
        double theta,
        double phi,
        std::complex<double> (&result)[2][2]) override
    {
        m_antenna_model->response(freq, theta, phi, std::complex<double> (&result)[2][2]);
    }

private:
    int m_antenna_id;
    AntennaModel::Ptr m_antenna_model;

};

#endif
