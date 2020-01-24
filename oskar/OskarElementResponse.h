#include "../ElementResponse.h"
#include "OskarSphericalWaveCoeff.h"

#include <memory>

class OSKARElementResponse : public ElementResponse
{
};

class OSKARElementResponseDipole : public OSKARElementResponse
{
public:
    virtual void element_response(
        double freq,
        double theta,
        double phi,
        std::complex<double> (&response)[2][2]) const final override;
};

class OskarElementResponseSphericalWave : public OSKARElementResponse
{
public:
    OskarElementResponseSphericalWave();

    virtual void element_response(
        double freq,
        double theta,
        double phi,
        std::complex<double> (&response)[2][2]) const final override;

protected:
    std::string get_path(const char*) const;

    std::unique_ptr<OskarSphericalWaveCoefficients> m_coeffs;

    mutable std::mutex m_mutex;
};