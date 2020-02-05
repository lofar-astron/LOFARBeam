#include "../ElementResponse.h"
#include "../Singleton.h"

#include "OskarSphericalWaveCoeff.h"

#include <memory>

namespace LOFAR {
namespace StationResponse {

class OSKARElementResponseDipole : public ElementResponse
{
public:
    virtual void response(
        double freq,
        double theta,
        double phi,
        std::complex<double> (&response)[2][2]) const final override;

    static std::shared_ptr<OSKARElementResponseDipole> getInstance()
    {
        return Singleton<OSKARElementResponseDipole>::getInstance();
    }
};

class OSKARElementResponseSphericalWave : public ElementResponse
{
public:

    static std::shared_ptr<OSKARElementResponseSphericalWave> getInstance()
    {
        return Singleton<OSKARElementResponseSphericalWave>::getInstance();
    }

    OSKARElementResponseSphericalWave();

    virtual void response(
        double freq,
        double theta,
        double phi,
        std::complex<double> (&response)[2][2]) const final override;


protected:
    std::string get_path(const char*) const;

    std::unique_ptr<OskarSphericalWaveCoefficients> m_coeffs;

    mutable std::mutex m_mutex;
};

} // namespace StationResponse
} // namespace LOFAR
