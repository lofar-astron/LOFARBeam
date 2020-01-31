#include "../AntennaModel.h"
#include "HamakerCoeff.h"

#include <memory>

namespace StationResponse {
class HamakerElementResponse : public AntennaModel
{
public:
    virtual void response(
        double freq,
        double theta,
        double phi,
        std::complex<double> (&response)[2][2]) const final override;

    static std::shared_ptr<HamakerElementResponse> getInstance(const std::string &name);

protected:
    std::string get_path(const char*) const;

    std::unique_ptr<HamakerCoefficients> m_coeffs;
};

class HamakerElementResponseHBA : public HamakerElementResponse
{
public:
    HamakerElementResponseHBA();
};

class HamakerElementResponseLBA : public HamakerElementResponse
{
public:
    HamakerElementResponseLBA();
};
}
