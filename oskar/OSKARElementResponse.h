#include "../ElementResponse.h"
#include "../Singleton.h"

#include <memory>

namespace LOFAR {
namespace StationResponse {

class OSKARElementResponse : public ElementResponse
{

};

class OSKARElementResponseDipole : public OSKARElementResponse
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

} // namespace StationResponse
} // namespace LOFAR
