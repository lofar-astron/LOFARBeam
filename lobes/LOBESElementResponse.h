#include "../AntennaModel.h"

#include <memory>

namespace StationResponse {

class LOBESElementResponse : public AntennaModel
{
public:

    LOBESElementResponse(std::string name);

    virtual void response(
        double freq,
        double theta,
        double phi,
        std::complex<double> (&response)[2][2]) const final override;

    static std::shared_ptr<LOBESElementResponse> getInstance(std::string name);
};

} // namespace StationResponse
