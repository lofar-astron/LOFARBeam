#include "../AntennaModel.h"

#include <memory>

namespace StationResponse {

class OSKARElementResponse : public AntennaModel
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
};

}
