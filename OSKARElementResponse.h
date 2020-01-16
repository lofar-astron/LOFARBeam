#include "ElementResponse.h"

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
