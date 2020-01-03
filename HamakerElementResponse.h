#include "ElementResponse.h"
#include "HamakerCoeff.h"

#include <memory>

class HamakerElementResponse : public ElementResponse
{
public:
    virtual void element_response(
        double freq,
        double theta,
        double phi,
        std::complex<double> (&response)[2][2]) const final override;

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
