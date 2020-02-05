//# HamakerElementResponse.cc:
//# Functions to compute the (idealized) response of a LOFAR
//# LBA or HBA dual dipole antenna.


#include<stdexcept>

#include "config.h"

#include "HamakerElementResponse.h"
#include "../Singleton.h"

namespace LOFAR {
namespace StationResponse {

// PI / 2.0
const double pi_2 = 1.570796326794896619231322;

std::shared_ptr<HamakerElementResponse> HamakerElementResponse::getInstance(const std::string &name)
{
    if (name.length() >= 3) {
        if (name.substr(name.length()-3, 3) == "LBA") {
            return Singleton<HamakerElementResponseLBA>::getInstance();
        }
        if (name.substr(name.length()-3, 3) == "HBA") {
            return Singleton<HamakerElementResponseHBA>::getInstance();
        }
    }
    throw std::invalid_argument("HamakerElementResponse::getInstance: name should end in either 'LBA' or 'HBA'");
}


std::string HamakerElementResponse::get_path(
    const char* filename) const
{
    std::stringstream ss;
    ss << LOFARBEAM_DATA_DIR << "/";
    ss << filename;
    return ss.str();
}

void HamakerElementResponse::response(
    double freq,
    double theta,
    double phi,
    std::complex<double> (&response)[2][2]) const
{
    // Initialize the response to zero.
    response[0][0] = 0.0;
    response[0][1] = 0.0;
    response[1][0] = 0.0;
    response[1][1] = 0.0;

    // Clip directions below the horizon.
    if(theta >= pi_2)
    {
        return;
    }

    const double freq_center       = m_coeffs->get_freq_center();
    const double freq_range        = m_coeffs->get_freq_range();
    const unsigned int nHarmonics  = m_coeffs->get_nHarmonics();
    const unsigned int nPowerTheta = m_coeffs->get_nPowerTheta();
    const unsigned int nPowerFreq  = m_coeffs->get_nPowerFreq();;

    // The model is parameterized in terms of a normalized frequency in the
    // range [-1, 1]. The appropriate conversion is taken care of below.
    freq = (freq - freq_center) / freq_range;

    // The variables sign and kappa are used to compute the value of kappa
    // mentioned in the description of the beam model [kappa = (-1)^k * (2 * k
    //+ 1)] incrementally.
    int sign = 1, kappa = 1;

    std::pair<std::complex<double>, std::complex<double>> P;
    std::pair<std::complex<double>, std::complex<double>> Pj;
    for(unsigned int k = 0; k < nHarmonics; ++k)
    {
        // Compute the (diagonal) projection matrix P for the current harmonic.
        // This requires the evaluation of two polynomials in theta and freq (of
        // degree nPowerTheta in theta and nPowerFreq in freq), one for each
        // element of P. The polynomials are evaluated using Horner's rule.

        // Horner's rule requires backward iteration of the coefficients, so
        // start indexing the block of coefficients at the last element

        // Evaluate the highest order term.
        P = m_coeffs->get_coeff(k, nPowerTheta-1, nPowerFreq-1);

        for(unsigned int i = 0; i < nPowerFreq - 1; ++i)
        {
            auto Pk = m_coeffs->get_coeff(k, nPowerTheta-1, nPowerFreq-i-2);
            P.first  = P.first  * freq + Pk.first;
            P.second = P.second * freq + Pk.second;
        }

        // Evaluate the remaining terms.
        for(unsigned int j = 0; j < nPowerTheta - 1; ++j)
        {
            Pj = m_coeffs->get_coeff(k, nPowerTheta-j-2, nPowerFreq-1);

            for(unsigned int i = 0; i < nPowerFreq - 1; ++i)
            {
                auto Pk = m_coeffs->get_coeff(k, nPowerTheta-j-2, nPowerFreq-i-2);
                Pj.first  = Pj.first  * freq + Pk.first;
                Pj.second = Pj.second * freq + Pk.second;
            }

            P.first  = P.first  * theta + Pj.first;
            P.second = P.second * theta + Pj.second;
        }

        // Compute the Jones matrix for the current harmonic, by rotating P over
        // kappa * az, and add it to the result.
        const double angle = sign * kappa * phi;
        const double caz = std::cos(angle);
        const double saz = std::sin(angle);

        response[0][0] += caz * P.first;
        response[0][1] += -saz * P.second;
        response[1][0] += saz * P.first;
        response[1][1] += caz * P.second;

        // Update sign and kappa.
        sign = -sign;
        kappa += 2;
    }
}

HamakerElementResponseHBA::HamakerElementResponseHBA()
{
    std::string path = get_path("HamakerHBACoeff.h5");
    m_coeffs.reset(new HamakerCoefficients(path));
}

HamakerElementResponseLBA::HamakerElementResponseLBA()
{
    std::string path = get_path("HamakerLBACoeff.h5");
    m_coeffs.reset(new HamakerCoefficients(path));
}

} // namespace StationResponse
} // namespace LOFAR
