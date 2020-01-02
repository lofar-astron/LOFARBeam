//# ElementResponse.cc: Functions to compute the (idealized) response of a LOFAR
//# LBA or HBA dual dipole antenna.
//#
//# Copyright (C) 2011
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$

#include "ElementResponse.h"

#include "config.h"

#include <cmath>
#include <sstream>

// PI / 2.0
const double pi_2 = 1.570796326794896619231322;

namespace LOFAR
{

void element_response_lba(double freq, double theta, double phi,
    std::complex<double> (&response)[2][2],
    HamakerCoefficients& lba_coeff)
{
    element_response(freq, theta, phi, response, lba_coeff);
}

void element_response_hba(double freq, double theta, double phi,
    std::complex<double> (&response)[2][2],
    HamakerCoefficients& hba_coeff)
{
    element_response(freq, theta, phi, response, hba_coeff);
}

void element_response(double freq, double theta, double phi,
    std::complex<double> (&response)[2][2],
    HamakerCoefficients& coeffs)
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


    const double freq_center       = coeffs.get_freq_center();
    const double freq_range        = coeffs.get_freq_range();
    const unsigned int nHarmonics  = coeffs.get_nHarmonics();
    const unsigned int nPowerTheta = coeffs.get_nPowerTheta();
    const unsigned int nPowerFreq  = coeffs.get_nPowerFreq();;

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
        P = coeffs.get_coeff(k, nPowerTheta-1, nPowerFreq-1);

        for(unsigned int i = 0; i < nPowerFreq - 1; ++i)
        {
            auto Pk = coeffs.get_coeff(k, nPowerTheta-1, nPowerFreq-i-2);
            P.first  = P.first  * freq + Pk.first;
            P.second = P.second * freq + Pk.second;
        }

        // Evaluate the remaining terms.
        for(unsigned int j = 0; j < nPowerTheta - 1; ++j)
        {
            Pj = coeffs.get_coeff(k, nPowerTheta-j-2, nPowerFreq-1);

            for(unsigned int i = 0; i < nPowerFreq - 1; ++i)
            {
                auto Pk = coeffs.get_coeff(k, nPowerTheta-j-2, nPowerFreq-i-2);
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

} //# namespace LOFAR
