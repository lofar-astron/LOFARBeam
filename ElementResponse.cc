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

#include <cmath>

// The coefficients are kept in an unnamed namespace which effectively makes
// them invisible outside this translation unit.
namespace
{
// PI / 2.0
const double pi_2 = 1.570796326794896619231322;

#include "DefaultCoeffLBA.cc"
#include "DefaultCoeffHBA.cc"
}

namespace LOFAR
{

void element_response_lba(double freq, double theta, double phi,
    std::complex<double> (&response)[2][2])
{
    element_response(freq, theta, phi, response, default_lba_freq_center,
        default_lba_freq_range, default_lba_coeff_shape, default_lba_coeff);
}

void element_response_hba(double freq, double theta, double phi,
    std::complex<double> (&response)[2][2])
{
    element_response(freq, theta, phi, response, default_hba_freq_center,
        default_hba_freq_range, default_hba_coeff_shape, default_hba_coeff);
}

void element_response(double freq, double theta, double phi,
    std::complex<double> (&response)[2][2], double freq_center,
    double freq_range, const unsigned int (&coeff_shape)[3],
    const std::complex<double> coeff[])
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

    const unsigned int nHarmonics  = coeff_shape[0];
    const unsigned int nPowerTheta = coeff_shape[1];
    const unsigned int nPowerFreq  = coeff_shape[2];
    const unsigned int nInner      = 2;

    // Define an multi-dimensional array to index the coefficients
    typedef std::complex<double> coeff_type[nHarmonics][nPowerTheta][nPowerFreq][nInner];
    const coeff_type *coeff_arr = (coeff_type *) coeff;

    // The model is parameterized in terms of a normalized frequency in the
    // range [-1, 1]. The appropriate conversion is taken care of below.
    freq = (freq - freq_center) / freq_range;

    // The variables sign and kappa are used to compute the value of kappa
    // mentioned in the description of the beam model [kappa = (-1)^k * (2 * k
    //+ 1)] incrementally.
    int sign = 1, kappa = 1;

    std::complex<double> P[2], Pj[2];
    for(unsigned int k = 0; k < nHarmonics; ++k)
    {
        // Compute the (diagonal) projection matrix P for the current harmonic.
        // This requires the evaluation of two polynomials in theta and freq (of
        // degree nPowerTheta in theta and nPowerFreq in freq), one for each
        // element of P. The polynomials are evaluated using Horner's rule.

        // Horner's rule requires backward iteration of the coefficients, so
        // start indexing the block of coefficients at the last element

        // Evaluate the highest order term.
        P[0] = (*coeff_arr)[k][nPowerTheta-1][nPowerFreq-1][0];
        P[1] = (*coeff_arr)[k][nPowerTheta-1][nPowerFreq-1][1];

        for(unsigned int i = 0; i < nPowerFreq - 1; ++i)
        {
            P[0] = P[0] * freq + (*coeff_arr)[k][nPowerTheta-1][nPowerFreq-i-2][0];
            P[1] = P[1] * freq + (*coeff_arr)[k][nPowerTheta-1][nPowerFreq-i-2][1];
        }

        // Evaluate the remaining terms.
        for(unsigned int j = 0; j < nPowerTheta - 1; ++j)
        {
            Pj[0] = (*coeff_arr)[k][nPowerTheta-j-2][nPowerFreq-1][0];
            Pj[1] = (*coeff_arr)[k][nPowerTheta-j-2][nPowerFreq-1][1];

            for(unsigned int i = 0; i < nPowerFreq - 1; ++i)
            {
                Pj[0] = Pj[0] * freq + (*coeff_arr)[k][nPowerTheta-j-2][nPowerFreq-i-2][0];
                Pj[1] = Pj[1] * freq + (*coeff_arr)[k][nPowerTheta-j-2][nPowerFreq-i-2][1];
            }

            P[0] = P[0] * theta + Pj[0];
            P[1] = P[1] * theta + Pj[1];
        }

        // Compute the Jones matrix for the current harmonic, by rotating P over
        // kappa * az, and add it to the result.
        const double angle = sign * kappa * phi;
        const double caz = std::cos(angle);
        const double saz = std::sin(angle);

        response[0][0] += caz * P[0];
        response[0][1] += -saz * P[1];
        response[1][0] += saz * P[0];
        response[1][1] += caz * P[1];

        // Update sign and kappa.
        sign = -sign;
        kappa += 2;
    }
}

} //# namespace LOFAR
