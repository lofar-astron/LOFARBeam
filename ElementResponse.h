//# ElementResponse.h: Functions to compute the (idealized) response of a LOFAR
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

#ifndef LOFAR_ELEMENTRESPONSE_H
#define LOFAR_ELEMENTRESPONSE_H

// \file
// Functions to compute the (idealized) response of a LOFAR LBA or HBA dual
// dipole antenna.

#include <complex>

namespace LOFAR
{

// \addtogroup ElementResponse
// @{

// Compute the response of an idealized LOFAR LBA dual dipole antenna to
// radiation at frequency freq (Hz) arriving from the direction given by theta,
// phi (rad). The antenna model is described in a spherical coordinate system
// with coordinates theta (zenith angle) and phi (azimith). The +X dipole is at
// azimuth zero, the +Y dipole is at azimuth PI / 2.0.
//
// Preconditions:
// --------------
// freq: Frequency in Hz in the range [10 MHz, 100 MHz].
// theta: Zenith angle in rad in the range [0.0, PI / 2.0].
// phi: Azimuth in rad in the range [0.0, 2.0 * PI].
//
void element_response_lba(double freq, double theta, double phi,
    std::complex<double> (&response)[2][2]);

// Compute the response of an idealized LOFAR HBA dual dipole antenna to
// radiation at frequency freq (Hz) arriving from the direction given by theta,
// phi (rad). The antenna model is described in a spherical coordinate system
// with coordinates theta (zenith angle) and phi (azimith). The +X dipole is at
// azimuth zero, the +Y dipole is at azimuth PI / 2.0.
//
// Preconditions:
// --------------
// freq: Frequency in Hz in the range [120 MHz, 240 MHz].
// theta: Zenith angle in rad in the range [0.0, PI / 2.0].
// phi: Azimuth in rad in the range [0.0, 2.0 * PI].
//
void element_response_hba(double freq, double theta, double phi,
    std::complex<double> (&response)[2][2]);

// Compute the response of an idealized LOFAR dual dipole antenna to radiation
// at frequency freq (Hz) arriving from the direction given by theta, phi (rad).
// The antenna model is described in a spherical coordinate system with
// coordinates theta (zenith angle) and phi (azimith). The +X dipole is at
// azimuth zero, the +Y dipole is at azimuth PI / 2.0.
//
// This function uses a set of user defined coefficients to evaluate the beam
// model. The coeff_shape parameter defines the shape of the coefficient array
// as no. of harmonics x degree in theta x degree in frequency x 2. The last
// dimension is implicit and always equal to 2. The coeff parameter points to an
// array of coefficients of the proper size, stored in row-major order
// ("C"-order). The freq_center and freq_range parameters define the frequency
// range covered by the model described by the set of coefficients.
//
// Preconditions:
// --------------
// freq: Frequency in Hz in the range [freq_center - freq_range,
//     freq_center + freq_range].
// theta: Zenith angle in rad in the range [0.0, PI / 2.0].
// phi: Azimuth in rad in the range [0.0, 2.0 * PI].
// freq_range, freq_center: Frequency center and range in Hz, should be > 0.
// coeff_shape: Shape of the coefficient array, all dimensions should be > 0.
//
void element_response(double freq, double theta, double phi,
    std::complex<double> (&response)[2][2], double freq_center,
    double freq_range, const unsigned int (&coeff_shape)[3],
    const std::complex<double> coeff[]);

// @}

} //# namespace LOFAR

#endif
