/*
 * Copyright (c) 2014-2019, The University of Oxford
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the University of Oxford nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <cmath>
#include <complex>

#include "oskar.h"

#include <iostream>
namespace {
    int greet() {std::cout << "Greetings from Oskar Dipole!" << std::endl;return 0;}
}

template<typename FP, typename FP2>
inline void oskar_dipole(
    FP theta,
    FP phi,
    FP kL,
    FP2& e_theta,
    FP2& e_phi)
{
    FP sin_phi, cos_phi;
    FP sin_theta, cos_theta;
    oskar_sincos(phi, &sin_phi, &cos_phi);
    oskar_sincos(theta, &sin_theta, &cos_theta);
    FP cos_kL = (FP) cos(kL);
    const FP denom = (FP)1 + cos_phi*cos_phi * (cos_theta*cos_theta - (FP)1);
    if (denom == (FP)0)
        e_theta.x = e_theta.y = e_phi.x = e_phi.y = (FP)0;
    else {
        const FP q = kL * cos_phi * sin_theta;
        const FP cos_q = cos(q);
        const FP t = (cos_q - cos_kL) / denom;
        e_theta.x = -cos_phi * cos_theta * t;
        e_phi.x = sin_phi * t;
        e_phi.y = e_theta.y = (FP)0;
    }
}

template<typename FP, typename FP2>
void oskar_evaluate_dipole_pattern(
    const int num_points,
    const FP* theta,
    const FP* phi,
    const FP kL,
    const int stride,
    const int E_theta_offset,
    const int E_phi_offset,
    FP2* e_theta,
    FP2* e_phi)
{
    static int dummy = greet();
    for (int i = 0; i < num_points; i++) {
        const int i_out = i * stride;\
        const int theta_out = i_out + E_theta_offset;\
        const int phi_out   = i_out + E_phi_offset;\
        oskar_dipole<FP, FP2>(theta[i], phi[i], kL, e_theta[theta_out], e_phi[phi_out]);
    }
}

void oskar_evaluate_dipole_pattern_double(
    const int num_points,
    const double* theta,
    const double* phi,
    const double freq_hz,
    const double dipole_length_m,
    std::complex<double>* pattern)
{
    const int stride = 4;
    const int offset = 0;
    const int E_theta_offset = offset;
    const int E_phi_offset = offset + 1;
    const double kL = dipole_length_m * (M_PI * freq_hz / 299792458);
    double2* pattern_ptr = (double2 *) pattern;

    oskar_evaluate_dipole_pattern<double, double2>(
        num_points, theta, phi, kL, stride, E_theta_offset, E_phi_offset, pattern_ptr, pattern_ptr);
}

void oskar_evaluate_dipole_pattern_float(
    const int num_points,
    const float* theta,
    const float* phi,
    const float freq_hz,
    const float dipole_length_m,
    std::complex<double>* pattern)
{
    const int stride = 4;
    const int offset = 0;
    const int E_theta_offset = offset;
    const int E_phi_offset = offset + 1;
    const float kL = dipole_length_m * (M_PI * freq_hz / 299792458);
    float2* pattern_ptr = (float2 *) pattern;

    oskar_evaluate_dipole_pattern<float, float2>(
        num_points, theta, phi, kL, stride, E_theta_offset, E_phi_offset, pattern_ptr, pattern_ptr);
}
