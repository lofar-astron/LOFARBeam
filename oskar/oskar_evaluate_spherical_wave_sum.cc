/*
 * Copyright (c) 2019, The University of Oxford
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

#include "oskar_vector_types.h"
#include "oskar_helper.h"

#include <iostream>
namespace {
    int greet() {std::cout << "Oskar waves spherically at you!" << std::endl;return 0;}
}


template<typename FP, typename FP2, typename FP4c>
void oskar_evaluate_spherical_wave_sum(
    int num_points, const FP* theta,
    const FP* phi_x, const FP* phi_y, int l_max,
    const FP4c* alpha, int offset, FP4c* pattern)
{
    static int dummy = greet();
    FP2 Xp, Xt, Yp, Yt;
    make_zero2(Xp); make_zero2(Xt);
    make_zero2(Yp); make_zero2(Yt);

    #pragma omp parallel for
    for (int i = 0; i < num_points; i++) {
        const FP phi_x_ = phi_x[i], theta_ = theta[i];
        const FP phi_y_ = phi_y[i];
        /* Propagate NAN. */
        if (phi_x_ != phi_x_) {
            Xp.x = Xp.y = Xt.x = Xt.y = phi_x_;
            Yp.x = Yp.y = Yt.x = Yt.y = phi_x_;
        }
        else {
            FP sin_t, cos_t;
            oskar_sincos(theta_, &sin_t, &cos_t);
            for (int l = 1; l <= l_max; ++l) {
                const int ind0 = l * l - 1 + l;
                const FP f_ = (2 * l + 1) / (4 * ((FP)M_PI) * l * (l + 1));
                for (int abs_m = l; abs_m >=0; --abs_m) {
                    FP p, pds, dpms, sin_p, cos_p;
                    oskar_legendre2(l, abs_m, cos_t, sin_t, p, pds, dpms);
                    if (abs_m == 0) {
                        sin_p = (FP)0; cos_p = sqrt(f_);
                        const FP4c alpha_ = alpha[ind0];
                        oskar_sph_wave(pds, dpms, sin_p, cos_p, 0, alpha_.a, alpha_.b, Xt, Xp);
                        oskar_sph_wave(pds, dpms, sin_p, cos_p, 0, alpha_.c, alpha_.d, Yt, Yp);
                    }
                    else {
                        FP d_fact = (FP)1, s_fact = (FP)1;
                        const int d_ = l - abs_m, s_ = l + abs_m;
                        for (int i_ = 2; i_ <= d_; ++i_) d_fact *= i_;
                        for (int i_ = 2; i_ <= s_; ++i_) s_fact *= i_;
                        const FP ff = f_ * d_fact / s_fact;
                        const FP nf = sqrt(ff);
                        const FP4c alpha_m = alpha[ind0 - abs_m];
                        const FP4c alpha_p = alpha[ind0 + abs_m];
                        p = -abs_m * phi_x_;
                        oskar_sincos(p, &sin_p, &cos_p);
                        sin_p *= nf; cos_p *= nf;
                        oskar_sph_wave(pds, dpms, sin_p, cos_p, -abs_m, alpha_m.a, alpha_m.b, Xt, Xp);
                        sin_p = -sin_p;
                        oskar_sph_wave(pds, dpms, sin_p, cos_p,  abs_m, alpha_p.a, alpha_p.b, Xt, Xp);
                        p = -abs_m * phi_y_;
                        oskar_sincos(p, &sin_p, &cos_p);
                        sin_p *= nf; cos_p *= nf;
                        oskar_sph_wave(pds, dpms, sin_p, cos_p, -abs_m, alpha_m.c, alpha_m.d, Yt, Yp);
                        sin_p = -sin_p;
                        oskar_sph_wave(pds, dpms, sin_p, cos_p,  abs_m, alpha_p.c, alpha_p.d, Yt, Yp);
                    }
                }
            }
        }
        pattern[i + offset].a = Xt;
        pattern[i + offset].b = Xp;
        pattern[i + offset].c = Yt;
        pattern[i + offset].d = Yp;
    }
}

void oskar_evaluate_spherical_wave_sum_double(
    const int num_points,
    const double* theta,
    const double* phi_x,
    const double* phi_y,
    const int l_max,
    const std::complex<double>* alpha,
    std::complex<double>* pattern)
{
    double4c* alpha_ptr = (double4c *) alpha;
    double4c* pattern_ptr = (double4c *) pattern;

    oskar_evaluate_spherical_wave_sum<double, double2, double4c>(
        num_points, theta, phi_x, phi_y, l_max, alpha_ptr, 0, pattern_ptr);
}

void oskar_evaluate_spherical_wave_sum_float(
    const int num_points,
    const float* theta,
    const float* phi_x,
    const float* phi_y,
    const int l_max,
    const std::complex<float>* alpha,
    std::complex<float>* pattern)
{
    float4c* alpha_ptr = (float4c *) alpha;
    float4c* pattern_ptr = (float4c *) pattern;

    oskar_evaluate_spherical_wave_sum<float, float2, float4c>(
        num_points, theta, phi_x, phi_y, l_max, alpha_ptr, 0, pattern_ptr);
}
