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

/* OUT_S += A * B */
template<typename FP>
inline void oskar_mul_add_complex(FP& out_s, FP A, FP B)
{
    out_s.x += A.x * B.x; out_s.x -= A.y * B.y;
    out_s.y += A.x * B.y; out_s.y += A.y * B.x;
}

/* OUT_S -= A * B */
template<typename FP>
inline void oskar_mul_sub_complex(FP& out_s, FP A, FP B)
{
    out_s.x -= A.x * B.x; out_s.x += A.y * B.y;
    out_s.y -= A.x * B.y; out_s.y -= A.y * B.x;
}

template<typename FP>
inline void make_zero2(FP& X)
{
    X.x = X.y = 0;
};

template<typename FP, typename FP2>
inline void oskar_sph_wave(
    FP pds,
    FP dpms,
    FP sin_p,
    FP cos_p,
    int M,
    FP2 a_tm,
    FP2 a_te,
    FP2& c_theta,
    FP2& c_phi)
{
    FP2 qq, dd;
    qq.x = -cos_p * dpms;
    qq.y = -sin_p * dpms;
    dd.x = -sin_p * pds * (M);
    dd.y = cos_p * pds * (M);
    oskar_mul_add_complex(c_phi, qq, a_tm);
    oskar_mul_sub_complex(c_phi, dd, a_te);
    oskar_mul_add_complex(c_theta, dd, a_tm);
    oskar_mul_add_complex(c_theta, qq, a_te);
}

/* OUT0 is P_l^m(cos_theta),
 * OUT1 is P_l^m(cos_theta) / sin_theta,
 * OUT2 is d/d(cos_theta){P_l^m(cos_theta)} * sin_theta. */
template<typename FP>
inline void oskar_legendre2(
    int l,
    int m,
    FP cos_theta,
    FP sin_theta,
    FP& out0,
    FP& out1,
    FP& out2)
{
    FP p0_ = (FP) 1, p1_;
    if (m > 0) {
        FP fact = (FP) 1;
        for (int i_ = 1; i_ <= m; ++i_) {
            p0_ *= (-fact) * sin_theta;
            fact += (FP) 2;
        }
    }
    out0 = cos_theta * (2 * m + 1) * p0_;
    if (l == m) {
        p1_ = out0;
        out0 = p0_;
    }
    else {
        p1_ = out0;
        for (int i_ = m + 2; i_ <= l + 1; ++i_) {
            out0 = p1_;
            p1_ = ((2 * i_ - 1) * cos_theta * out0 - (i_ + m - 1) * p0_) / (i_ - m);
            p0_ = out0;
        }
    }
    if (sin_theta != (FP) 0) {
        /* BOTH of these are divides. */
        out1 = out0 / sin_theta;
        out2 = (cos_theta * out0 * (l + 1) - p1_ * (l - m + 1)) / sin_theta;
    }
    else {
        out1 = out2 = (FP) 0;
    }
}

inline void oskar_sincos(float x, float* s, float* c) { sincosf(x, s, c); }
inline void oskar_sincos(double x, double* s, double* c) { sincos(x, s, c); }
