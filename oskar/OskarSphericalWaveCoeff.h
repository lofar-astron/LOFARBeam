#ifndef OSKAR_SPHERICAL_WAVE_COEFF_H
#define OSKAR_SPHERICAL_WAVE_COEFF_H

#include <iostream>
#include <string>
#include <complex>
#include <cassert>
#include <vector>
#include <cstring>

#include <H5Cpp.h>

class OskarSphericalWaveCoefficients {
    public:
        // Constructor for reading coeff from file
        OskarSphericalWaveCoefficients(
            std::string& filename);

        // Get
        size_t get_l_max() const;

        std::complex<double> get_alpha();

        // HDF5 I/O
        void read_coeffs(
            std::string& filename);

        // Debugging
        void print_coeffs();

    private:
        // Parameters
        double m_l_max;
        double m_m_max;

        // Data
        std::vector<std::complex<double>> m_coeff;
};

#endif
