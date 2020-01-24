#ifndef OSKAR_SPHERICAL_WAVE_COEFF_H
#define OSKAR_SPHERICAL_WAVE_COEFF_H

#include <iostream>
#include <string>
#include <complex>
#include <cassert>
#include <vector>
#include <cstring>
#include <memory>
#include <mutex>

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

        bool read_frequency(
            const unsigned int freq);

        // Debugging
        void print_coeffs();

    private:
        // Parameters
        double m_l_max;
        double m_m_max;
        unsigned int m_freq = 0;
        bool m_dataset_available = false;
        std::unique_ptr<H5::H5File> m_h5_file;
        std::unique_ptr<H5::DataSet> m_h5_dataset;

        // Data
        std::vector<std::complex<double>> m_coeff;
};

#endif
